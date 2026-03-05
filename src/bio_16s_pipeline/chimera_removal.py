import argparse
import time
import multiprocessing
from pathlib import Path
from collections import Counter
from .db_utils import open_kmer_db, fetch_all_kmer_labels
from .common_func import read_fasta, save_files, yield_kmers
from .logs import start_log_listener, stop_log_listener, init_worker_logger, get_logger


logger = get_logger()


### ====== base functions ====== ###
def build_kmer_cache(seqs, k):
    seq_kmers = {}
    kmer_counter = Counter()
    for header, seq in seqs:
        kmers = list(yield_kmers(seq, k))
        seq_kmers[(header, seq)] = kmers
        kmer_counter.update(kmers)
    return seq_kmers, kmer_counter


def dominant_label(ref_kmer, kmers, kmer_counter):
    label_counter = Counter()
    for km in kmers:
        labels = ref_kmer.get(km, ())
        weight = 1 if kmer_counter[km] <= 10 else 1 / kmer_counter[km] ** 0.5
        for l in labels:
            label_counter[l] += weight
    if not label_counter:
        return None, 0
    return label_counter.most_common(1)[0]


def chimera_check(kmers, ref_kmer, kmer_counter, min_match_thresh):
    mid = len(kmers) // 2
    if mid < 1 or len(kmers) - mid < 1:
        return False, None, None
    l_kmers, r_kmers = kmers[:mid], kmers[mid:]
    l_label, l_cnt = dominant_label(ref_kmer, l_kmers, kmer_counter)
    r_label, r_cnt = dominant_label(ref_kmer, r_kmers, kmer_counter)
    l_ratio = l_cnt / len(l_kmers) if len(l_kmers) > 0 else 0
    r_ratio = r_cnt / len(r_kmers) if len(r_kmers) > 0 else 0
    is_chimera = (l_label != r_label and l_ratio >= min_match_thresh and r_ratio >= min_match_thresh)
    return is_chimera, l_label, r_label


def process_sample(args):
    fasta_path, db_path, k, min_match_thresh = args
    worker_logger = get_logger()
    
    try:
        seqs = list(read_fasta(fasta_path))
        if not seqs:
            return []

        seq_kmers, kmer_counter = build_kmer_cache(seqs, k)

        with open_kmer_db(db_path) as conn:
            ref_kmer = fetch_all_kmer_labels(conn, k)

        results = []
        for (header, seq), kmers in seq_kmers.items():
            is_chimera, l_lab, r_lab = chimera_check(kmers, ref_kmer, kmer_counter, min_match_thresh)
            results.append((header, seq, is_chimera))
        
        return results

    except Exception as e:
        worker_logger.error(f"Error processing {fasta_path}: {e}")
        return []


def run_pipeline(args):

    ctx = multiprocessing.get_context("spawn")
    q = ctx.Queue()
    start_log_listener(queue=q)

    start_time = time.time()
 
    try:
        base_path = Path(args.input)
        fasta_files = list(base_path.rglob("*_merged.fasta.gz")) + list(base_path.rglob("*_filtered.fasta.gz"))
 
        if not fasta_files:
            logger.error(f"No processed fasta.gz files found in {args.input}")
            return

        logger.info(f"Found {len(fasta_files)} samples for chimera removal.")

        tasks = []
        for f in fasta_files:
            output_name = f.name.replace(".fasta.gz", "_no_chi.fasta.gz")
            output_file = Path(args.output) / f.parent.name / output_name

            if output_file.exists():
                logger.info(f"Skip sample {f.parent.name}, output file already exists.")
                continue

            tasks.append((f, args.db, args.k, args.thresh))

        if not tasks:
            logger.info("All samples are already processed. No action needed.")
            return

        with ctx.Pool(processes=args.threads, initializer=init_worker_logger, initargs=(q,)) as pool:
            all_results = pool.map(process_sample, tasks)

        for i, (f, _, _, _) in enumerate(tasks):
            results = all_results[i]
            if not results: 
                continue

            output_name = f.name.replace(".fasta.gz", "_no_chi.fasta.gz")
            output_file = Path(args.output) / f.parent.name / output_name
            output_file.parent.mkdir(parents=True, exist_ok=True)

            non_chimeras = [(h, s) for h, s, is_chi in results if not is_chi]
            save_files(output_file, non_chimeras)
            logger.info(f"Sample {f.parent.name} done: {len(non_chimeras)}/{len(results)} kept.")

    except Exception as e:
        logger.error(f"Critical error: {e}", exc_info=True)
    finally:
        logger.info(f"Chimera removal finished in {time.time() - start_time:.2f}s.")
        stop_log_listener()


### ====== main function ====== ###
def main():
    parser = argparse.ArgumentParser(description="16S rRNA Pipeline - Chimera Removal", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    group_base = parser.add_argument_group("Global Settings")
    group_base.add_argument("-i", "--input", required=True, help="Root path to search for filtered.fasta or merged.fasta")    
    group_base.add_argument("-o", "--output", required=True, help="Output directory")
    group_base.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads for parallel processing")    

    group_chi = parser.add_argument_group("Chimera Settings")
    group_chi.add_argument("--db", required=True, help="Path to K-mer SQLite database")
    group_chi.add_argument("-k", type=int, default=25, help="K-mer length (must match database)")
    group_chi.add_argument("--thresh", type=float, default=0.2, help="Chimera detection threshold")

    args = parser.parse_args()
    run_pipeline(args)

if __name__ == '__main__':
    main()

