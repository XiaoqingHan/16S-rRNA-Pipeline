import time
import argparse
import multiprocessing
import pandas as pd
from pathlib import Path
from collections import Counter
from .common_func import read_fasta, yield_kmers
from .db_utils import open_kmer_db, fetch_all_kmer_labels
from .logs import get_logger, init_worker_logger, start_log_listener, stop_log_listener


logger = get_logger()


# get taxon of each seq by kmer hits
def classify(seq, k, kmer_map):
    votes = Counter()
    for kmer in yield_kmers(seq, k):
        if kmer in kmer_map:
            votes.update(kmer_map[kmer])
    if not votes:
        return "Unclassified"
    return votes.most_common(1)[0][0]


# calculate abundance
def compute_abundance(seqs, sample_id, k, kmer_map):
    counter = Counter()
    for seq in seqs:
        taxon = classify(seq, k, kmer_map)
        counter[taxon] += 1
    total = sum(counter.values()) or 1
    rows = []
    for taxon, cnt in counter.items():
        # consider formats: [Genus] species, Genus_species, Genus_1
        genus = taxon.strip("[]").replace('_', ' ').split()[0]
        rows.append({
            'Sample_ID': sample_id,
            'Taxon': taxon,
            'Genus': genus,
            'Absolute_Abundance': cnt,
            'Relative_Abundance': (cnt / total) * 100
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values('Absolute_Abundance', ascending=False).reset_index(drop=True)
    return df


def process_sample(args):
    sample_dir, db_path, k, file_suffix, out_path = args
    logger = get_logger()
    sample_id = sample_dir.name

    suffixes = [file_suffix, "_merged_no_chi.fasta.gz", "_filtered_no_chi.fasta.gz"]
    input_file = None
    for suffix in suffixes:
        if not suffix:
            continue
        candidate = sample_dir / f"{sample_id}{suffix}"
        if candidate.exists():
            input_file = candidate
            break

    if input_file is None:
        logger.error(f"[{sample_id}] Required non-chimera file not found in {sample_dir}.")
        return None

    try:
        logger.info(f"Processing abundance for sample: {sample_id} (k={k})")
        seqs = [seq for _, seq in read_fasta(input_file)]
        if not seqs:
            logger.info(f"Sample {sample_id} has no sequences.")
            return sample_id

        with open_kmer_db(db_path) as conn:
            # Note: k is now passed directly instead of best_k
            kmer_map = fetch_all_kmer_labels(conn, k)

        abundance_df = compute_abundance(seqs, sample_id, k, kmer_map)
        # translate species abundance to genus level for each sample
        genus_df = abundance_df.groupby(['Sample_ID', 'Genus']).agg({'Absolute_Abundance': 'sum', 'Relative_Abundance': 'sum'}).reset_index()
        genus_df = genus_df.sort_values('Absolute_Abundance', ascending=False)

        sample_dir = Path(out_path) / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        abundance_df.to_csv(sample_dir / f"{sample_id}_species_abundance.csv", index=False)
        genus_df.to_csv(sample_dir / f"{sample_id}_genus_abundance.csv", index=False)

        return sample_id

    except Exception as e:
        logger.error(f"Error in processing {sample_id}: {e}", exc_info=True)
        return None


def run_pipeline(args):
    ctx = multiprocessing.get_context("spawn")
    q = ctx.Queue()
    start_log_listener(queue=q)
    
    start_time = time.time()
    logger.info("Starting Taxonomy Abundance Analysis.")

    try:
        base_path = Path(args.input)
        all_sequence_files = list(base_path.rglob("*.fasta.gz"))
        sample_dirs = sorted(list(set(f.parent for f in all_sequence_files)))
 
        if not sample_dirs:
            logger.error(f"No sequence files found in {args.input}. Please check the path.")
            return

        tasks = []
        for d in sample_dirs:
            tasks.append((d, args.db, args.k, args.suffix, args.output))

        with ctx.Pool(processes=args.threads, 
                      initializer=init_worker_logger, 
                      initargs=(q,)) as pool:
            pool.map(process_sample, tasks)

    except Exception as e:
        logger.error(f"Critical error in Abundance pipeline: {e}", exc_info=True)
    finally:
        logger.info(f"Abundance analysis finished in {time.time() - start_time:.2f}s.")
        stop_log_listener()


### ====== main ====== ###
# The main function now accepts 'k' directly instead of 'best_k' being derived.
def main():
    parser = argparse.ArgumentParser(description="16S rRNA Pipeline - Taxonomy Abundance", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    group_base = parser.add_argument_group("Global Settings")
    group_base.add_argument("-i", "--input", required=True, help="Root directory containing sample folders")
    group_base.add_argument("-o", "--output", required=True, help="Output root directory")
    group_base.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads")

    group_data = parser.add_argument_group("Data & Database Settings")
    group_data.add_argument("--db", required=True, help="Path to K-mer SQLite database")
    group_data.add_argument("-k", type=int, default=25, help="K-mer length to use from database")
    group_data.add_argument("--suffix", default=None, help="Explicitly specify the input file suffix. If not provided, the script automatically searches for standard names: '_merged_no_chi.fasta.gz' or '_filtered_no_chi.fasta.gz")

    args = parser.parse_args()
    run_pipeline(args)


if __name__ == '__main__':
    main()
