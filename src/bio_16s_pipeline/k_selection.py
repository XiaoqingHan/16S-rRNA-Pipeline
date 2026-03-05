import random
import time
import json
import argparse
import multiprocessing
import pandas as pd
import numpy as np
from typing import TextIO
from itertools import combinations
from collections import Counter
from pathlib import Path
from multiprocessing import Pool, cpu_count
from .db_utils import open_kmer_db, fetch_all_kmer_labels
from .common_func import read_fasta
from .logs import get_logger, start_log_listener, stop_log_listener, init_worker_logger

###############
# this script is run for multiple samples in multiple groups
# - inputs: ref.fasta, merged.fasta, unmerged.fasta
# construct kmer database - build_kmer_db()
# select samples randomly - select_samples()
# test kmer influence - kmer_impact()
# - outputs: top_n.csv, kmer_jaccard.csv
###############


### ====== main functions ====== ###
# calculate jaccard similarity
def jaccard(results):
    ks = sorted(results.keys())
    jacc_map = {k: {} for k in ks}
    for k1, k2 in combinations(ks, 2):
        set1 = set(results[k1]['top_taxon'])
        set2 = set(results[k2]['top_taxon'])
        inter = len(set1 & set2)
        union = len(set1 | set2)
        j = inter / union if union > 0 else 0.0
        jacc_map[k1][k2] = j
        jacc_map[k2][k1] = j
    for k in ks:
        jacc_map[k][k] = 1.0
    return jacc_map


# parallel
def process_batch(args):
    logger = get_logger()
    batch, k_list, kmer_maps = args
    taxon_counter = {k: Counter() for k in k_list}
    times = {k: 0.0 for k in k_list}
    # max k for dealing with multiple k
    max_k = max(k_list)
    for header, seq in batch:
        votes_per_k = {k: Counter() for k in k_list}
        n = len(seq)
        for i in range(n - max_k + 1):
            window = seq[i:i + max_k]
            for k in k_list:
                t0 = time.time()
                km = window[:k]
                labels = kmer_maps[k].get(km, [])
                if labels:
                    weight = 1.0 / len(labels)
                    for label in labels:
                        votes_per_k[k][label] += weight
                times[k] += time.time() - t0
        for k in k_list:
            dominant = votes_per_k[k].most_common(1)[0][0] if votes_per_k[k] else "Unclassified"
            taxon_counter[k][dominant] += 1
    return taxon_counter, times


def kmer_impact(seqs, kmer_maps, config):
    logger = get_logger()
    seq_list = list(seqs)
    k_list = config['klist']
    
    batches = [(seq_list[i:i + config['batch']], k_list, kmer_maps) for i in range(0, len(seq_list), config['batch'])]
    num_processes = min(cpu_count(), len(batches)) if batches else 1

    ctx = multiprocessing.get_context("spawn")
    q = ctx.Queue() 

    with ctx.Pool(processes=num_processes, initializer=init_worker_logger, initargs=(q,)) as pool:
        batch_results = pool.map(process_batch, batches)

    # aggregate partial counts and timings from all batches
    merged = {k: Counter() for k in k_list}
    k_times = {k: 0.0 for k in k_list}
    for br, bt in batch_results:
        for k in k_list:
            merged[k].update(br[k])
            k_times[k] += bt[k]

    results = {}
    for k in k_list:
        top = merged[k].most_common(config['top'])
        results[k] = {
            'top_taxon': [g for g, _ in top],
            'top_counts': [c for _, c in top],
            'unique_taxon': len(merged[k])
        }
        logger.info(f"Analysis for k={k} finished. Found {results[k]['unique_taxon']} taxon.")

    return results, k_times


def run_pipeline(config):
    ctx = multiprocessing.get_context("spawn")
    q = ctx.Queue()
    start_log_listener(queue=q)
    logger = get_logger()
    logger.info("Starting K-selection analysis.")

    start_time = time.time()
    res_dir = Path(config['output'])
    res_dir.mkdir(parents=True, exist_ok=True)

    try:
        logger.info(f"K-selection list to be evaluated: {config['klist']}")

        kmer_maps = {}
        with open_kmer_db(config['db']) as conn:
            for k in config['klist']:
                table_name = f"kmer_labels_k{k}"
                kmer_maps[k] = fetch_all_kmer_labels(conn, k, table_name)

        if config['mode'] == "single":
            target_file = "*_filtered.fasta.gz"
        else:
            target_file = "*_merged.fasta.gz"

        input_path = Path(config['input'])
        found_files = list(input_path.rglob(target_file))

        if not found_files:
            logger.error("No valid sequence files found. Please check input path and mode.")
            return

        # determine sampling iterations, if robust is False, use a list with only the base seed
        if not config['robust']:
            iteration_seeds = [config['seed']]
            logger.info("Mode: Simple sampling (Single seed)")
        else:
            random.seed(config['seed'])
            iteration_seeds = [random.randint(0, 10000) for _ in range(config['n_seeds'])]
            logger.info(f"Mode: Robust evaluation ({config['n_seeds']} seeds)")

        all_metrics = []
        for idx, s in enumerate(iteration_seeds):
            logger.info(f"Processing seed set {idx + 1}/{len(iteration_seeds)}: {s}")
            random.seed(s)
            test_files = random.sample(found_files, min(config['n_sample'], len(found_files)))

            for f_path in test_files:
                sample_name = f_path.name.split('_')[0]
                logger.info(f"Running k-impact for sample: {sample_name}")
            
                seqs = read_fasta(f_path)
                results, k_times = kmer_impact(seqs, kmer_maps, config)
    
                for k in config['klist']:
                    all_metrics.append({
                        "seed": s,
                        "k": k,
                        "sample": sample_name,
                        "unique_taxa": results[k]["unique_taxon"],
                        "time_sec": k_times[k]
                    })

        df_res = pd.DataFrame(all_metrics)
        df_res.to_csv(res_dir / "k_selection_results.csv", index=False)
    
        summary = {}
        for k in config['klist']:
            avg_t = df_res[df_res['k'] == k]['time_sec'].mean()
            summary[k] = {"avg_time_sec": round(avg_t, 2)}
            logger.info(f"k={k} average computation time: {avg_t:.2f}s")

        with open(res_dir / "k_selection_summary.json", "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=4)

    except Exception as e:
        logger.error(f"Error during pipeline execution: {e}")
        raise
    finally:
        elapsed_time = time.time() - start_time
        logger.info(f"K-selection completed in {elapsed_time:.2f}s")
        stop_log_listener()


def main():
    parser = argparse.ArgumentParser(description="K-mer optimal length selection", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Core paths
    parser.add_argument("-i", "--input", required=True, help="Input directory")
    parser.add_argument("-o", "--output", default="./kmer_results", help="Output directory")
    parser.add_argument("-d", "--db", required=True, help="K-mer database path")
    
    # Evaluation settings
    parser.add_argument("-k", "--klist", type=int, nargs='+', default=[21, 25, 31], help="k-mer lengths")
    parser.add_argument("-n", "--n_sample", type=int, default=5, help="Number of samples to test per iteration")
    parser.add_argument("--top", type=int, default=50, help="Number of top abundant taxa to consider for Jaccard similarity")

    # Robustness toggle
    parser.add_argument("--robust", action="store_true", help="Enable multi-seed robust evaluation (True/False)")
    parser.add_argument("-s", "--n_seeds", type=int, default=3, help="Number of random seeds if robust is enabled")
    parser.add_argument("--seed", type=int, default=42, help="Initial random seed for reproducibility")
    
    # Processing settings
    parser.add_argument("--mode", choices=['single', 'paired'], default='paired', help="Sequencing mode")
    parser.add_argument("--batch", type=int, default=5000, help="Number of sequences processed per parallel batch")

    args = parser.parse_args()
    run_pipeline(vars(args))


if __name__ == '__main__':
    main()
