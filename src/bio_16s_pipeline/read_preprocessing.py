import time
import re
import multiprocessing
import argparse
from pathlib import Path
from .common_func import read_fastq, save_files, reverse_complement
from .logs import get_logger, start_log_listener, stop_log_listener, init_worker_logger


###############
# this script is run for all samples in a group
# - inputs: fastq.gz
# remove primer - trim()
# quality control - qc()
# merge two sequences - merge_paired_reads()
# qc merged seqs - merge_qc()
# - outputs: merge.fasta, unmerged.fasta
###############


### ====== main functions ====== ###
# calculate average phred score per read, phred=ascii-33 (Illumina >= 1.8)
def phred_score(base_quality):
    if not base_quality or len(base_quality) == 0:
        return 0
    return sum(ord(i) - 33 for i in base_quality) / len(base_quality)


# calculate the mismatch number of base
def hamming_distance(seq1, seq2):
    count = 0
    if len(seq1) == len(seq2):
        for base1, base2 in zip(seq1, seq2):
            if base1 != base2:
                count += 1
    return count


def primer_to_regex(primer):
    iupac_map = {
        'A':'A', 'C':'C', 'G':'G', 'T':'T', 'U':'T',
        'R':'[AG]', 'Y':'[CT]', 'S':'[GC]', 'W':'[AT]',
        'K':'[GT]', 'M':'[AC]', 'B':'[CGT]', 'D':'[AGT]',
        'H':'[ACT]', 'V':'[ACG]', 'N':'[ACGT]'
    }
    return "".join([iupac_map.get(b.upper(), b.upper()) for b in primer])


# remove primer, no need to consider 3'(seq[-len(primer):]) cuz primer usually happens in 5', allow 2 mismatch
# primer not always starts from first position, need slide
def detect_primer(seq, qual, primer, max_offset):
    seq = seq.upper()
    primer = primer.upper()
    search_limit = max_offset + len(primer)
    
    # use full length of primer
    full_pattern = primer_to_regex(primer)
    match = re.search(full_pattern, seq[:search_limit])
 
    if match:
        cut_site = match.end()
        return seq[cut_site:], qual[cut_site:], True
 
    # detect the tail of primer, usually 12-15bp, conserved seq
    anchor_primer = primer[-12:]
    anchor_pattern = primer_to_regex(anchor_primer)
    match_a = re.search(anchor_pattern, seq[:search_limit])
    
    if match_a:
        cut_site = match_a.end()
        return seq[cut_site:], qual[cut_site:], True

    # if search nothing, return original seq  
    return seq, qual, False


# trim primer
def trim(reads, primer, max_offset):
    trimmed_reads = []
    trimmed_cnt = 0

    for header, seq, qual in reads:
        s, q, trimmed_flag = detect_primer(seq, qual, primer, max_offset)
        if trimmed_flag:
            trimmed_cnt += 1
            trimmed_reads.append((header, s, q))

    return trimmed_cnt, trimmed_reads


# quality control, use Q30 strategy
def qc(reads, min_len, min_q, max_n, label):
    logger = get_logger()
    filtered_reads = []
    passed = 0
    filter_reasons = {"low_quality": 0, "short_length": 0, "high_N": 0}

    for idx, (header, seq, qual) in enumerate(reads):
        n_content = seq.count('N') / len(seq)
        avg_q = phred_score(qual)
        if len(seq) >= min_len and avg_q >= min_q and n_content <= max_n:
            filtered_reads.append((header, seq, qual))
            passed += 1
        else:
            if avg_q < min_q:
                filter_reasons["low_quality"] += 1
            elif len(seq) < min_len:
                filter_reasons["short_length"] += 1
            elif n_content > max_n:
                filter_reasons["high_N"] += 1

        # debug
        if idx < 5:
            logger.debug(f"Read {idx}: len={len(seq)}, avg_q={avg_q:.1f}, n_content={n_content:.3f}")

    logger.info(f"{label} QC breakdown: {filter_reasons}")
    return passed, filtered_reads


# normalize headers
def normal_header(h):
    h = h.lstrip('@>').split()[0]
    return h


# merge reads1 and reads2, allow 1 mismatch, R1+reverse_complement(R2)
def merge_paired_reads(reads1, reads2, min_overlap, max_merge_rate):
    logger = get_logger()
    reads1 = list(reads1)
    reads2 = list(reads2)
    header_map1 = {normal_header(h): h for h, _, _ in reads1}
    header_map2 = {normal_header(h): h for h, _, _ in reads2}

    seq1 = {}
    for h, r, _ in reads1:
        norm_h = normal_header(h)
        if norm_h not in seq1:
            seq1[norm_h] = r.upper()
    seq2 = {}
    for h, r, _ in reads2:
        norm_h = normal_header(h)
        if norm_h not in seq2:
            seq2[norm_h] = r.upper()

    common_headers = []
    for norm_h in seq1:
        if norm_h in seq2:
            orig_h1 = header_map1[norm_h]
            orig_h2 = header_map2[norm_h]
            if orig_h1.split()[0] == orig_h2.split()[0]:
                common_headers.append(norm_h)
    for norm_h in common_headers[:3]:
        logger.debug("Normalized common headers: %r ↔ %r", header_map1[norm_h], header_map2[norm_h])

    merged_seqs = {}
    unmerged1 = {}
    unmerged2 = {}
    for h in common_headers:
        s1, s2 = seq1[h], seq2[h]
        merged = False
        max_search_overlap = int(min(len(s1), len(s2)))
        for i in range(max_search_overlap, min_overlap - 1, -1):
            if hamming_distance(s1[-i:], s2[:i]) / i <= max_merge_rate:
                merged_seqs[f">{h}"] = s1 + s2[i:]
                merged = True
                break
        if not merged:
            unmerged1[f">{h}/1"] = s1
            #unmerged2[f">{h}/2"] = reverse_complement(s2)

    merged_list = [(header, seq) for header, seq in merged_seqs.items()]
    unmerged1_list = [(header, seq) for header, seq in unmerged1.items()]
    unmerged2_list = [(header, seq) for header, seq in unmerged2.items()]

    return merged_list, unmerged1_list, unmerged2_list


# check the quality of merged seqs
def merge_qc(merged_seqs, min_len_m, max_len_m, max_n_m):
    filtered = []
    for header, seq in merged_seqs:
        seq_len = len(seq)
        n_content = seq.count('N') / seq_len
        if (min_len_m <= seq_len <= max_len_m) and (n_content <= max_n_m):
            filtered.append((header, seq))
    return filtered


# process a single sample, r1 and r2 for paired-end, if single-end, r2 is None.
def process_sample(r1_path, r2_path, out_dir, config):
    logger = get_logger()
    out_dir = Path(out_dir)

    sample_name = out_dir.name
    lp = f"[{sample_name}]"

    trim_primer_flag = config["trim_primer"]

    fastq_r1 = Path(r1_path)
    if not fastq_r1.exists():
        logger.error(f"{lp} R1 file does not exist.")
        return
    if fastq_r1.stat().st_size == 0:
        logger.warning(f"{lp} R1 file is empty.")  
        return

    is_paired_end = bool(r2_path)
    if is_paired_end:
        fastq_r2 = Path(r2_path)
        if not fastq_r2.exists():
            logger.error(f"{lp} R2 file does not exist.")
            return
        if fastq_r2.stat().st_size == 0:
            logger.warning(f"{lp} R2 file is empty.")
            return
  
    out_dir.mkdir(parents=True, exist_ok=True)

    if is_paired_end:
        out_merged = out_dir / f"{sample_name}_merged.fasta.gz"
        out_unmerged1 = out_dir / f"{sample_name}_unmerged_1.fasta.gz"
        out_unmerged2 = out_dir / f"{sample_name}_unmerged_2.fasta.gz"

        # skip if file exists	
        if out_merged.exists():
            logger.info(f"{sample_name}: {out_merged.name} exists, skip.")
            return

    else:
        out_qc = out_dir / f"{sample_name}_filtered.fasta.gz"
        if out_qc.exists():
            logger.info(f"{sample_name}: {out_qc.name} exists, skip.")
            return

    reads1 = list(read_fastq(fastq_r1))
    count1 = len(reads1)

    if is_paired_end:
        reads2 = list(read_fastq(fastq_r2))
        count2 = len(reads2)
        logger.info(f"{lp} Processed R1: {count1}, R2: {count2} read pairs.")
    else:
        logger.info(f"{lp} Processed: {count1} reads.")

    # trim primers if enabled
    if trim_primer_flag:
        if is_paired_end:    
            trim_cnt1, trim1 = trim(reads1, config["p1"], config["max_offset"])
            trim_cnt2, trim2 = trim(reads2, config["p2"], config["max_offset"])
            logger.info(f"{lp} R1 trimmed: {trim_cnt1} reads; R2 trimmed: {trim_cnt2} reads.")
        else:
            trim_cnt1, trim1 = trim(reads1, config["p1"], config["max_offset"])
            logger.info(f"{lp} trimmed: {trim_cnt1} reads.")

    else:
        if is_paired_end:
            trim1, trim2 = reads1, reads2
        else:
            trim1 = reads1
        logger.info(f"{lp} Skipping primer trimming.")

    # quality control
    if is_paired_end:
        passed1, filt1 = qc(trim1, config["min_len"], config["min_q"], config["max_n"], label=f"{lp} R1")
        logger.info(f"{lp} R1 QC Result: {passed1} passed.")
        passed2, filt2 = qc(trim2, config["min_len"], config["min_q"], config["max_n"], label=f"{lp} R2")
        logger.info(f"{lp} R2 QC Result: {passed2} passed.")        

        # reverse complement R2 for merging
        filt2_rev = [(h, reverse_complement(s), q) for h, s, q in filt2]

        # merge paired reads
        merged, un1, un2 = merge_paired_reads(filt1, filt2_rev, config["min_overlap"], config["max_merge_rate"])
        logger.info(f"{lp} Merge success: {len(merged)} pairs.")

        # post-merge QC
        qc_merged = merge_qc(merged, config["min_len_m"], config["max_len_m"], config["max_n_m"])
        logger.info(f"{lp} Post-merge QC success: {len(qc_merged)} pairs.")

        if qc_merged:
            lengths = [len(seq) for _, seq in qc_merged]
            logger.info(f"{lp} Final length stats: min={min(lengths)}, max={max(lengths)}, avg={sum(lengths) / len(lengths):.1f}bp.")
        else:
            logger.info(f"{lp} No merged sequences passed post-merge QC.")

        # save output files
        save_files(out_merged, qc_merged)
        #save_files(out_unmerged1, un1)
        #save_files(out_unmerged2, un2)

    else:
        # QC
        passed, filt = qc(trim1, config["min_len"], config["min_q"], config["max_n"], lp)
        logger.info(f"{lp} QC Result: {passed} reads passed.")        

        #fasta_reads = [(r[0].replace('@', '>'), seq) for h, seq, _ in filt]
        fasta_reads = [(f">{normal_header(h)}", seq) for h, seq, _ in filt]
        save_files(out_qc, fasta_reads)

    logger.info(f"{lp} Done.")


# paralleling
def process_wrapper(r1, r2, out_dir, config):
    logger = get_logger()
    try:
        process_sample(r1, r2, out_dir, config)
    except Exception as e:
        logger.error(f"Sample processing failed: {out_dir}, Error: {e}", exc_info=True)


def run_pipeline(config):
    ctx = multiprocessing.get_context("spawn")
    q = ctx.Queue()
    start_log_listener(queue=q)
    logger = get_logger()

    logger.info("Starting 16S rRNA read preprocessing analysis.")
    start = time.time()

    try:
        input_path = Path(config["input"]) 
        output_path = Path(config["output"])
        output_path.mkdir(parents=True, exist_ok=True)
	 
        seq_mode = config["mode"]
        n_procs = config["threads"]

        task_args = []
	
        # find all .fastq.gz files
        if seq_mode == "single":
            extensions = ("*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz")
            sample_files = []
            for ext in extensions:
                sample_files.extend(input_path.rglob(ext))
            sample_files = sorted(list(set(sample_files)))            
            logger.info(f"Single-end mode identification complete: {len(sample_files)} files found.")
            for f in sample_files:
                if f.name.endswith(".gz"):
                    # sample.fastq.gz -> sample.fastq -> sample
                    sample_name = Path(f.stem).stem  
                else:
                    sample_name = f.stem                
                out_dir = output_path / sample_name
                task_args.append((str(f), None, str(out_dir), config))

        # find _1.fastq.gz first, then find corresponding _2.fastq.gz
        elif seq_mode == "paired":
            r1_patterns = ["*_1.fastq.gz", "*_1.fastq", "*_1.fq.gz", "*_1.fq"]
            r1_files = []
            for pattern in r1_patterns:
                r1_files.extend(input_path.rglob(pattern))
            r1_files = sorted(list(set(r1_files)))
            for f1 in r1_files:
                # f1.suffixes can give all suffixes, ['.fastq', '.gz'] or ['.fastq']
                suffix_str = "".join(f1.suffixes)
                f2_str = str(f1).replace(f"_1{suffix_str}", f"_2{suffix_str}")
                f2 = Path(f2_str)
                if f2.exists():
                    sample_name = f1.name.replace(f"_1{suffix_str}", "")
                    out_dir = output_path / sample_name
                    task_args.append((str(f1), str(f2), str(out_dir), config))
                else:
                    logger.warning(f"Missing R2 file for {f1.name}, skipping sample.")
            logger.info(f"Paired-end mode identification complete: {len(task_args)} pairs found.")

        if not task_args:
            logger.error("No valid samples found. Please check your input path and mode.")
            return

        with ctx.Pool(processes=n_procs, initializer=init_worker_logger, initargs=(q,)) as pool:
            pool.starmap(process_wrapper, task_args)
    
    except Exception as e:
        logger.error(f"Main process failed: {e}", exc_info=True)
    finally:
        logger.info(f"All samples completed in {time.time() - start:.2f}s.")
        stop_log_listener()

def main():

    # create Parser
    parser = argparse.ArgumentParser(description="16S rRNA Pipeline - Trim, QC and Merge", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # add parameters, required=True means mandatory
    # --basics--
    group_base = parser.add_argument_group("Global Settings")
    group_base.add_argument("-i", "--input", required=True, help="Root path to Fastq files")
    group_base.add_argument("-m", "--mode", choices=["single", "paired"], default="paired", help="Sequencing mode")
    group_base.add_argument("-o", "--output", default="Preprocessing_results", help="Output root directory")
    group_base.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads for parallel processing")
    
    # --primer--
    group_primer = parser.add_argument_group("Primer Settings")
    # mutually exclusive group, decide execute or skip
    trim_action = group_primer.add_mutually_exclusive_group()
    trim_action.add_argument("--trim_primer", action="store_true", dest="trim_primer", default=True, help="Trimming primer")
    trim_action.add_argument("--skip_trim", action="store_false", dest="trim_primer", help="Skip trimming primer (default: False)")
    group_primer.add_argument("--p1", help="Forward Primer / Single-end Primer")
    group_primer.add_argument("--p2", help="Reverse Primer")
    group_primer.add_argument("--max_offset", type=int, default=25, help="Max displacement of primer search")

    # --QC--
    group_qc = parser.add_argument_group("Quality Control")
    group_qc.add_argument("--min_len", type=int, default=75, help="Min read length")
    group_qc.add_argument("--min_q", type=int, default=30, help="Min average phred score")
    group_qc.add_argument("--max_n", type=float, default=0.05, help="Max 'N' content rate")

    # --Merge--
    group_merge = parser.add_argument_group("Merge Settings")
    group_merge.add_argument("--min_overlap", type=int, default=20, help="Min overlap for merging")
    group_merge.add_argument("--max_merge_rate", type=float, default=0.05, help="Max mismatch rate in overlap")
    group_merge.add_argument("--min_len_m", type=int, default=200, help="Min length of merged sequence")
    group_merge.add_argument("--max_len_m", type=int, default=500, help="Max length of merged sequence")
    group_merge.add_argument("--max_n_m", type=float, default=0.01, help="Max 'N' rate of merged sequence")

    args = parser.parse_args()

    if args.trim_primer:
        if not args.p1:
            parser.error("Error: Trimming primer, please provide --p1 parameter.")
        if args.mode == "paired" and not args.p2:
            parser.error("Error: Paired-end mode, please provide --p2 parameter.")
    
    current_config = vars(args)

    run_pipeline(current_config)    


if __name__ == '__main__':
    main()

