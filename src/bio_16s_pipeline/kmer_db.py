import time
import gzip
import sqlite3
import argparse
from pathlib import Path
from .common_func import read_fasta, reverse_sim, yield_kmers
from .logs import get_logger, start_log_listener, stop_log_listener



###############
# this script is used for constructing kmer database by using reference data
###############


logger = get_logger()


### ====== main functions ====== ###
# get kmer table for single sequence
def process_sequence(seq, source, label, k_list, cursor):
    if seq and source and label:
        seq_variants = [seq.upper(), reverse_sim(seq)]
        for k in k_list:
            table_name = f"kmer_labels_k{k}"
            kmers = set()
            for s in seq_variants:
                kmers.update(yield_kmers(s, k))
            if kmers:
                cursor.executemany(
                    f'INSERT OR IGNORE INTO {table_name} VALUES (?,?,?)',
                    [(kmer, label, source) for kmer in kmers])


# create table structure, do not create index 
def create_table(cursor, k_list):
    for k in k_list:
        table_name = f"kmer_labels_k{k}"
        cursor.execute(f"""
            CREATE TABLE IF NOT EXISTS {table_name}(
                kmer TEXT NOT NULL,
                label TEXT NOT NULL, 
                source TEXT,
                PRIMARY KEY (kmer, label))
        """)


# construct kmer SQLite database, easy to query
def build_kmer_db(ref_fasta, kmer_db, k_list):
    start_log_listener()
    start_time = time.time()

    db_path = Path(kmer_db)
    if not db_path.parent.exists():
        logger.info(f"Creating directory: {db_path.parent}")
        db_path.parent.mkdir(parents=True, exist_ok=True)
   
    try:
        with sqlite3.connect(kmer_db) as conn:
            cursor = conn.cursor()
            cursor.execute("PRAGMA synchronous = OFF")
            cursor.execute("PRAGMA journal_mode = MEMORY")

            create_table(cursor, k_list)

            if ref_fasta:
                fasta_path = Path(ref_fasta)
                for header, seq in read_fasta(fasta_path):
                    parts = header.split()
                    if len(parts) >= 3:
                        source = parts[0]
                        label = ' '.join(parts[2:])
                        process_sequence(seq, source, label, k_list, cursor)
                    else:
                        logger.error(f"Invalid header skipped: {header}")
                conn.commit()

            for k in k_list:
                table_name = f"kmer_labels_k{k}"
                logger.info(f"Creating K-mer index...")
                cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_kmer_{k} ON {table_name}(kmer)")

                cursor.execute(f"SELECT COUNT(DISTINCT kmer), COUNT(DISTINCT label) FROM kmer_labels_k{k}")
                kmer_count, label_count = cursor.fetchone()
                logger.info(f"k={k}: {kmer_count:,} kmers -> {label_count:,} labels")

            conn.commit()

    finally:
        logger.info(f"K-mer database construction finished in {time.time() - start_time:.2f}s")
        stop_log_listener()


def main():

    parser = argparse.ArgumentParser(description="Construct k-mer database from reference FASTA")
    parser.add_argument("-i", "--input", required=True, help="Reference FASTA file (.fa, .fasta, .gz)")
    parser.add_argument("-d", "--db", required=True, help="Output SQLite database path")
    parser.add_argument("-k", "--klist", type=int, nargs='+', default=[25], help="List of k values (e.g., -k 21 31)")
    
    args = parser.parse_args()

    build_kmer_db(args.input, args.db, args.klist)


if __name__ == '__main__':
    main()

