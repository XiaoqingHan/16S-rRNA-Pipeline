# db_utils.py
import sqlite3
from contextlib import contextmanager
from collections import defaultdict
from .logs import get_logger

@contextmanager
def open_kmer_db(path):
    logger = get_logger()
    conn = None
    try:
        conn = sqlite3.connect(path)
        yield conn
    except sqlite3.Error as e:
        logger.error(f"Database error: {e}")
        raise
    finally:
        if conn:
            conn.close()

def fetch_all_kmer_labels(conn, k, table=None):
    if table is None:
        table = f'kmer_labels_k{k}'
    cursor = conn.cursor()
    cursor.execute(f"SELECT kmer, label FROM {table}")
    kmer_map = defaultdict(list)
    for kmer, label in cursor.fetchall():
        kmer_map[kmer].append(label.strip())
    return dict(kmer_map)

def fetch_labels_for_kmer(conn, k, kmer):
    table = f'kmer_labels_k{k}'
    cursor = conn.cursor()
    cursor.execute(f"SELECT label FROM {table} WHERE kmer = ?", (kmer,))
    return [row[0] for row in cursor.fetchall()]
