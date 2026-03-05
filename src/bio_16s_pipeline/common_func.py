import gzip
from pathlib import Path
from .logs import get_logger

logger = get_logger()


###########################
def read_fasta(filepath):
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"the input file {filepath} doesn't exist!")
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        header, seq = '', ''
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    yield header, seq
                header, seq = line, ''
            else:
                seq += line
        if header:
            yield header, seq


def read_fastq(filepath):
    path = Path(filepath)
    opener = gzip.open if filepath.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            if not header.startswith('@'):
                raise ValueError(f"Invalid FASTQ header: '{header}' (must start with '@')")
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            # if file is incomplete
            if not qual:
                break
            yield header, seq, qual


def save_files(outfile, reads):
    outfile = Path(outfile)
    is_gz = outfile.suffix == '.gz'
    opener = gzip.open if is_gz else open
    mode = 'wt' if is_gz else 'w'

    with opener(outfile, mode) as f:
        for read in reads:
            # fastq format
            if len(read) == 3:
                header, seq, qual = read
                f.write(f"{header}\n{seq}\n+\n{qual}\n")
            # fasta format
            else:
                header, seq = read
                header = '>' + header.lstrip('@>').strip()
                f.write(f"{header}\n{seq}\n")


# generate kmers
def yield_kmers(seq, k):
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if 'N' not in kmer:
            yield kmer


# get reversed, complementary sequence of reads 2, including IUPAC
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
                  'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B'}
    # let abnomal base be N 
    return ''.join(complement.get(base.upper(), 'N') for base in reversed(seq))


# get reversed seq for simplified seqs
def reverse_sim(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base.upper(), 'N') for base in reversed(seq))

  
# generate a group-sample dictionary mapping
def build_group_mapping(project_path, groups):
    project_path = Path(project_path)
    if not groups:
        raise ValueError("No groups specified.")

    mapping = {}
    for group in groups:
        group_dir = project_path / group
        if not group_dir.exists() or not group_dir.is_dir():
            logger.warning(f"Group directory {group_dir} not found, skipping.")
            continue

        for sample_dir in sorted(group_dir.iterdir()):
            if sample_dir.is_dir():
                sample_name = sample_dir.name
                mapping[sample_name] = group

    return mapping

