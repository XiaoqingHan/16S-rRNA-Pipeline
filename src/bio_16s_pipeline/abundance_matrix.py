import time
import logging
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path
from .logs import get_logger


logger = get_logger()


# transform abundance table, combine all samples into one matrix
def abundance_matrix(dir_path):
    base = Path(dir_path)
    all_files = list(base.rglob("*_genus_abundance.csv"))
    if not all_files:
        logger.error(f"No genus_abundance files in {dir_path}.")
        return pd.DataFrame()
    logger.info(f"Found {len(all_files)} genus_abundance files, starting constructing abundance matrix...")

    matrix_list = []   
    for f in all_files:
        try:
            df = pd.read_csv(f)
            if df.empty: 
                continue
            sid = df['Sample_ID'].iloc[0]
            sample_series = df.set_index('Genus')['Relative_Abundance']
            sample_series.name = sid
            matrix_list.append(sample_series)
        except Exception as e:
            logger.error(f"File error {f.name}: {e}")

    final_matrix = pd.concat(matrix_list, axis=1, sort=True).fillna(0)
    final_matrix.index.name = 'Genus'
    
    return final_matrix


# only taxa with non-zero counts in at least min_fraction of the samples are retained
def filter_min_samples(df, min_fraction):
    if df.empty:
        return df.copy(), 0
    n_samples = df.shape[1]
    min_count = max(1, int(n_samples * min_fraction))
    keep = df.gt(0).sum(axis=1) >= min_count
    n_filtered = (~keep).sum()
    return df.loc[keep].copy(), int(n_filtered)


# plot the composition of each sample
def barplot(df, top_n):
    if df.empty:
        return None

    mean_abundance = df.mean(axis=1).sort_values(ascending=False)
    keep_genera = mean_abundance.head(top_n).index.tolist()

    df_top = df.loc[keep_genera].copy()
    others_sum = df.drop(keep_genera).sum(axis=0)
    if not others_sum.empty:
        df_top.loc['Others'] = others_sum

    df_plot = df_top.T

    width = max(12.0, len(df.columns) * 0.5)
    fig, ax = plt.subplots(figsize=(width, 8))
    df_plot.plot(kind='bar', stacked=True, colormap='tab20', edgecolor='white', width=0.8, alpha=0.9, ax=ax)
    ax.set_title("Genus Composition")
    ax.set_ylabel('Relative Abundance (%)', fontsize=8)
    plt.xticks(rotation=45, ha='right', fontsize=8)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.02, 1), loc='upper left', title="Genus", fontsize=8)
    plt.tight_layout()
    
    return fig


def run_analysis(args):
    out_path = Path(args.output)
    out_path.mkdir(parents=True, exist_ok=True)
    
    # construct matrix
    df_gen = abundance_matrix(args.input)
    if df_gen.empty: 
        return
    df_gen.to_csv(out_path / "all_samples_genus_matrix.csv", index=True)

    # filter low abundance
    filtered_df, n_filt = filter_min_samples(df_gen, args.min_fraction)
    logger.info(f"Filtered out {n_filt} low-frequency genera.")

    fig_bar = barplot(filtered_df, args.top_n)
    if fig_bar:
        fig_bar.savefig(out_path / "composition.png", dpi=300, bbox_inches='tight')
        plt.close(fig_bar)
    

def main():
    start_time = time.time()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    parser = argparse.ArgumentParser(description="16S rRNA Pipeline - Genus abundance and Visualization", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", required=True, help="Root directory containing genus abundance files")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("--top_n", type=int, default=20, help="Top N genera for composition plots")
    parser.add_argument("--min_fraction", type=float, default=0.1, help="Min sample fraction to retain genus")

    args = parser.parse_args()

    try:
        run_analysis(args)
        duration = time.time() - start_time
        logger.info(f"Analysis complete. Total time elapsed: {duration:.2f} seconds.")
    except Exception as e:
        logger.error(f"Pipeline failed due to an error: {e}")


if __name__ == "__main__":
    main()
