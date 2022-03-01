import scanpy as sc
from anndata import read_h5ad
import pandas as pd
import numpy as np
import scipy as sp
import os
import time
import argparse
from statsmodels.stats.multitest import multipletests
import scdrs


def main(args):
    sys_start_time = time.time()

    ###########################################################################################
    ######                                    Parse Options                              ######
    ###########################################################################################
    H5AD_FILE = args.h5ad_file
    FLAG_FILTER = args.flag_filter == "True"
    FLAG_RAW_COUNT = args.flag_raw_count == "True"
    OUT_FILE = args.out_file

    ###########################################################################################
    ######                                     Load data                                 ######
    ###########################################################################################
    print("Load data:")

    # Load .h5ad file
    adata = read_h5ad(H5AD_FILE)
    if FLAG_FILTER:
        sc.pp.filter_cells(adata, min_genes=250)
        sc.pp.filter_genes(adata, min_cells=50)
    if FLAG_RAW_COUNT:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
    print(
        "--h5ad_file loaded: n_cell=%d, n_gene=%d (sys_time=%0.1fs)"
        % (adata.shape[0], adata.shape[1], time.time() - sys_start_time)
    )

    ###########################################################################################
    ######                                  Computation                                  ######
    ###########################################################################################

    # Preprocess
    scdrs.preprocess(adata)
    df_stats = adata.uns["SCDRS_PARAM"]["GENE_STATS"].copy()
    print("df_stats computed (sys_time=%0.1fs)" % (time.time() - sys_start_time))

    # zero_prop
    v = 1 - (adata.X != 0).mean(axis=0)
    temp_df = pd.DataFrame(
        index=adata.var_names, data={"zero_prop": np.array(v).flatten()}
    )
    df_stats = df_stats.join(temp_df)
    print(
        "zero_prop added to df_stats (sys_time=%0.1fs)" % (time.time() - sys_start_time)
    )

    # Write file
    df_stats.to_csv(OUT_FILE, sep="\t", index=True)
    print(
        "df_stats written: n_gene=%d (sys_time=%0.1fs)"
        % (df_stats.shape[0], time.time() - sys_start_time)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compute score")

    parser.add_argument("--h5ad_file", type=str, required=True)

    parser.add_argument(
        "--flag_filter",
        type=str,
        required=False,
        default="True",
        help="If to apply cell and gene filters to the h5ad_file data",
    )
    parser.add_argument(
        "--flag_raw_count",
        type=str,
        required=False,
        default="True",
        help="If True, apply size factor normalization and log1p transformation",
    )
    parser.add_argument(
        "--out_file",
        type=str,
        required=True,
        help="Save file at out_file",
    )

    args = parser.parse_args()

    main(args)
