#! /usr/bin/env python
import argparse as ap
import logging
import os
import sys

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.decomposition import TruncatedSVD


def get_logger(name, path=None):
    logger = logging.getLogger(name)
    if not logger.handlers:
        # Prevent logging from propagating to the root logger
        logger.propagate = 0
        console = logging.StreamHandler()
        logger.addHandler(console)

        log_format = "[%(asctime)s - %(levelname)s] %(message)s"
        date_format = "%Y-%m-%d %H:%M:%S"
        formatter = logging.Formatter(fmt=log_format, datefmt=date_format)
        console.setFormatter(formatter)

        if path is not None:
            disk_log_stream = open("{}.log".format(path), "w")
            disk_handler = logging.StreamHandler(disk_log_stream)
            logger.addHandler(disk_handler)
            disk_handler.setFormatter(formatter)

    return logger


def tsvd(dat, k, fast, seed, log):
    if fast:
        log.info("Begin TSVD using sklearn")
        matt = TruncatedSVD(n_components=k, n_iter=20, random_state=seed)
        US = matt.fit_transform(csr_matrix(dat))
        U = US / (matt.singular_values_)
        V = (matt.components_).T
        D = matt.singular_values_
    else:
        log.info("Begin TSVD using built-in svd")
        U, D, Vh = np.linalg.svd(dat, full_matrices=False)
        U = U[:, 0:k]
        D = D[0:k]
        V = Vh.T
        V = V[:, 0:k]
        matt = None

    return U, D, V, matt


def read_dat(dat_path, divideN, N_path, center, scale, log):
    # input data is: SNP, trait1, trait2,...,traitN (pxn)
    df_z = pd.read_csv(dat_path, delimiter="\t", header=0)
    p_snps, n_studies = df_z.shape
    log.info(f"Found N = {n_studies-1} studies, P = {p_snps} SNPs")

    snp_col = df_z.columns[0]
    # df_z = df_z[df_z[snp_col].isin(df_s2g['rsid'].values)]
    # p, n = df_z.shape
    # log.info(f"Keep variants with S2G, now dfZ matrice shape: {p}x{n}")

    # drop the first column (axis = 1)
    df_z.drop(labels=[snp_col], axis=1, inplace=True)
    df_z = df_z.astype("float")  # pxn

    if divideN is True:
        # read sample size file and convert str into numerics, convert to nxp matrix
        df_N = pd.read_csv(N_path, delimiter="\t", header=0)
        df_N = df_N.astype("float")
        # convert sampleN (a file with one column and header)to arrays
        N_col = df_N.columns[0]
        sampleN = df_N[N_col].values
        sampleN_sqrt = np.sqrt(sampleN)
        df_z = df_z / sampleN_sqrt

    # center data
    if center == "traits":
        df_z = df_z.subtract(df_z.mean())
        log.info("center within each trait")
    elif center == "snps":
        df_z_t = df_z.T
        df_z = (df_z_t.subtract(df_z_t.mean())).T
        log.info("center within each snp")
    else:
        pass

    if scale == "traits":
        df_z = df_z.divide(df_z.std())
        log.info("scale within each trait")
    elif scale == "snps":
        df_z_t = df_z.T
        df_z = (df_z_t.divide(df_z_t.std())).T
        log.info("scale within each snp")
    else:
        pass

    return df_z


def R2(B, W_m, Z_m):
    # import pdb; pdb.set_trace()
    n, p = B.shape
    _, k = Z_m.shape

    tss = np.sum(B * B)

    sse = np.zeros((k,))
    for i in range(n):
        WZ = W_m * Z_m[i]  # pxk
        res = B[i][:, None] - WZ  # pxk
        tmp = np.sum(res * res, axis=0)  # (k,)
        sse += tmp
    r2 = 1.0 - sse / tss

    return r2


def main(args):
    argp = ap.ArgumentParser(description="")  # create an instance
    argp.add_argument("Zscore_path")
    argp.add_argument("N_path")
    argp.add_argument(
        "-k", type=int, default=10
    )  # "-" must only has one letter like "-k", not like "-knum"
    argp.add_argument(
        "--center",
        choices=["traits", "snps", "none"],
        default="none",
        help="How to center data: [traits, snps, or none]. Default is none",
    )
    argp.add_argument(
        "--scale",
        choices=["traits", "snps", "none"],
        default="none",
        help="How to standardize data: [traits, snps, or none]. Default is none",
    )
    argp.add_argument(
        "--sklearn",
        action="store_true",
        help="How to run tsvd: either the built-in function(for small data), or fast (sklean.decomposition.TruncatedSVD)",
    )
    argp.add_argument(
        "--divideN",
        action="store_true",
        help="divide Z score by square root of N",
    )
    argp.add_argument(
        "--seed",
        default=24983,
        type=int,
        help="Maximum number of iterations to learn parameters",
    )
    argp.add_argument(
        "-o", "--output", type=str, default="TSVDres", help="Prefix path for output"
    )

    args = argp.parse_args(args)  # a list a strings

    log = get_logger(__name__, args.output)
    log.setLevel(logging.INFO)

    Zdf = read_dat(
        args.Zscore_path, args.divideN, args.N_path, args.center, args.scale, log
    )

    Zarr = (np.array(Zdf)).T
    log.info(f"Finished pre-process with {Zarr.shape} data matrice.")

    check_na_inf = np.sum(np.sum(np.asarray_chkfinite(Zarr)))
    if np.isinf(check_na_inf) or np.isnan(check_na_inf):
        sys.exit("Processed data contains NA or Inf.")

    U, D, V, matt_res = tsvd(Zarr, int(args.k), args.sklearn, args.seed, log)
    log.info("Writing results.")

    vbr2 = R2(Zarr, V, U * D)

    np.savetxt(
        f"{args.output}.U.divideN${args.divideN}.tsv.gz", U, fmt="%s", delimiter="\t"
    )
    np.savetxt(
        f"{args.output}.D.divideN${args.divideN}.tsv.gz", D, fmt="%s", delimiter="\t"
    )
    np.savetxt(
        f"{args.output}.V.divideN${args.divideN}.tsv.gz", V, fmt="%s", delimiter="\t"
    )
    np.savetxt(
        f"{args.output}.r2.divideN${args.divideN}.tsv.gz",
        matt_res.explained_variance_ratio_,
        fmt="%s",
        delimiter="\t",
    )
    np.savetxt(
        f"{args.output}.vbr2.divideN${args.divideN}.tsv.gz",
        vbr2,
        fmt="%s",
        delimiter="\t",
    )

    log.info("Finished writing.")

    return 0


# user call this script will treat it like a program
if __name__ == "__main__":
    sys.exit(
        main(sys.argv[1:])
    )  # grab all arguments; first arg is alway the name of the script
