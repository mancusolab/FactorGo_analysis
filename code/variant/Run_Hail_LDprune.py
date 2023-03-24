# PYSPARK_SUBMIT_ARGS='--driver-memory 400g --executor-memory 400g pyspark-shell'
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-dat", help="Prune which dat?")
parser.add_argument("-chr", help="Prune which chr?")
args = parser.parse_args()

outdir = "./data/PanUKBB/3_ldprune/" + args.dat + "/pruned/"
indir = "./data/PanUKBB/3_ldprune/" + args.dat + "/toprune/"
chr = args.chr

import time
from math import ceil

# from ukbb_pan_ancestry import *
import hail as hl
import numpy as np
import pandas as pd
from hail.linalg import BlockMatrix


def tie_breaker(left, right):
    return hl.sign(right.maf - left.maf)


## create windown size of locus 250kb
def ld_prune(ht_idx, r2_bm, bp_window_size, keep_higher_maf=True):
    ht_idx = ht_idx.add_index()
    _, stops = hl.linalg.utils.locus_windows(ht_idx.locus, bp_window_size)
    print("done assigning locus window")
    ## expand entry of bm matrix into (i,j)
    entries = r2_bm.sparsify_row_intervals(
        range(stops.size), stops, blocks_only=True
    ).entries(keyed=False)
    print("done sparsify")
    ## filter to variant pairs with r2 higher than torelant
    entries = entries.filter((entries.entry >= 0.3) & (entries.i < entries.j))
    ## remove r2 entry column
    entries = entries.select(i=hl.int32(entries.i), j=hl.int32(entries.j))

    if keep_higher_maf:
        fields = ["maf", "locus"]
    else:
        fields = ["locus"]

    info = ht_idx.aggregate(
        hl.agg.collect(ht_idx.row.select("idx", *fields)), _localize=False
    )
    info = hl.sorted(info, key=lambda x: x.idx)
    entries = entries.annotate_globals(info=info)

    entries = entries.filter(
        (entries.info[entries.i].locus.contig == entries.info[entries.j].locus.contig)
        & (
            entries.info[entries.j].locus.position
            - entries.info[entries.i].locus.position
            <= bp_window_size
        )
    )

    entries = entries.annotate(
        i=hl.struct(idx=entries.i, maf=entries.info[entries.i].maf),
        j=hl.struct(idx=entries.j, maf=entries.info[entries.j].maf),
    )
    print("done preparing entries")

    variants_to_remove = hl.maximal_independent_set(
        entries.i, entries.j, keep=False, tie_breaker=tie_breaker, keyed=False
    )
    print("done find variants to remove")

    ht_idx = ht_idx.annotate_globals(
        variants_to_remove=variants_to_remove.aggregate(
            hl.agg.collect_as_set(variants_to_remove.node.idx), _localize=False
        )
    )

    ## remove rows that is pruned out
    ht_idx_pruned = ht_idx.filter(
        ht_idx.variants_to_remove.contains(hl.int32(ht_idx.idx)), keep=False
    ).select()

    print("Preparing pd to output.")
    ht_idx_pruned = ht_idx_pruned.to_pandas()
    if ht_idx_pruned.shape[0] < 1:
        ht_idx_pruned = pd.DataFrame()
    else:
        ht_idx_pruned["alleles"] = ht_idx_pruned["alleles"].str.join(":")

    return ht_idx_pruned


# for chr_idx in range(20):
hl.init(
    spark_conf={
        "spark.hadoop.fs.gs.requester.pays.mode": "CUSTOM",
        "spark.hadoop.fs.gs.requester.pays.buckets": "ukb-diverse-pops-public, pan-ukb-us-east-1",
        "spark.hadoop.fs.gs.requester.pays.project.id": "panukbb-2022",
        "spark.driver.memory": "400g",
        "spark.executor.memory": "400g",
    },
    default_reference="GRCh37",
    min_block_size=128,
)

# LD block matrix table for EUR
ldblock = "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm"
# this is pearson R (not R2!)
bm = BlockMatrix.read(ldblock)
## Variant index Hail Table for EUR: which row/column corresponds to which variant
ldvar = "s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.variant.ht"
ht_idx = hl.read_table(ldvar)
print(f"Process for chr {chr}")

## choose a short chr
ht = hl.import_table(indir + "chr" + chr + ".toprune")
ht = ht.transmute(**hl.parse_variant(ht.variant)).key_by("locus", "alleles")
# ht = ht.annotate(LDidx=hl.int32(ht.LDidx))
ht = ht.annotate(maf=hl.float32(ht.maf_EUR))
# biallelic_ht = ht.filter((hl.len(ht.alleles[0]) == 1) & (hl.len(ht.alleles[1]) == 1))

ht_idx = ht_idx.join(ht, "inner")

### !! restrict to biallelic variants
LDidx = ht_idx.idx.collect()  # take idx values of genetic variants
# checkidx = (len(LDidx) == len(np.unique(LDidx)))
# print(f"LD idx is unique: {checkidx}")

bm = bm.filter(LDidx, LDidx)
bm = bm ** 2
print(f"done preparing LD bm matrix with {len(LDidx)}")

pruned_ht = ld_prune(ht, bm, bp_window_size=250000, keep_higher_maf=True)
pruned_ht["variant"] = pruned_ht["locus"].astype("str") + ":" + pruned_ht["alleles"]
pruned_ht["variant"].to_csv(outdir + "chr" + chr + ".pruned", sep="\t", index=False)
print(f"done with chr {chr}")
hl.stop()
time.sleep(5)
