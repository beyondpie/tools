"""
A temporary demo for splitting one bam file into multiple single-cell bam files.


# Need to set:
# ulimit -n 50000
# instead of 1200
# to allow multiple open files
"""
import os
import sys
from typing import List, Dict, Tuple, Any
from collections import OrderedDict
import string
import random
import re
import logging
from multiprocessing import Pool
import functools


import pandas as pd
import pysam

from sinto import filterbarcodes
from sinto import utils

class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    Ref:
    https://stackoverflow.com/questions/19425736/how-to-redirect-stdout-and-stderr-to-logger-in-python
    """

    def __init__(self, logger, level):
        self.logger = logger
        self.level = level
        self.linebuf = ""

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.level, line.rstrip())

    def flush(self):
        pass


def set_file_logger(
    fnm: str, fmode: str = "a", name: str = "sa2_pp",
    log_level: int = logging.DEBUG
) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(log_level)
    fh = logging.FileHandler(filename=fnm, mode=fmode)
    fm = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    fh.setFormatter(fm)
    logger.addHandler(fh)
    return logger
def get_barcodes_from_intervals_async(
        intervals: List[Tuple[str, int, int]],
        bam_fnm: str,
        index_fnm: str,
        barcode_pattern: str) -> List[str]:
    with pysam.AlignmentFile(bam_fnm,
                               index_filename = index_fnm,
                               mode = "rb") as abam:
        barcode = [re.search(
                    barcode_pattern, r.query_name).group().strip(":")
                for i in intervals
                for r in abam.fetch(contig = i[0], start = i[1], stop = i[2])]
    return barcode

def flatten_list(xss: List[List]) -> List:
    return [x for xs in xss for x in xs]
def unique_list(xs: List) -> List:
    return list(set(xs))


# * from snakemake or define by ourselves
logfnm: str = snakemake.log[0]
bam_fnm: str = snakemake.input['bam_fnm']
index_fnm: str = snakemake.input['index_fnm']
cell_meta_fnm: str = snakemake.input['cell_meta_fnm']
sublib_id: str = snakemake.params['sublib_id']
out_bam_dir: str = snakemake.params['out_bam_dir']
nthread: int = snakemake.threads
debug: int = False
overwrite: int = True

print("Under normal mode ...")
verbose = True

# * set logger
logger = set_file_logger(logfnm, name = "split bam")
## works in Linux, but have some troubles in mac
sys.stdout = StreamToLogger(logger = logger, level = logging.INFO)
sys.stderr = StreamToLogger(logger = logger, level = logging.ERROR)

# * functions
## FIXME: async function may not work in the same single file.
# If error happened, put this function into another file.
def get_barcodes_from_bam_async(
        bam_fnm:str,
        index_fnm: str,
        barcode_pattern: str = ':\w\w:\w\w:\w\w:',
        nproc: int = 2) -> List[str]:
    with pysam.AlignmentFile(
        bam_fnm, index_filename = index_fnm, mode = "rb") as abam:
        ## TODO: can we rewrite chunk_bam instead of using sinto package?
        intervals: Dict[int, List[Tuple[str, int, int]]] = utils.chunk_bam(
            abam, nproc, unmapped = True)
    p = Pool(nproc)
    barcodes = p.map_async(
        functools.partial(get_barcodes_from_intervals_async,
                          bam_fnm = bam_fnm,
                          index_fnm = index_fnm,
                          barcode_pattern = barcode_pattern),
        intervals.values()
    ).get(9999999)
    barcodes = flatten_list(barcodes)
    return(barcodes)

def write_oneread(oneread,
                  sublib_id,
                  barcodes: set,
                  outfhandles: Dict[str, Any],
                  logger: None|logging.Logger = None,
                  barcode_pattern: str = ':\w\w:\w\w:\w\w:',
                  verbose: bool = False) -> None:
    """
    1. get barcode, then get outfnm
    2. update_query_name by adding sublib_id
    3. write to the file
    """
    b = re.search(barcode_pattern,
                  oneread.query_name).group().strip(":")
    b2 = f"{sublib_id}:{b}"
    if b2 not in barcodes:
        if verbose:
            l = f"{b2} is not in the barcode list, ignored."
            if logger is None:
                print(l)
            else:
                logger.info(l)
        return None
    q = re.sub(
        barcode_pattern, f":{b2}:", oneread.query_name)
    if verbose:
        l = f"query name from {oneread.query_name} to {q}."
        if logger is None:
            print(l)
        else:
            logger.info(l)
    oneread.query_name = q
    try:
        outfhandles[b2].write(oneread)
    except KeyError:
        l = f"{b2} is not found in outfhandles."
        if logger is None:
            print(l)
        else:
            logger.info(l)
    finally:
        return None

# * main
barcode_pattern = ':\w\w:\w\w:\w\w:'
cell_meta = pd.read_csv(cell_meta_fnm,
                        sep = ",", header = 0,
                        names = None)
cell_meta = cell_meta[
    cell_meta['barcode'].str.contains(sublib_id, na = False)]
barcode2outdir: Dict[str, str] = OrderedDict(
    [(row['barcode'] , os.path.join(out_bam_dir, "{r}_{h}_{s}".format(
        r = row['brainregion'],
        h = row['modality'],
        s = row['rep']))) for _, row in cell_meta.iterrows()]
)

uoutdirs = list(set(barcode2outdir.values()))
# for d in uoutdirs:
#     os.makedirs(d, exist_ok = True)

barcode2outfnm: Dict[str, str] = OrderedDict(
    [k, os.path.join(v, f"{k}.srt.rmdup.bam")]
    for k, v in barcode2outdir.items()
)

barcodes: List[str] = get_barcodes_from_bam_async(
    bam_fnm = bam_fnm,
    index_fnm = index_fnm,
    barcode_pattern = barcode_pattern,
    nproc = nthread
)
logger.info(f"{len(barcodes)} from {sublib_id}.")
ubarcodes = list(set(barcodes))
logger.info(f"unique {len(ubarcodes)} from {sublib_id}.")

uubarcodes = [f"{sublib_id}:{u}" for u in ubarcodes]
uubarcodes = [u for u in uubarcodes if u in barcode2outfnm]
assert len(uubarcodes) == len(set(uubarcodes))
logger.info(f"{len(uubarcodes)} passed QC.")

if overwrite >= 1:
    logger.info("Overwrite existed bam files.")
else:
    logger.info("Only bam files not found will be generated.")
    uubarcodes = [u for u in uubarcodes
                  if not os.path.exists(barcode2outfnm[u])]
    if len(uubarcodes) < 1:
        logger.info("All the barcodes have the corresponding files.")
        sys.exit(status = 0)
    logger.info(f"{len(uubarcodes)} will be processed.")
    
# get new head
with pysam.AlignmentFile(
    bam_fnm, index_filename = index_fnm, mode = "rb") as abam:
    header = abam.header.to_dict()
    validKeys = [x for x in ["HD", "SQ", "RG"] if x in header.keys()]
    newhead = dict((k, header[k]) for k in validKeys)

# set up file handlers
logger.info(
    f"create {len(uubarcodes)} out bam filehandles for {sublib_id}.")
barcode2outfhs: Dict[str, Any] = OrderedDict(
    [(k,
      pysam.AlignmentFile(barcode2outfnm[k], 'wb', header = newhead))
     for k in uubarcodes]
)
logger.info(
    f"start to split bam for {sublib_id} with {len(uubarcodes)} barcodes.")
with pysam.AlignmentFile(bam_fnm, index_filename = index_fnm,
                         mode = 'rb') as abam:
    # for r in abam.fetch('chr1',0 , 195471971):
    for r in abam.fetch():
        write_oneread(r, sublib_id = sublib_id,
                      barcodes = set(uubarcodes),
                      outfhandles = barcode2outfhs,
                      logger = logger,
                      barcode_pattern = ':\w\w:\w\w:\w\w:',
                      verbose = verbose)
logger.info(f"end of splitting bam for {sublib_id}.")
for _, fh in barcode2outfhs.items():
    fh.close()
logger.info(f"close all the handles in {sublib_id}.")


