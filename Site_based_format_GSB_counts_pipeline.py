import concurrent.futures
import csv
import heapq
import io
import itertools
import os
import re
import shutil
import sys
import threading
import time

import numpy as np
import pandas as pd

"""
Speed-optimised rewrite - same arguments, same methodology.

Requires csv file (meta.csv) with meta data, sdrf location (eg. "PRIDE/SDRFs/"), gold count and silver count - if meta or sdrf files not available, use "NA"
Alternatively, can use simple meta data files ("simple_meta/") when the meta data is the same for all files in a data set
Can also accept optional decoy and contam prefix as well as modification:target:decoy, if not specified "DECOY" and "CONTAM" prefixes will be used as default
and phospho:STY:A used as search mod.

Environment variables (all optional):
  GSB_WORKERS        number of worker processes (default: available CPUs, capped at the number of tasks)
  GSB_CSV_ENGINE     "pyarrow" (default) or "c" (pandas' reader, as used by the original script)
  GSB_ARROW_THREADS  pyarrow threads per worker process (default 1)


usage: python Site_based_format_GSB_counts_pipeline_optimised.py file_list.txt meta.csv PRIDE/SDRFs/ NA 5 2 [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM] [optional: modification:target:decoy ie. phospho:STY:A]
OR
python Site_based_format_GSB_counts_pipeline_optimised.py file_list.txt NA NA simple_meta/ 5 2 [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM] [optional: modification:target:decoy ie. phospho:STY:A]
OR
python Site_based_format_GSB_counts_pipeline_optimised.py file_list.txt NA NA NA 5 2 [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM] [optional: modification:target:decoy ie. phospho:STY:A]
"""

FLR_FILTER = 0.05

CSV_ENGINE = os.environ.get("GSB_CSV_ENGINE", "pyarrow")
try:
    import pyarrow as pa
    import pyarrow.csv as pa_csv
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False
if CSV_ENGINE == "pyarrow" and not HAVE_ARROW:
    print("WARNING: pyarrow not installed - falling back to pandas' CSV reader and writer (slower)", flush=True)
    CSV_ENGINE = "c"


_BLOCK_SIZE_LOCK = threading.Lock()
LARGE_BLOCK_SIZES = [16 << 20, 64 << 20, 256 << 20, 1024 << 20, (2 << 30) - 1]


def read_csv(path_or_buffer, **kwargs):
    try:
        return pd.read_csv(path_or_buffer, engine=CSV_ENGINE, **kwargs)
    except pd.errors.ParserError as e:
        # pyarrow reads in 1 MB blocks and a single row must fit in one block; files with a very long row (e.g. a
        # peptidoform with thousands of USIs in All_USI) fail with "straddling object straddles two block
        # boundaries". Retry with larger blocks - the block size does not change the values that are read.
        if CSV_ENGINE != "pyarrow" or "straddl" not in str(e):
            raise
        return _read_csv_large_rows(path_or_buffer, **kwargs)


def _read_csv_large_rows(path_or_buffer, **kwargs):
    """pandas does not expose pyarrow's block_size, so pyarrow.csv.ReadOptions (which pandas' pyarrow engine calls)
    is temporarily wrapped to request a larger block size."""
    name = path_or_buffer if isinstance(path_or_buffer, str) else "in-memory table"
    original = pa_csv.ReadOptions
    last_error = None
    for size in LARGE_BLOCK_SIZES:
        if hasattr(path_or_buffer, "seek"):
            path_or_buffer.seek(0)

        def read_options(*args, _size=size, **opts):
            opts.setdefault("block_size", _size)
            return original(*args, **opts)

        with _BLOCK_SIZE_LOCK:
            pa_csv.ReadOptions = read_options
            try:
                df = pd.read_csv(path_or_buffer, engine="pyarrow", **kwargs)
            except pd.errors.ParserError as e:
                if "straddl" not in str(e):
                    raise
                last_error = e
                continue
            finally:
                pa_csv.ReadOptions = original
        print(f"Note: {name} has very long rows - read with pyarrow block size {size >> 20} MB", flush=True)
        return df
    raise last_error


def init_worker():
    """One pyarrow thread per worker process (otherwise every worker starts a thread pool the size of the node)."""
    if HAVE_ARROW:
        n = int(os.environ.get("GSB_ARROW_THREADS", 1))
        pa.set_cpu_count(n)
        pa.set_io_thread_count(n)


def available_cpus():
    try:
        n = len(os.sched_getaffinity(0))  # Linux: respects the SLURM / cgroup CPU allocation
    except AttributeError:
        n = os.cpu_count() or 1
    slurm = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm and slurm.isdigit():
        n = min(n, int(slurm))
    return max(1, n)


# ----------------------------------------------------------------------------------------------------------------
# arguments / file names
# ----------------------------------------------------------------------------------------------------------------

def get_pxd(i):
    PXD_index = [idx for idx, s in enumerate(i.split("/")) if 'PXD' in s][0]
    return i.split("/")[PXD_index]


def parse_args(argv):
    with open(argv[1], "r") as fh:
        folder_list = [f for f in fh.read().replace('\n', ';').split(";") if f.strip() != ""]

    dataset_list = []
    for i in folder_list:
        if get_pxd(i) not in dataset_list:
            dataset_list.append(get_pxd(i))

    cfg = {
        'folder_list': folder_list,
        'dataset_list': dataset_list,
        'meta_all': argv[2],
        'sdrf_dir': argv[3],
        'simple_meta_dir': argv[4],
        'gold_count': int(argv[5]),
        'silver_count': int(argv[6]),
        'dict_all': None,
    }
    if cfg['meta_all'] != "NA":
        with open(cfg['meta_all'], 'r') as f:
            reader = csv.reader(f)
            cfg['dict_all'] = {rows[0]: rows[1:] for rows in reader}

    if len(argv) > 7:
        cfg['decoy_prefix'] = argv[7]
        cfg['contam_prefix'] = argv[8]
    else:
        cfg['decoy_prefix'] = "DECOY"
        cfg['contam_prefix'] = "CONTAM"
    if len(argv) > 9:
        cfg['search_mod'] = argv[9].split(":")[0]
        cfg['target_aa'] = argv[9].split(":")[1]
        cfg['decoy_aa'] = argv[9].split(":")[2]
    else:
        cfg['decoy_aa'] = "A"
        cfg['search_mod'] = "Phospho"
        cfg['target_aa'] = "STY"
    return cfg


def site_file(output_location, i, p):
    return output_location + "/" + i.replace("/", "_") + "_Site_" + p + "_centric.tsv"


def merged_file(output_location, dataset, p):
    return output_location + "/" + dataset + "_merged_Site_" + p + "_centric.tsv"


def all_merged_file(output_location, p):
    return output_location + "/all_datasets_merged_Site_" + p + "_centric.tsv"


def dataset_folders(cfg):
    # same grouping as the original (folder belongs to a dataset if the dataset ID is in its path)
    return {d: [x for x in cfg['folder_list'] if d in x] for d in cfg['dataset_list']}


def input_size(folder, p, decoy_method):
    name = "binomial" if p == "PSM" else "binomial_peptidoform_collapsed_FLR"
    try:
        return os.path.getsize(folder + "/FDR_0.01/" + name + decoy_method + ".csv")
    except OSError:
        return 0


# ----------------------------------------------------------------------------------------------------------------
# writing
# ----------------------------------------------------------------------------------------------------------------

def _finish(tmp, path):
    os.replace(tmp, path)


def write_pandas(df, path):
    tmp = path + ".tmp"
    df.to_csv(tmp, sep="\t", index=False)
    _finish(tmp, path)


def _arrow_column(s):
    if s.dtype == bool:
        return pa.array(np.where(s.to_numpy(), "True", "False"))
    if s.dtype == object:
        try:
            arr = pa.array(s, from_pandas=True)
            if pa.types.is_null(arr.type):   # column with only missing values
                return pa.nulls(len(s), type=pa.string())
            if not pa.types.is_boolean(arr.type):
                return arr
        except (pa.ArrowInvalid, pa.ArrowTypeError):
            pass
        # mixed types: text as pandas would write it
        return pa.array([None if (v is None or (isinstance(v, float) and v != v)) else str(v) for v in s.tolist()],
                        type=pa.string())
    return pa.array(s, from_pandas=True)


WRITE_BLOCK_ROWS = 100_000  # rows converted to pyarrow at a time (keeps the extra memory for writing small)


def write_tsv_fast(df, path):
    """Tab-separated output via pyarrow (values identical to pandas' to_csv, formatting can differ slightly).
    Falls back to pandas if pyarrow is unavailable or a value contains a tab / quote / newline."""
    if not HAVE_ARROW:
        return write_pandas(df, path)
    tmp = path + ".tmp"
    names = [str(c) for c in df.columns]
    options = pa_csv.WriteOptions(include_header=False, delimiter="\t", quoting_style="none")
    try:
        header = io.StringIO()
        csv.writer(header, delimiter="\t", lineterminator="\n").writerow(names)
        with open(tmp, "wb") as fh:
            fh.write(header.getvalue().encode("utf-8"))
            for start in range(0, len(df), WRITE_BLOCK_ROWS):
                block = df.iloc[start:start + WRITE_BLOCK_ROWS]
                table = pa.table([_arrow_column(block[c]) for c in df.columns], names=names)
                pa_csv.write_csv(table, fh, options)
                del table
    except (pa.ArrowInvalid, pa.ArrowTypeError, pa.ArrowNotImplementedError):
        return write_pandas(df, path)
    _finish(tmp, path)


def concat_tsv(paths, out):
    """Concatenate tab-separated files with identical headers (header written once). Falls back to reading and
    re-writing if the headers differ (e.g. files left over from a run of another script version)."""
    headers = []
    for p in paths:
        with open(p, "rb") as fh:
            headers.append(fh.readline().rstrip(b"\r\n"))
    if any(h != headers[0] for h in headers):
        dfs = [read_csv(p, sep="\t") for p in paths]
        return write_tsv_fast(pd.concat(dfs) if len(dfs) > 1 else dfs[0], out)
    tmp = out + ".tmp"
    with open(tmp, "wb") as fo:
        for k, p in enumerate(paths):
            with open(p, "rb") as fi:
                header = fi.readline()
                if k == 0:
                    fo.write(header)
                shutil.copyfileobj(fi, fo, 16 * 1024 * 1024)
            with open(p, "rb") as fi:
                fi.seek(0, os.SEEK_END)
                if fi.tell() > 0:
                    fi.seek(-1, os.SEEK_END)
                    if fi.read(1) != b"\n":
                        fo.write(b"\n")
    _finish(tmp, out)


# ----------------------------------------------------------------------------------------------------------------
# per-folder site tables
# ----------------------------------------------------------------------------------------------------------------

def site_protein_positions(proteins, prot_positions, ptm_positions, ptm_pos):
    """Keep only the protein positions of the site being reported (same string logic as the original loop)."""
    target = int(ptm_pos)
    ptm_split = ptm_positions.split(";")
    protein_PTM = ""
    for z, y in zip(proteins.split(":"), prot_positions.split(":")):
        for a, b in zip(ptm_split, y.split(";")):
            if int(a) == target:
                protein_PTM += b + ";"
        protein_PTM = protein_PTM[:-1] + ":"
    return protein_PTM[:-1]


def remap_protein_positions(df):
    pp = df['All_PTM_protein_positions'].astype(str).tolist()
    cache = {}
    out = []
    for key in zip([str(v) for v in df['All_Proteins'].tolist()], pp,
                   [str(v) for v in df['All_PTM_positions'].tolist()], df['PTM positions'].tolist()):
        v = cache.get(key)
        if v is None:
            v = cache[key] = site_protein_positions(*key)
        out.append(v)
    df['All_PTM_protein_positions'] = out


def map_unique(series, func, dtype):
    """Apply func (list of distinct values -> list of results) once per distinct value and map the results back to
    every row. Returns None if the column has missing values (callers then use the original row-by-row code)."""
    codes, uniques = pd.factorize(series)
    if (codes < 0).any():
        return None
    results = func(list(uniques))
    if dtype is object:
        arr = np.empty(len(results), dtype=object)
        for k, v in enumerate(results):
            arr[k] = v
    else:
        arr = np.array(results, dtype=dtype)
    return arr[codes]


def clean_source(source):
    return str(source).replace("_raw", "").replace("(", "").replace(")", "")


def collapse_metadata(per_file):
    """Single value if all source files agree, otherwise one value per source file, ';'-separated."""
    if all(v == per_file[0] for v in per_file):
        return per_file[0]
    return ";".join(per_file)


PSM_COLS = ['All_Proteins', 'All_PTM_protein_positions', 'All_PTM_positions', 'All_PTMs', 'All_PTM_scores',
            'All_USI', 'Peptide', 'Peptide_mod', 'Protein position', 'PTM Score', 'Score',
            'PTM positions', 'PTM_info', 'Binomial_final_score', 'Source', 'USI']
PEP_COLS = ['All_Proteins', 'All_PTM_protein_positions', 'All_PTM_positions', 'All_PTMs', 'All_PTM_scores',
            'All_USI', 'All_Source', 'Peptide', 'Protein', 'Protein position', 'PTM Score', 'Score', 'PTM positions',
            'PTM_info', 'Binomial_final_score', 'Peptide_mod', 'Peptide_mod_pos', 'USI']
PEP_COUNT_COLS = ['0.05FLR_threshold_count', '0.01FLR_threshold_count',
                  '0.01<P<=0.05_count', '0.05<P<=0.19_count', '0.19<P<=0.81_count', '0.81<P<=0.95_count',
                  '0.95<P<0.99_count', 'P>=0.99_count']


def build_site_table(i, p, decoy_method, cfg):
    """PSM- or Peptidoform-centric site table for one folder (as the original per_folder, before writing)."""
    decoy_aa = cfg['decoy_aa']
    decoy_prefix = cfg['decoy_prefix']
    contam_prefix = cfg['contam_prefix']
    if p == "PSM":
        file = i + "/FDR_0.01/binomial" + decoy_method + ".csv"
        cols = PSM_COLS + ['p' + decoy_aa + '_q_value']
    else:
        file = i + "/FDR_0.01/binomial_peptidoform_collapsed_FLR" + decoy_method + ".csv"
        cols = PEP_COLS + ['p' + decoy_aa + '_q_value_BA'] + PEP_COUNT_COLS
    df = read_csv(file, usecols=cols)[cols].copy()
    df = df.reset_index(drop=True)

    df['Modification'] = cfg['search_mod']
    PXD = get_pxd(i)

    remap_protein_positions(df)
    # Simple meta data mode
    if cfg['simple_meta_dir'] != "NA":
        simple_meta = pd.read_csv(cfg['simple_meta_dir'] + "/" + PXD + ".tsv", sep="\t")
        df['PubMedIDs'] = simple_meta["PubMedIDs"].values[0]
        df['Sample ID'] = simple_meta["Sample ID"].values[0]
        df['Source Dataset Identifier'] = PXD
        df['Reanalysis Dataset Identifier'] = "NA"
        df['Organism'] = simple_meta["Organism"].values[0]
        df['Organism Part'] = simple_meta["Organism Part"].values[0]
        df['Cell Line'] = simple_meta["Cell Line"].values[0]
        df['Disease Information'] = simple_meta["Disease Information"].values[0]
    # SDRF mode
    elif cfg['sdrf_dir'] != "NA":
        meta_all = cfg['meta_all']
        dict_all = cfg['dict_all']
        meta = pd.read_csv(cfg['sdrf_dir'] + "/" + PXD + ".sdrf.tsv", sep="\t")
        meta['comment[data file]'] = meta['comment[data file]'].str.split(".").str[0]
        meta['comment[data file]'] = meta['comment[data file]'].str.replace(" ", "_")
        meta.set_index('comment[data file]', inplace=True)
        meta_dict = meta.to_dict()

        if p == "PSM":
            source_new = pd.Series([clean_source(s) for s in df['Source'].tolist()], index=df.index)
            df['Sample ID'] = source_new.map(meta_dict['source name'])
            df['Organism'] = source_new.map(meta_dict['characteristics[organism]'])
            df['PubMedIDs'] = dict_all[PXD][0] if meta_all != "NA" else "NA"
            df['Source Dataset Identifier'] = PXD
            df['Reanalysis Dataset Identifier'] = "NA"
            try:
                df['Organism Part'] = source_new.map(meta_dict["characteristics[organism part]"])
            except KeyError:
                df['Organism Part'] = "NA"
            try:
                df['Cell Line'] = source_new.map(meta_dict["characteristics[cell line]"])
            except KeyError:
                df['Cell Line'] = "NA"
            try:
                df['Disease Information'] = source_new.map(meta_dict["characteristics[disease]"])
            except KeyError:
                df['Disease Information'] = "NA"
        else:
            # metadata for every source file the peptidoform was identified in (All_Source): per column a single
            # value if all of the peptidoform's source files agree, otherwise one value per file, ';'-separated
            meta_columns = [('Sample ID', "source name", False),
                            ('Organism', "characteristics[organism]", False),
                            ('Organism Part', "characteristics[organism part]", True),
                            ('Cell Line', "characteristics[cell line]", True),
                            ('Disease Information', "characteristics[disease]", True)]
            values = {out_col: [] for out_col, _, _ in meta_columns}
            for source_all in df['All_Source'].tolist():
                sources = [clean_source(s) for s in str(source_all).split(";")]
                for out_col, sdrf_col, optional in meta_columns:
                    if optional and sdrf_col not in meta_dict:
                        per_file = ["NA"] * len(sources)
                    else:
                        per_file = [str(meta_dict[sdrf_col][s]) for s in sources]
                    values[out_col].append(collapse_metadata(per_file))

            df['PubMedIDs'] = dict_all[PXD][0] if meta_all != "NA" else "NA"
            df['Sample ID'] = values['Sample ID']
            df['Source Dataset Identifier'] = PXD
            df['Reanalysis Dataset Identifier'] = "NA"
            df['Organism'] = values['Organism']
            df['Organism Part'] = values['Organism Part']
            df['Cell Line'] = values['Cell Line']
            df['Disease Information'] = values['Disease Information']
    # No meta data mode
    else:
        df['PubMedIDs'] = "NA"
        df['Sample ID'] = "NA"
        df['Source Dataset Identifier'] = PXD
        df['Reanalysis Dataset Identifier'] = "NA"
        df['Organism'] = "NA"
        df['Organism Part'] = "NA"
        df['Cell Line'] = "NA"
        df['Disease Information'] = "NA"

    df['Pep_pos'] = df['Peptide'] + "-" + df['PTM positions'].astype(str)
    q_col = 'p' + decoy_aa + '_q_value' if p == "PSM" else 'p' + decoy_aa + '_q_value_BA'
    df['Site Passes Threshold [0.05]'] = np.where(df[q_col] <= 0.05, 1, 0)
    df['Site Passes Threshold [0.01]'] = np.where(df[q_col] <= 0.01, 1, 0)

    # protein / decoy / contaminant counts, once per distinct protein list
    decoy_re = re.compile(decoy_prefix)    # str.count() in the original treats the prefix as a regex
    contam_re = re.compile(contam_prefix)

    def protein_stats(protein_lists):
        out = []
        for prots in protein_lists:
            lst = prots.split(":")
            dc = sum(pr.startswith(decoy_prefix) or pr.startswith(contam_prefix) for pr in lst)
            out.append((len(lst), len(decoy_re.findall(prots)), len(contam_re.findall(prots)),
                        1 if len(lst) == dc else 0))
        return out

    arr = map_unique(df['All_Proteins'], protein_stats, np.int64) \
        if all(isinstance(v, str) for v in df['All_Proteins'].tolist()) else None
    if arr is not None:
        df['Protein_count'] = arr[:, 0]
        df['Decoy_count'] = arr[:, 1]
        df['Contam_count'] = arr[:, 2]
        df['Decoy / Contaminant Peptide'] = arr[:, 3]
    else:  # original code path
        proteins_list = df['All_Proteins'].str.split(':')
        df['Protein_count'] = proteins_list.str.len()
        df['Decoy_count'] = df['All_Proteins'].str.count(decoy_prefix)
        df['Contam_count'] = df['All_Proteins'].str.count(contam_prefix)
        dcc = np.array([sum(pr.startswith(decoy_prefix) or pr.startswith(contam_prefix) for pr in prots)
                        for prots in proteins_list], dtype=np.int64)
        df['Decoy / Contaminant Peptide'] = np.where(df['Protein_count'] == dcc, 1, 0)

    def decoy_site(x, y):
        return 0 if x[int(y) - 1] != decoy_aa or int(y) == 0 else 1

    def decoy_sites(pep_pos_values):
        return [decoy_site(*pp.split("-")[:2]) for pp in pep_pos_values]

    ds = map_unique(df['Pep_pos'], decoy_sites, np.int64)
    df['Decoy Modification Site'] = np.array(ds if ds is not None else decoy_sites(df['Pep_pos'].tolist()),
                                             dtype=np.int64)
    if p == "PSM":
        df['PSM Site ID'] = df.index
        df['PSM Count Passing Threshold [0.05]'] = "1"
    else:
        df['Peptidoform Site ID'] = df.index

    # peptidoform without reagent labels, once per distinct peptide
    def strip_labels(peptides):
        s = pd.Series(peptides, dtype=object)
        s = s.str.replace(r'\[TMT(.*?)\]', '', regex=True)
        return s.str.replace(r'\[iTRAQ.*?\]', '', regex=True).tolist()

    stripped = map_unique(df['Peptide'], strip_labels, object)
    if stripped is not None:
        df['Peptide'] = stripped
    else:
        df['Peptide'] = df['Peptide'].str.replace(r'\[TMT(.*?)\]', '', regex=True)
        df['Peptide'] = df['Peptide'].str.replace(r'\[iTRAQ.*?\]', '', regex=True)

    # df.replace("not applicable", "NA"), restricted to the non-numeric columns (the only ones it can change).
    # Text columns can be "object" (pandas 2) or the string dtype (pandas 3 / future.infer_string), so only numeric
    # and boolean columns are skipped.
    for c in df.columns:
        if pd.api.types.is_numeric_dtype(df[c]) or pd.api.types.is_bool_dtype(df[c]):
            continue
        mask = df[c].eq("not applicable").fillna(False).to_numpy(dtype=bool)
        if mask.any():
            col = df[c].copy()
            col[mask] = "NA"
            df[c] = col
    return df


# columns of the merged peptidoform tables that the Uniprot format and the GSB step use
UNIPROT_RENAME = {'All_Proteins': 'Proteins',
                  'All_PTM_protein_positions': 'Protein Modification Positions',
                  'Peptide_mod': 'Peptidoform', 'Peptide': 'Unmodified Sequence',
                  'Score': 'PSM Probability',
                  'PTM Score': 'PTM Probability',
                  'PTM positions': 'Peptide Modification Position',
                  'Binomial_final_score': 'Final Probability',
                  'All_PTMs': 'Modifications',
                  'All_PTM_scores': 'Modification probabilities',
                  'All_PTM_positions': 'Modification positions',
                  '0.05FLR_threshold_count': 'PSM Count Passing Threshold [0.05]',
                  '0.01FLR_threshold_count': 'PSM Count Passing Threshold [0.01]',
                  '0.01<P<=0.05_count': 'opt_PSM count 0.01<P<=0.05',
                  '0.05<P<=0.19_count': 'opt_PSM count 0.05<P<=0.19',
                  '0.19<P<=0.81_count': 'opt_PSM count 0.19<P<=0.81',
                  '0.81<P<=0.95_count': 'opt_PSM count 0.81<P<=0.95',
                  '0.95<P<0.99_count': 'opt_PSM count 0.95<P<0.99',
                  'P>=0.99_count': 'opt_PSM count P>=0.99',
                  'USI': 'Universal Spectrum Identifier'}
UNIPROT_COLS = ['Peptidoform Site ID', 'Proteins', 'Unmodified Sequence', 'Peptidoform', 'Modification', 'Peptide Modification Position',
                'Protein Modification Positions', 'PSM Probability', 'PTM Probability', 'Final Probability', 'Site Q-Value', 'Site Passes Threshold [0.05]',
                'Site Passes Threshold [0.01]', 'Decoy / Contaminant Peptide', 'Decoy Modification Site', 'PSM Count Passing Threshold [0.05]',
                'PSM Count Passing Threshold [0.01]', 'Source Dataset Identifier', 'Reanalysis Dataset Identifier', 'PubMedIDs', 'Sample ID', 'Organism',
                'Organism Part', 'Cell Line', 'Disease Information', 'Universal Spectrum Identifier', 'opt_PSM count 0.01<P<=0.05', 'opt_PSM count 0.05<P<=0.19',
                'opt_PSM count 0.19<P<=0.81', 'opt_PSM count 0.81<P<=0.95', 'opt_PSM count 0.95<P<0.99', 'opt_PSM count P>=0.99']
GSB_COLS = ['Peptide_mod_pos', 'Protein', 'Protein position', 'All_Proteins', 'All_PTM_protein_positions', 'Peptide',
            'PTM positions', 'Binomial_final_score', '0.05FLR_threshold_count']


def uniprot_rename(decoy_aa):
    return dict(UNIPROT_RENAME, **{'p' + decoy_aa + '_q_value_BA': 'Site Q-Value'})


def pep_keep_cols(decoy_aa):
    inverse = {v: k for k, v in uniprot_rename(decoy_aa).items()}
    cols = [inverse.get(c, c) for c in UNIPROT_COLS]
    return cols + [c for c in GSB_COLS + ['p' + decoy_aa + '_q_value_BA'] if c not in cols]


def collapse_peptidoforms(df, decoy_aa):
    """Per-dataset merge: keep one row per peptidoform site (lowest q-value, then highest binomial score)."""
    df = df.sort_values(['Peptide_mod_pos', 'p' + decoy_aa + '_q_value_BA', 'Binomial_final_score'],
                        ascending=[True, False, True])
    df = df.drop_duplicates(subset=('Peptide_mod_pos'), keep="last", inplace=False)
    df = df.reset_index(drop=True)
    df['Peptidoform Site ID'] = df.index
    return df


def write_merged_peptidoform(df, path, cfg):
    """Write a dataset's merged peptidoform table with pandas (as the original - its formatting determines how the
    table reads back, which the Uniprot output depends on) and return the table as it reads back from that file."""
    text = df.to_csv(sep="\t", index=False)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)
    _finish(tmp, path)
    return read_csv(io.BytesIO(text.encode("utf-8")), sep="\t", usecols=pep_keep_cols(cfg['decoy_aa']))


def read_merged_peptidoform(path, cfg):
    return read_csv(path, sep="\t", usecols=pep_keep_cols(cfg['decoy_aa']))


# ----------------------------------------------------------------------------------------------------------------
# worker tasks
# ----------------------------------------------------------------------------------------------------------------

def folder_task(i, p, decoy_method, output_location, cfg, only_folder):
    """Site table for one folder. If the folder is its dataset's only folder, the dataset's merged file is written
    directly (no per-folder file). Returns the merged peptidoform table (as read back) in that case."""
    start_time = time.time()
    dataset = get_pxd(i)
    out = merged_file(output_location, dataset, p) if only_folder else site_file(output_location, i, p)
    if os.path.isfile(out):
        print(f"{p}-Centric format for {i} already exists ({os.path.basename(out)})", flush=True)
        if only_folder and p == "Peptidoform":
            return read_merged_peptidoform(out, cfg)
        return None
    print(f"Creating {p}-Centric format for {i}", flush=True)
    df = build_site_table(i, p, decoy_method, cfg)
    result = None
    if p == "PSM":
        write_tsv_fast(df, out)
    elif only_folder:
        result = write_merged_peptidoform(collapse_peptidoforms(df, cfg['decoy_aa']), out, cfg)
    else:
        write_pandas(df, out)
    print(f"{p}-Centric format for {i} complete{' (merged file for ' + dataset + ')' if only_folder else ''}"
          f" --- {time.time() - start_time:.1f} seconds ---", flush=True)
    return result


def dataset_peptidoform_task(dataset, folders, output_location, cfg):
    """Merged peptidoform file for a dataset with several folders (from the per-folder files, as the original)."""
    start_time = time.time()
    out = merged_file(output_location, dataset, "Peptidoform")
    if os.path.isfile(out):
        print(f"Merged Peptidoform files for {dataset} already exists", flush=True)
        return read_merged_peptidoform(out, cfg)
    print(f"creating merged Peptidoform files for {dataset}", flush=True)
    dfs = [read_csv(site_file(output_location, loc, "Peptidoform"), sep="\t",
                    dtype={'All_PTM_protein_positions': str}) for loc in folders]
    df = dfs[0] if len(dfs) == 1 else pd.concat(dfs)
    result = write_merged_peptidoform(collapse_peptidoforms(df, cfg['decoy_aa']), out, cfg)
    print(f"Merged Peptidoform files for {dataset} complete --- {time.time() - start_time:.1f} seconds ---",
          flush=True)
    return result


# ----------------------------------------------------------------------------------------------------------------
# all-dataset outputs
# ----------------------------------------------------------------------------------------------------------------

def finish_psm_files(jobs, groups, only_folder, output_location, dataset_list):
    """Wait for the PSM site files, then build the merged PSM files by concatenation."""
    start_time = time.time()
    for dataset in dataset_list:
        folders = groups[dataset]
        for f in folders:
            jobs[(f, "PSM")].result()
        out = merged_file(output_location, dataset, "PSM")
        if not (len(folders) == 1 and only_folder[folders[0]]):
            if os.path.isfile(out):
                print(f"Merged PSM files for {dataset} already exists", flush=True)
            else:
                concat_tsv([site_file(output_location, f, "PSM") for f in folders], out)
                print(f"Merged PSM files for {dataset} complete", flush=True)
    out = all_merged_file(output_location, "PSM")
    if os.path.isfile(out):
        print("Merged PSM-centric files for all Datasets already exists", flush=True)
    else:
        # original concatenates as pd.concat([new, accumulated]) -> last dataset first
        concat_tsv([merged_file(output_location, d, "PSM") for d in dataset_list[::-1]], out)
        print(f"Merged PSM-centric files for all Datasets complete --- {time.time() - start_time:.1f} seconds"
              f" after the PSM site files were requested ---", flush=True)


def write_all_peptidoform(output_location, dataset_list):
    out = all_merged_file(output_location, "Peptidoform")
    if os.path.isfile(out):
        print("Merged Peptidoform-centric files for all Datasets already exists", flush=True)
        return
    concat_tsv([merged_file(output_location, d, "Peptidoform") for d in dataset_list[::-1]], out)
    print("Merged Peptidoform-centric files for all Datasets complete", flush=True)


def write_uniprot(pep_tables, output_location, cfg):
    """Uniprot format, byte-identical to the original's: the original reads the all-datasets peptidoform file that
    pandas wrote; that read-back is reproduced in memory from the columns the Uniprot format uses."""
    out = output_location + "/all_datasets_merged_Site_Peptidoform_centric_Uniprot.tsv"
    if os.path.isfile(out):
        print("Uniprot format for merged Peptidoform-centric files for all Datasets already exists", flush=True)
        return
    print("Creating Uniprot Peptidoform-centric merged format for all Datasets", flush=True)
    rename = uniprot_rename(cfg['decoy_aa'])
    inverse = {v: k for k, v in rename.items()}
    src_cols = [inverse.get(c, c) for c in UNIPROT_COLS]
    tables = [t[src_cols] for t in pep_tables[::-1]]   # last dataset first, as the original
    df_all = tables[0] if len(tables) == 1 else pd.concat(tables)
    text = df_all.to_csv(sep="\t", index=False)
    df_all = read_csv(io.BytesIO(text.encode("utf-8")), sep="\t")
    df_all = df_all.rename(columns=rename)[UNIPROT_COLS]
    write_pandas(df_all, out)


def uniprot_task(output_location, cfg):
    """Worker task: Uniprot format from the merged peptidoform files (runs alongside the GSB step)."""
    start_time = time.time()
    if os.path.isfile(output_location + "/all_datasets_merged_Site_Peptidoform_centric_Uniprot.tsv"):
        print("Uniprot format for merged Peptidoform-centric files for all Datasets already exists", flush=True)
        return
    tables = [read_merged_peptidoform(merged_file(output_location, d, "Peptidoform"), cfg)
              for d in cfg['dataset_list']]
    write_uniprot(tables, output_location, cfg)
    print(f"Uniprot Peptidoform-centric merged format for all Datasets complete"
          f" --- {time.time() - start_time:.1f} seconds ---", flush=True)


def ratio(peptides, targets, decoy):
    STY_count = 0
    for target in list(targets):
        STY_count += peptides.str.count(target).sum()
    A_count = peptides.str.count(decoy).sum()
    return STY_count / A_count


def init_plotter():
    """Initializer of the plotting process: import matplotlib once, at the start of the run (the import can take
    many seconds on the HPC, where the Python packages are on network storage)."""
    t = time.time()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot  # noqa: F401
    print(f"[timing] plotting process: matplotlib {matplotlib.__version__} imported in {time.time() - t:.1f} s "
          f"(cache directory: {matplotlib.get_cachedir()})", flush=True)


def plotter_ready():
    return True


def plot_categories(counts, path):
    """Bar charts of the G/S/B counts per residue (runs in the plotting process)."""
    t = time.time()
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(ncols=3)
    for k, (name, group) in enumerate(counts.groupby("PTM_FLR_category_FLR", sort=False)):
        axes[k].set_title(name)
        group.plot(kind="bar", x="PTM_residue", y="Count", ax=axes[k], legend=False)
        axes[k].set_ylabel("count")
        axes[k].set_xlabel("")
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close(fig)
    print(f"[timing] plotting process: {os.path.basename(path)} drawn in {time.time() - t:.1f} s", flush=True)


def gsb_counts(pep_tables, output_location, decoy_method, cfg, plotter):
    """GSB counts / mapping files. The charts are drawn by the plotting process; returns their futures."""
    plots = []
    dataset_list = cfg['dataset_list']
    decoy_aa = cfg['decoy_aa']
    gold_count = cfg['gold_count']
    silver_count = cfg['silver_count']
    q_col = 'p' + decoy_aa + '_q_value_BA'
    flr_filter = FLR_FILTER
    prefix = output_location + "/G" + str(gold_count) + "S" + str(silver_count) + "B_" + str(flr_filter)

    print("PSM and Peptidoform formats done, creating GSB counts")
    print("Gold threshold: " + str(gold_count) + "\nSilver threshold: " + str(silver_count))

    # target:decoy ratio over all rows of the all-datasets peptidoform table (as the original)
    STY_ratio = ratio(pd.concat([t['Peptide'] for t in pep_tables]), cfg['target_aa'], decoy_aa)
    print(f"Using target:decoy ratio of: {STY_ratio}")

    prepared = {}
    for dataset, df in zip(dataset_list, pep_tables):
        df = df.loc[df[q_col] <= flr_filter]
        PSM_threshold = dict(df.groupby('Peptide_mod_pos')['0.05FLR_threshold_count'].sum())
        df = df.sort_values(['Peptide_mod_pos', q_col, 'Binomial_final_score'], ascending=[True, False, True])
        df = df.drop_duplicates(subset=('Peptide_mod_pos'), keep="last", inplace=False)
        df['0.05FLR_threshold_count'] = df['Peptide_mod_pos'].map(PSM_threshold)
        if len(df) == 0:
            print(f"WARNING: {dataset} has no peptidoform sites with {q_col} <= {flr_filter} - it contributes no "
                  f"sites to the G/S/B counts (its columns in the mapping files are empty)", flush=True)
        prepared[dataset] = df

    for m in ["single", "all"]:
        print(f"Creating GSB format mapped to {m} proteins")
        counts_list = []
        df_final = None
        for dataset in dataset_list:
            df = prepared[dataset].copy()
            if len(df) >= 1:
                df['PTM_residue'] = [pep[pos - 1] for pep, pos in zip(df['Peptide'].tolist(), df['PTM positions'].tolist())]
                if m == "all":
                    df['All_PTM_protein_positions'] = df['All_PTM_protein_positions'].astype(str)
                    df['All_PTM_protein_positions'] = df['All_PTM_protein_positions'].str.split(':')
                    df['All_Proteins'] = df['All_Proteins'].str.split(':')
                    df = df.explode(['All_PTM_protein_positions', 'All_Proteins'])
                    df["Protein-pos"] = df['All_Proteins'] + "-" + df['All_PTM_protein_positions']
                else:
                    df['Protein-pos'] = df['Protein'] + "-" + df['Protein position'].astype(str)
            else:
                # no sites at 5% in this dataset (the original crashed here with KeyError: 'Protein-pos')
                df['PTM_residue'] = pd.Series(dtype=object)
                df['Protein-pos'] = pd.Series(dtype=object)
            df['PTM_residue'] = np.where(df['PTM positions'] == 0, "N-term", df['PTM_residue'])
            df = df.sort_values(['Protein-pos', q_col, 'Binomial_final_score'], ascending=[True, False, True])
            df['Protein_pos_res'] = df['Protein-pos'] + "_" + df['PTM_residue']
            PSM_threshold_2 = dict(df.groupby('Protein_pos_res')['0.05FLR_threshold_count'].sum())
            df = df.drop_duplicates(subset=('Protein_pos_res'), keep="last", inplace=False)
            df['Sum_of_PSM_counts(5%FLR)'] = df['Protein_pos_res'].map(PSM_threshold_2)

            df['PXD'] = dataset
            counts_list.append(df[['Protein_pos_res', 'PXD']])

            df_temp = df[['Peptide_mod_pos', q_col, 'Binomial_final_score', 'Protein_pos_res',
                          'Sum_of_PSM_counts(5%FLR)', '0.05FLR_threshold_count']]
            df_temp = df_temp.rename(columns={q_col: dataset + "_FLR",
                                              'Binomial_final_score': dataset + '_BinomialScore',
                                              'Peptide_mod_pos': dataset + "_peptide_mod_pos",
                                              'Sum_of_PSM_counts(5%FLR)': dataset + "_Sum_of_PSM_counts(5%FLR)",
                                              '0.05FLR_threshold_count': dataset + "_peptidoform_PSMcount(5%FLR)"})
            df_temp = df_temp.set_index('Protein_pos_res')
            # pairwise concat kept so the row order matches the original exactly
            df_final = df_temp if df_final is None else pd.concat([df_final, df_temp], axis=1)

        column_list = [PXD + "_Sum_of_PSM_counts(5%FLR)" for PXD in dataset_list]
        df_final["Sum_of_PSM_counts(5%FLR)"] = df_final[column_list].sum(axis=1)
        cols = df_final.columns.tolist()
        df_final = df_final[[cols[-1]] + cols[:-1]]
        df_final = df_final.drop(column_list, axis=1)

        df_final = df_final[~df_final.index.str.contains(cfg['decoy_prefix'])]
        df_final = df_final[~df_final.index.str.contains(cfg['contam_prefix'])]

        df_counts = pd.concat(counts_list).reset_index(drop=True)
        df_counts['PTM_residue'] = df_counts['Protein_pos_res'].str.rsplit("_", n=1).str[-1]
        counts_res = pd.crosstab(df_counts['PTM_residue'], df_counts['PXD']).replace(0, np.nan).stack().reset_index().rename(
            columns={0: 'Count'})
        print(counts_res)
        counts_res.to_csv(prefix + "_Residue_counts_" + m + ".csv", index=False)

        # GSB classification: number of datasets in which the site has FLR <= 1%
        index_values = df_final.index.tolist()
        df_final['Protein'] = [s.rsplit("-", 1)[0] for s in index_values]
        df_final['Protein_pos'] = [s.rsplit("-", 1)[1].split("_")[0] for s in index_values]
        df_final['PTM_residue'] = [s.rsplit("_", 1)[-1] for s in index_values]
        flr1_count = np.zeros(len(df_final), dtype=np.int64)
        for dataset in dataset_list:
            flr = pd.to_numeric(df_final[dataset + "_FLR"], errors="coerce")
            flr1_count += (flr <= 0.01).to_numpy(dtype=np.int64)
        df_final['PTM_FLR_category'] = np.where(flr1_count >= gold_count, "Gold",
                                                np.where(flr1_count >= silver_count, "Silver", "Bronze"))

        cols = [c for c in df_final.columns if c not in ('Protein', 'Protein_pos', 'PTM_residue', 'PTM_FLR_category')]
        df_final['Decoy_mod'] = np.where(df_final['PTM_residue'] == decoy_aa, 1, 0)
        df_final = df_final[['Protein', 'Protein_pos', 'PTM_residue', 'Decoy_mod', 'PTM_FLR_category'] + cols]

        # replace 0 PSM counts -> 1
        column_list = [PXD + "_peptidoform_PSMcount(5%FLR)" for PXD in dataset_list]
        column_list.append("Sum_of_PSM_counts(5%FLR)")
        df_final[column_list] = df_final[column_list].replace(0, 1)

        df_final.to_csv(prefix + "_protein_pos" + decoy_method + "_" + m + "_prot_mapping.csv", index=False)

        counts = pd.crosstab(df_final.PTM_residue, df_final.PTM_FLR_category).replace(0, np.nan).stack().reset_index().rename(
            columns={0: 'Count'})

        # FLR estimate per category using the target:decoy ratio
        f = lambda x, y, z: np.nan if x != decoy_aa else str(
            round(z * STY_ratio / counts.loc[counts['PTM_FLR_category'] == y, 'Count'].sum() * 100, 2)) + "%"
        counts['FLR'] = counts.apply(lambda x: f(x.PTM_residue, x.PTM_FLR_category, x.Count), axis=1)
        counts['FLR'] = counts.groupby('PTM_FLR_category')['FLR'].ffill()
        counts['PTM_FLR_category_FLR'] = counts['PTM_FLR_category'] + " - " + counts['FLR']
        counts['PTM_FLR_category'] = pd.Categorical(counts['PTM_FLR_category'], ["Bronze", "Silver", "Gold"])
        counts = counts.sort_values("PTM_FLR_category")
        print(counts, flush=True)

        plots.append(plotter.submit(plot_categories, counts,
                                    prefix + "_protein_pos_categories" + decoy_method + "_" + m + "_prot_mapping.png"))
    return plots


# ----------------------------------------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------------------------------------

def check_inputs(folder_list, decoy_methods):
    checked = [folder_list[0] + "/FDR_0.01/binomial" + dm + ".csv" for dm in decoy_methods]
    if not any(os.path.isfile(f) for f in checked):
        print("ERROR: no input found - none of these files exist:", file=sys.stderr)
        for f in checked:
            print(f"   {f!r}", file=sys.stderr)
        print(f"Current directory: {os.getcwd()}", file=sys.stderr)
        print("Paths in the file list are relative to the current directory. A '\\r' in the paths above means "
              "the file list has Windows line endings (fix with: dos2unix <file_list>).", file=sys.stderr)
        sys.exit(1)
    missing = [fo for fo in folder_list if not any(os.path.isfile(fo + "/FDR_0.01/binomial" + dm + ".csv")
                                                   for dm in decoy_methods)]
    if missing:
        print(f"ERROR: input files not found for {len(missing)} folder(s) in the file list: {missing!r}",
              file=sys.stderr)
        print(f"Current directory: {os.getcwd()}", file=sys.stderr)
        sys.exit(1)


class PriorityScheduler:
    """Hands the pool at most n_workers tasks at a time, always the waiting task with the lowest priority value.
    The pool itself runs tasks in submission order, so without this a dataset merge queued once its folders are done
    would wait behind every folder task still queued. submit() returns a Future that resolves with the task."""

    def __init__(self, pool, n_workers):
        self.pool, self.n_workers = pool, n_workers
        self.cond = threading.Condition()
        self.heap = []
        self.seq = itertools.count()   # ties keep submission order
        self.running = 0
        self.closed = False
        self.thread = threading.Thread(target=self._dispatch, daemon=True)
        self.thread.start()

    def submit(self, priority, fn, *args):
        fut = concurrent.futures.Future()
        with self.cond:
            heapq.heappush(self.heap, (priority, next(self.seq), fn, args, fut))
            self.cond.notify()
        return fut

    def close(self):
        """Stops dispatching; tasks still waiting fail with CancelledError."""
        with self.cond:
            self.closed = True
            waiting, self.heap = self.heap, []
            self.cond.notify()
        for *_, fut in waiting:
            fut.set_exception(concurrent.futures.CancelledError())

    def _finished(self, real, fut):
        with self.cond:
            self.running -= 1
            self.cond.notify()
        if real.cancelled():
            fut.set_exception(concurrent.futures.CancelledError())
        elif real.exception() is not None:
            fut.set_exception(real.exception())
        else:
            fut.set_result(real.result())

    def _dispatch(self):
        while True:
            with self.cond:
                while not self.closed and (not self.heap or self.running >= self.n_workers):
                    self.cond.wait()
                if self.closed:
                    return
                _, _, fn, args, fut = heapq.heappop(self.heap)
                self.running += 1
            try:
                real = self.pool.submit(fn, *args)
            except BaseException as e:   # e.g. BrokenProcessPool: fail this task, keep draining the rest
                with self.cond:
                    self.running -= 1
                fut.set_exception(e)
                continue
            real.add_done_callback(lambda r, fut=fut: self._finished(r, fut))


def run_decoy_method(decoy_method, cfg, n_cpus):
    folder_list = cfg['folder_list']
    dataset_list = cfg['dataset_list']
    output_location = "All_site_formats_Optimised" + decoy_method
    if not os.path.exists(output_location):
        os.mkdir(output_location)

    groups = dataset_folders(cfg)
    only_folder = {f: groups[get_pxd(f)] == [f] for f in folder_list}
    tasks = [(f, p) for p in ["PSM", "Peptidoform"] for f in folder_list]
    tasks.sort(key=lambda t: input_size(t[0], t[1], decoy_method), reverse=True)   # largest first
    n_workers = int(os.environ.get("GSB_WORKERS", 0)) or min(len(tasks), n_cpus)
    print(f"Using {n_workers} worker processes ({n_cpus} CPUs available)", flush=True)
    multi = [d for d in dataset_list if not (len(groups[d]) == 1 and only_folder[groups[d][0]])]
    if multi:
        print(f"Datasets with several folders (per-folder and merged files): {', '.join(multi)}", flush=True)

    t0 = time.time()
    # ProcessPoolExecutor (not multiprocessing.Pool): if a worker dies, e.g. killed for exceeding the job's memory
    # limit, the run stops with an error instead of waiting forever for the lost task
    pool = concurrent.futures.ProcessPoolExecutor(max_workers=n_workers, initializer=init_worker)
    # separate single process for the charts: imports matplotlib while the site files are being built
    plotter = concurrent.futures.ProcessPoolExecutor(max_workers=1, initializer=init_plotter)
    try:
        run_tasks(pool, n_workers, plotter, tasks, groups, only_folder, decoy_method, output_location, cfg, t0)
    except concurrent.futures.process.BrokenProcessPool:
        pool.shutdown(wait=False, cancel_futures=True)
        plotter.shutdown(wait=False, cancel_futures=True)
        print("ERROR: a worker process was killed. The usual cause on the HPC is the job running out of memory - "
              "request more memory (#SBATCH --mem) or use fewer workers (export GSB_WORKERS=...). Completed "
              "outputs are kept; re-running skips them.", file=sys.stderr, flush=True)
        raise
    except BaseException:
        pool.shutdown(wait=False, cancel_futures=True)
        plotter.shutdown(wait=False, cancel_futures=True)
        raise
    pool.shutdown(wait=True)
    plotter.shutdown(wait=True)


def run_tasks(pool, n_workers, plotter, tasks, groups, only_folder, decoy_method, output_location, cfg, t0):
    sched = PriorityScheduler(pool, n_workers)
    try:
        _run_tasks(sched, plotter, tasks, groups, only_folder, decoy_method, output_location, cfg, t0)
    finally:
        sched.close()


# priorities (lower runs first): Uniprot / dataset merges, then folder tasks; largest input first within each
PRIORITY_MERGE, PRIORITY_FOLDER = 0, 1


def _run_tasks(sched, plotter, tasks, groups, only_folder, decoy_method, output_location, cfg, t0):
    dataset_list = cfg['dataset_list']
    errors = []
    # tasks are already sorted largest input first; the sequence number keeps that order within the priority
    jobs = {(f, p): sched.submit((PRIORITY_FOLDER, 0), folder_task, f, p, decoy_method, output_location, cfg,
                                 only_folder[f])
            for f, p in tasks}
    plotter_up = plotter.submit(plotter_ready)   # starts the plotting process (matplotlib import) now

    # the PSM merges (file concatenation) run in a thread alongside
    def psm_worker():
        try:
            finish_psm_files(jobs, groups, only_folder, output_location, dataset_list)
            print(f"[timing] PSM-centric files (per-folder + merged) done at {time.time() - t0:.1f} s",
                  flush=True)
        except BaseException as e:  # re-raised in the main thread
            errors.append(e)

    # daemon threads: if the main thread fails, the pool is terminated and these must not keep the process alive
    psm_thread = threading.Thread(target=psm_worker, daemon=True)
    psm_thread.start()

    # peptidoform chain: per-dataset merged tables -> all-datasets file + Uniprot -> GSB
    # each dataset's merge is queued as soon as its own folders are done (whatever the order of the file list), ahead
    # of folder tasks still waiting; larger datasets first if several become ready together
    merge_jobs = {}
    pep_tables = {}
    pending = {d: set(groups[d]) for d in dataset_list}

    def dataset_ready(d, last_result):
        folders = groups[d]
        if len(folders) == 1 and only_folder[folders[0]]:
            pep_tables[d] = last_result
        else:
            size = sum(input_size(f, "Peptidoform", decoy_method) for f in folders)
            merge_jobs[d] = sched.submit((PRIORITY_MERGE, -size), dataset_peptidoform_task, d, folders,
                                         output_location, cfg)

    for d in dataset_list:
        if not pending[d]:   # no folders in the list for this dataset (as the original: merge of nothing)
            dataset_ready(d, None)
    pep_futures = {jobs[(f, "Peptidoform")]: f for f in cfg['folder_list']}
    for fut in concurrent.futures.as_completed(pep_futures):
        f = pep_futures[fut]
        result = fut.result()
        for d in dataset_list:
            if f in pending[d]:
                pending[d].discard(f)
                if not pending[d]:
                    dataset_ready(d, result)
    for d in dataset_list:
        job = merge_jobs.get(d)
        if job is None:
            continue
        pep_tables[d] = job.result()
    pep_list = [pep_tables[d] for d in dataset_list]
    print(f"[timing] merged Peptidoform-centric files per dataset done at {time.time() - t0:.1f} s", flush=True)

    def all_pep_worker():
        try:
            write_all_peptidoform(output_location, dataset_list)
        except BaseException as e:
            errors.append(e)

    all_pep_thread = threading.Thread(target=all_pep_worker, daemon=True)
    all_pep_thread.start()
    uniprot_job = sched.submit((PRIORITY_MERGE, 0), uniprot_task, output_location, cfg)

    t = time.time()
    plots = gsb_counts(pep_list, output_location, decoy_method, cfg, plotter)
    print(f"[timing] GSB counts: {time.time() - t:.1f} s (done at {time.time() - t0:.1f} s)", flush=True)
    plotter_up.result()
    for job in plots:
        job.result()
    print(f"[timing] GSB charts done at {time.time() - t0:.1f} s", flush=True)

    uniprot_job.result()
    print(f"[timing] Uniprot format done at {time.time() - t0:.1f} s", flush=True)
    all_pep_thread.join()
    psm_thread.join()
    if errors:
        raise errors[0]


def main():
    t0 = time.time()
    cfg = parse_args(sys.argv)
    print(f"Using decoy prefix: {cfg['decoy_prefix']}")
    print(f"Using contam prefix: {cfg['contam_prefix']}")
    print(f"Using modification: {cfg['search_mod']}")
    print(f"Using target amino acids: {cfg['target_aa']}")
    print(f"Using decoy amino acid: {cfg['decoy_aa']}")
    print(f"CSV reader engine: {CSV_ENGINE}; pyarrow writer: {'yes' if HAVE_ARROW else 'no'}", flush=True)

    decoy_methods = ["", "_peptidoform_decoy", "_site_decoy"]
    check_inputs(cfg['folder_list'], decoy_methods)
    n_cpus = available_cpus()
    for decoy_method in decoy_methods:
        if os.path.isfile(cfg['folder_list'][0] + "/FDR_0.01/" + "binomial" + decoy_method + ".csv"):
            run_decoy_method(decoy_method, cfg, n_cpus)
    print(f"[timing] TOTAL: {time.time() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    main()
