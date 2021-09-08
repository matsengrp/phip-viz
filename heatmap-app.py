#!/usr/bin/env python3

from plotnine import *
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
import os
import pickle
import xarray as xr
import phippery
from phippery.collapse import collapse_groups


# TITLE OF THE APP
st.title('PhIPseq Viz - Sample Enrichment Heatmap')

# SELECT WHICH FILE TO VIEW
input_file_list = [
    fp for fp in os.listdir(".")
    if fp.endswith('.phip')
]
selected_input_file = st.sidebar.selectbox(
    "Input File",
    input_file_list
)


# this is a messy way to access the data from my pickle dump binary
# phippery includes some niceer wrapper functions for this stuff
# but we'll forget about that for now

@st.cache(suppress_st_warning=True)
@st.cache(hash_funcs={xarray.core.dataset.Dataset: }, suppress_st_warning=True)
def load_data(input_file_path, collapse=False, **kwargs):
    ###
    # load the pickled binary dataset
    ###
    return phippery.load(input_file_path)
    #if collapse == False:
    #    return ds
    #else:
    #    return collapse_groups(ds, **kwargs)

    return sample_table

ds = load_data(selected_input_file, collapse=False)


enrichment_options = [dt for dt in set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])]
enrichment = st.sidebar.selectbox(
    "Enrichment Stat",
    enrichment_options
)
enrichments = ds[f"{enrichment}"]

# 
window_size_options = list(range(1, len(ds.peptide_id.values)))
window_size = st.sidebar.selectbox(
    "Peptide Window Size",
    window_size_options
)

# set a max limit on enrichment so outliers don't obfuscate signal
max_enrichment_value_options = list(range(1, 1000))
max_enrichment_value = st.sidebar.selectbox(
    "Max Enrichment Value",
    max_enrichment_value_options
)

# grab s and p tables
s_table = ds.sample_table.to_pandas().infer_object()
p_table = ds.peptide_table.to_pandas().infer_onjects()
p_table["Loc"] = p_table["Loc"].astype(int) 

# Sum up the enrichments within windows
windows = list(range(0, max(p_table["Loc"]), window_size))
enrichment_columns = []
for l in range(len(windows)-1):
    start = windows[l]
    end = windows[l+1]
    pep_id = p_table.loc[p_table["Loc"].isin(range(start, end)),:].index
    epitope_enrichment = enrichments.loc[pep_id, :]
    enrichment_columns.append(f"[{start}, {end})")
    if ag == "sum":
        agg = epitope_enrichment.sum(axis=0).values
    else:
        agg = epitope_enrichment.mean(axis=0).values
    s_table[f"[{start}, {end})"] = [min(max_enrichment_value, v) for v in agg]

print(s_table)
