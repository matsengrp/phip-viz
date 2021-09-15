#!/usr/bin/env python3

from plotnine import *
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import seaborn as sns
import streamlit as st
import os
import pickle
import xarray as xr
import pandas as pd
import dask
import phippery
from phippery.collapse import collapse_groups


# TITLE OF THE APP
st.sidebar.title('PhIPseq Viz - Sample Enrichment Heatmap')

# SELECT WHICH FILE TO VIEW
input_file_list = [
    f"_data/{fp}" for fp in os.listdir("_data/")
    if fp.endswith('.phip')
]
selected_input_file = st.sidebar.selectbox(
    "Input File",
    input_file_list
)


# this is a messy way to access the data from my pickle dump binary
# phippery includes some niceer wrapper functions for this st.sidebar.ff
# but we'll forget about that for now

#@st.sidebar.cache(suppress_st.sidebar.warning=True)
@st.cache(hash_funcs={xr.core.dataset.Dataset: dask.base.tokenize}, suppress_st_warning=True)
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

#@st.cache(hash_funcs={xr.core.dataset.Dataset: dask.base.tokenize}, suppress_st_warning=True)
#def compute_intervals(
#        ds, 
#        enrichment,
#        agg_func,
#        window_size, 
#        sample_groupby, 
#):
    



def floor_ceil_id(minv, maxv, v):
    if v < minv:
        return minv
    elif v > maxv:
        return maxv
    else:
        return v

# TODO not sure if this should be in a container
c, p = 0,0


#with st.sidebar.sidebar.form("calculate"):
#    st.sidebar.write(f"Calculate: {c}")
#slider_val = st.sidebar.slider("Form slider")
#checkbox_val = st.sidebar.checkbox("Form checkbox")

# Every form must.sidebar.have a submit button.
ds = load_data(selected_input_file, collapse=False)


##################################################
# Data
##################################################
enrichment_options = [dt for dt in set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])]
enrichment = st.sidebar.selectbox(
    "Normalization layer",
    enrichment_options
)


window_size_options = list(range(1, len(ds.peptide_id.values)))
window_size = st.sidebar.selectbox(
    "Peptide Window Size",
    window_size_options
)


agg_choices = ['sum', 'median', 'mean']
agg_func = st.sidebar.selectbox(
    "Window Aggregation Function",
    agg_choices
)


#print(agg_samples)
#submitted = st.sidebar.form_submit_button("Submit")
#if submitted:
#    st.sidebar.write("Done")
#    c+=1

#st.sidebar.write("Outside the form")

#with st.sidebar.form("plot"):
#    st.sidebar.write(f"Plot: {p}")



##################################################
# Aest.sidebar.
##################################################


#s_e[s_e<=values[0]]=values[0]
#s_e[s_e>=values[1]]=values[1]
sample_groupby_choices = ds.sample_metadata.values
#print(sample_groupby_choices)
sample_groupby = st.sidebar.selectbox(
    "Sample Order By",
    sample_groupby_choices,
    index=0
)

enrichments = ds[f"{enrichment}"]

# grab s and p tables
s_table = ds.sample_table.to_pandas().infer_objects()
p_table = ds.peptide_table.to_pandas().infer_objects()
p_table["Loc"] = p_table["Loc"].astype(int) 

# Sum up the enrichments within windows
windows = list(range(0, max(p_table["Loc"]), window_size))
enrichment_columns = []
aggregated_enrichment = {}
for l in range(len(windows)-1):
    start = windows[l]
    end = windows[l+1]
    pep_id = p_table.loc[p_table["Loc"].isin(range(start, end)),:].index
    epitope_enrichment = enrichments.loc[pep_id, :]
    enrichment_columns.append(f"[{start}, {end})")
    if agg_func == "sum":
        agg = epitope_enrichment.sum(axis=0).values
    elif agg_func == "median":
        agg = epitope_enrichment.median(axis=0).values
    elif agg_func == "mean":
        agg = epitope_enrichment.mean(axis=0).values
    aggregated_enrichment[f"[{start}, {end})"] = agg 

agg_samples = pd.concat(
        [s_table, pd.DataFrame(aggregated_enrichment, index=s_table.index)], 
        axis=1
)

sample_order = []
for g, g_df in agg_samples.groupby(
        [sample_groupby]
        ):
    sample_order.extend(g_df.index)    
agg_samples = agg_samples.reindex(sample_order)
s_e = agg_samples[enrichment_columns]

a = int(s_e.values.flatten().min())
b = int(s_e.values.flatten().max())
values = st.sidebar.slider(
    "Min Max Enrichment",
    a,
    b,
    (a,b)
)

heatcolormap_choices = [
        "RdPu",
        "bwr",
        "bone",
        "pink"

]
heatcolormap = st.sidebar.selectbox(
    "Colormap",
    heatcolormap_choices,
    0
)
norm = st.sidebar.selectbox(
    "Color Scale",
    ["Log", "Linear"],
    0
)

hcmap = getattr(cm, heatcolormap)
fig, ax = plt.subplots(figsize=[10, 14])
cbar_kws = dict(use_gridspec=False, location="right", label=enrichment)
#sns.heatmap(s_e, ax = ax, cmap=hcmap, cbar_kws=cbar_kws, norm=LogNorm())
if norm == "Log":
    sns.heatmap(s_e, ax = ax, vmin=values[0], vmax=values[1], cmap=hcmap, cbar_kws=cbar_kws, norm=LogNorm())
else:
    sns.heatmap(s_e, ax = ax, vmin=values[0], vmax=values[1], cmap=hcmap, cbar_kws=cbar_kws)

#base=0
#for g, g_df in s_table.groupby():
#    
#    height = len(g_df.index)
#
#    rect_v = patches.Rectangle(
#            (-10, base),
#            width=7,
#            height=height,
#            clip_on=False, 
#            linewidth=1, 
#            edgecolor='black',
#            facecolor="None"
#    )
#    axd["A"].axhline(base + height, lw=2, color="black")
#    axd["A"].text(-6.0, (base+(base+height))/2, label, rotation=90, va="center", ha="center", size=14)
#    axd["A"].add_patch(rect_v)
#    #axd["A"].set_xlabel("Amino acid position")
#    base = base + height
st.write(fig)






