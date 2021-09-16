#!/usr/bin/env python3

import copy
from plotnine import *
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import altair as alt
import seaborn as sns
import streamlit as st
import os
import pickle
import xarray as xr
import pandas as pd
import dask
import phippery
from phippery.tidy import tidy_ds

st.write("EXECUTED FROM TOP")

@st.cache(allow_output_mutation=True)
def get_sample_queries():
    return []

if 'already_subset_flag' not in st.session_state:
    st.session_state.already_subset_flag = 0 
    st.write('first_exec')
else:
    st.session_state.already_subset_flag += 1
    st.write(f'exec # {st.session_state.already_subset_flag}')
    st.write(f"sample_queries: {get_sample_queries()}")


# TODO Docstrings
@st.cache(hash_funcs={xr.core.dataset.Dataset: dask.base.tokenize}, suppress_st_warning=True)
def load_data(input_file_path: str, **kwargs):

    # no queries exist
    if len(get_sample_queries()) == 0:
        ds = phippery.load(input_file_path)
        return ds
        #return ds.loc[dict(
        #    sample_id=st.session_state.sample_id,
        #    peptide_id=st.session_state.peptide_id
        #)]

    else:
        ds = phippery.load(input_file_path)
        st.write(f"using sample queries: {get_sample_queries()}")
        return ds
        # TODO actually query
        #return ds.loc[dict(
        #    sample_id=st.session_state.sample_id,
        #    peptide_id=st.session_state.peptide_id
        #)]


@st.cache(hash_funcs={xr.core.dataset.Dataset: dask.base.tokenize}, suppress_st_warning=True)
def compute_intervals(
    ds, 
    enrichment,
    window_size,
    aff_func, 
):

    my_bar = st.progress(0)
    enrichments = ds[f"{enrichment}"]
    
    # grab s and p tables
    s_table = ds.sample_table.to_pandas().infer_objects()
    p_table = ds.peptide_table.to_pandas().infer_objects()
    p_table["Loc"] = p_table["Loc"].astype(int) 
    
    # Sum up the enrichments within windows
    windows = list(range(0, max(p_table["Loc"]), window_size))
    enrichment_columns = []
    aggregated_enrichment = {}

    # TODO I'm sure there's a better way to vectorize.
    # but I kinda like the progress bar
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

        percent_complete = int((l / len(windows))*100)
        my_bar.progress(percent_complete)
    my_bar.empty()
    agg_samples = pd.concat(
            [s_table, pd.DataFrame(aggregated_enrichment, index=s_table.index)], 
            axis=1
    )
    
    return (agg_samples, enrichment_columns)

# TITLE OF THE APP
# TODO decorate
st.sidebar.title('PhIPseq Viz')

# TODO
# Long description.
#expander = st.expander("FAQ")
#expander.write("Here you could put in some really, really long explanations...")

unique_types = ["Heatmap"]
selected_types = st.sidebar.multiselect(
    "Types",
    unique_types,
    default=unique_types
)

def update_sample_query(*args, **kwargs):
    get_sample_queries().append(st.session_state[args[0]])


with st.sidebar.expander('Sample Query Condtions'):
    st.write(f"{get_sample_queries()}")

    num_queries = len(get_sample_queries())
    inp = st.text_input(
        label=f"Condition: {num_queries+1}", 
        key=f"sq_{num_queries+1}", 
        on_change=update_sample_query,
        args=tuple([f"sq_{num_queries+1}"])
    )

# SELECT WHICH FILE TO VIEW
input_file_list = [
    f"_data/{fp}" for fp in os.listdir("_data/")
    if fp.endswith('.phip')
]
selected_input_file = st.sidebar.selectbox(
    "Input File",
    input_file_list
)

ds = load_data(selected_input_file)
    
##################################################
# Data
##################################################
if "Heatmap" in selected_types:

    with st.form("dt"):
        st.write(f"Data Transform")

        @st.cache(hash_funcs={xr.core.dataset.Dataset: dask.base.tokenize}, suppress_st_warning=True)
        def get_valid_enr(ds):
            enrichment_options = []
            for dt in set(list(ds.data_vars)) - set(["sample_table", "peptide_table"]):
                if ds[dt].values.flatten().min() != ds[dt].values.flatten().max():
                    enrichment_options.append(dt)
            return enrichment_options


        
        enrichment_options = get_valid_enr(ds)
        enrichment = st.selectbox(
            "Normalization layer",
            enrichment_options
        )
        #to_drop = set(enrichment_options) - set([enrichment])
        
        window_size_options = list(range(1, len(ds.peptide_id.values)))
        window_size = st.selectbox(
            "Peptide Window Size",
            window_size_options
        )
        
        agg_choices = ['sum', 'median', 'mean']
        agg_func = st.selectbox(
            "Window Aggregation Function",
            agg_choices
        )
        
        #print(agg_samples)

####    ##############################################
# Ae    st.sidebar.
####    ##############################################

        agg_samples, enrichment_columns = compute_intervals(
            ds, 
            enrichment,
            window_size,
            agg_func
        )
        submitted = st.form_submit_button("Submit")

    with st.form("p"):
        st.write(f"Plot Aesthetics")

        # How to group the samples
        sample_groupby_choices = ds.sample_metadata.values
        sample_groupby = st.selectbox(
            "Sample Order By",
            sample_groupby_choices,
            index=0
        )

        sample_order = []
        for g, g_df in agg_samples.groupby(
                [sample_groupby]
                ):
            sample_order.extend(g_df.index)    
        agg_samples = agg_samples.reindex(sample_order)

        ## Range of values to viz
        s_e = copy.deepcopy(agg_samples[enrichment_columns])

        a = int(s_e.values.flatten().min())
        b = int(s_e.values.flatten().max())
        values = st.slider(
            "Min Max Enrichment",
            a,
            b,
            (a,b)
        )
        s_e[s_e < values[0]] = values[0]
        s_e[s_e > values[1]] = values[1]

        # Colormap
        heatcolormap_choices = [
                "RdPu",
                "bwr",
                "bone",
                "pink"
        
        ]
        heatcolormap = st.selectbox(
            "Colormap",
            heatcolormap_choices,
            0
        )
        hcmap = getattr(cm, heatcolormap)

        norm = st.selectbox(
            "Color Scale",
            ["Log", "Linear"],
            0
        )

        submitted = st.form_submit_button("Submit")

    fig, ax = plt.subplots(figsize=[10, 14])
    cbar_kws = dict(use_gridspec=False, location="right", label=enrichment)
    sns.heatmap(
            s_e, 
            ax = ax, 
            cmap=hcmap, 
            cbar_kws=cbar_kws, 
            norm=LogNorm() if norm=="Log" else None
    )
    st.write(fig)

    #else:
    #    sns.heatmap(s_e, ax = ax, vmin=values[0], vmax=values[1], cmap=hcmap, cbar_kws=cbar_kws)

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






