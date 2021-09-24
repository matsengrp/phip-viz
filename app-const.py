#!/usr/bin/env python3

import copy
import os
import json

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import seaborn as sns

import streamlit as st
import xarray as xr
import pandas as pd
import numpy as np
import dask
import phippery
from phippery.tidy import tidy_ds

st.set_page_config(layout='wide')

if 'query_key_index' not in st.session_state:
    st.session_state.query_key_index = 0 

if 'drop_query_key_index' not in st.session_state:
    st.session_state.drop_query_key_index = 0 

if 'view_annotations' not in st.session_state:
    st.session_state.view_samples = False

if 'config' not in st.session_state:
    config = json.load(open("config.json", "r"))
    st.session_state.config = config

if 'queries' not in st.session_state:
    st.session_state['queries'] = pd.DataFrame({
            "qkey": [],
            "Type":[], 
            "Condition": []
        }).set_index("qkey")



@st.cache(
    hash_funcs={xr.core.dataset.Dataset: dask.base.tokenize}, 
    suppress_st_warning=True
)
def load_data(input_file_path: str, df: pd.DataFrame, **kwargs):
    st.write("Cache miss")

    ds = phippery.load(input_file_path)
    sid, pid = phippery.id_coordinate_from_query(ds, df)
    return ds.loc[dict(
        sample_id=sid,
        peptide_id=pid
    )]


@st.cache(
    hash_funcs={xr.core.dataset.Dataset: dask.base.tokenize}, 
    suppress_st_warning=True
)
def compute_intervals(
    ds, 
    enrichment,
    window_size,
    agg_func,
    loc_feature
):

    my_bar = st.progress(0)

    enrichments = ds[f"{enrichment}"]
    # grab s and p tables
    s_table = ds.sample_table.to_pandas().replace(np.nan, "None")
    p_table = ds.peptide_table.to_pandas().replace(np.nan, "None")
    p_table[loc_feature] = p_table[loc_feature].astype(int) 
    
    # Sum up the enrichments within windows
    windows = list(range(min(p_table[loc_feature]), max(p_table[loc_feature]), window_size))
    enrichment_columns = []
    aggregated_enrichment = {}

    # TODO I'm sure there's a better way to vectorize.
    # but I kinda like the progress bar

    # TODO is this really the best way to compute windows??
    for l in range(len(windows)-1):
        start = windows[l]
        end = windows[l+1]
        pep_id = p_table.loc[p_table[loc_feature].isin(range(start, end)),:].index
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

ph = st.sidebar.empty()
# TITLE OF THE APP
# TODO decorate
st.sidebar.title('PhIPseq Viz')

# TODO
# Long description.
#expander = st.expander("FAQ")
#expander.write("Here you could put in some really, really long explanations...")

# SELECT WHICH FILE TO VIEW
input_file_list = [
    f"_data/{fp}" for fp in os.listdir("_data/")
    if fp.endswith('.phip')
]
selected_input_file = st.sidebar.selectbox(
    "Input File",
    input_file_list
)


def add_query_condition(*args, **kwargs):

    st.session_state.query_key_index += 1
    try:
        if args[1] == 'sample':
            l = st.session_state.sample_table.query(st.session_state[args[0]])
        else:
            l = st.session_state.peptide_table.query(st.session_state[args[0]])
        if len(l) == 0:
            raise ValueError(f'Zero-length Array')

    except Exception as e:
        st.warning(f"Error: '{e}' -- Condition not applied")
        return

    st.session_state.queries.loc[args[0]] = [args[1], st.session_state[args[0]]]


def drop_query_condition(*args, **kwargs):

    st.session_state.drop_query_key_index += 1
    to_drop = st.session_state[args[0]]
    existing_keys = st.session_state.queries.index.values
    if to_drop not in existing_keys:
        st.warning(f'{to_drop} does not exist with any of the condition keys, non-operation. Existing Keys include: {existing_keys}')
    else:
        st.session_state.queries.drop(to_drop, axis=0, inplace=True)

qtype = st.sidebar.selectbox(
        'query type', 
        ["sample", "peptide"],
)


num_queries = st.session_state.query_key_index


num_dropped_queries = st.session_state.drop_query_key_index

# Add Query
st.sidebar.text_input(
    label=f"Query Condition",
    key=f"q{num_queries+1}",
    on_change=add_query_condition,
    args=tuple([f"q{num_queries+1}", qtype])
)


st.sidebar.dataframe(st.session_state.queries)

# Remove Query
st.sidebar.text_input(
    label=f"Remove Condition (by key)",
    key=f"rm_key_{num_dropped_queries}",
    on_change=drop_query_condition,
    args=tuple([f"rm_key_{num_dropped_queries}"])
)


# Load data (cached if no change)
df = copy.deepcopy(st.session_state.queries)
ds = load_data(selected_input_file, df)

#if 'sample_table' not in st.session_state:
st.session_state['sample_table'] = copy.deepcopy(
        ds.sample_table.to_pandas().infer_objects()
)
st.session_state['peptide_table'] = copy.deepcopy(
        ds.peptide_table.to_pandas().infer_objects()
)


unique_types = ["Heatmap"]
selected_types = st.sidebar.multiselect(
    "Types",
    unique_types,
    default=unique_types
)


with st.expander('Sample Table', expanded=False):
    st.write(st.session_state.sample_table)
    st.write('Juicy deets')

with st.expander('Peptide Table', expanded=False):
    st.write(st.session_state.peptide_table)
    st.write('Juicy deets')

#if "Heatmap" in selected_types and not st.session_state.clicked_query_type_flag:
if "Heatmap" in selected_types: 


    with st.expander('Heatmap Settings', expanded=True):
        left_column, right_column = st.columns(2)
        with left_column:
            with st.form("dt"):
                st.write(f"Transform Data")

                enrichment_options = []
                for dt in set(list(ds.data_vars)) - set(["sample_table", "peptide_table"]):
                    if ds[dt].values.flatten().min() != ds[dt].values.flatten().max():
                        enrichment_options.append(dt)

                enrichment = st.selectbox(
                    "Normalization layer",
                    enrichment_options
                )
                
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

            
                # TODO Make true checks
                # random thought - this could be something like virus??
                types = st.session_state.peptide_table.dtypes
                int_columns = types[types==np.int64].index.values
                default_loc = st.session_state.config["Locus"]

                if default_loc not in int_columns:
                    st.warning(f"At least one integer peptide feature which specifies the peptide location information is required for PhIP-Viz. Currently, config.json tells us that feature is named '{default_loc}', but does not exist in the peptide table. Be sure to select from the list of valid Loc features below. For more on what defines a locus feature, see: TODO")

                loc_feature = st.selectbox(
                    "Loc",
                    int_columns,
                    index = int(np.where(int_columns==default_loc)[0][0])
                )
                
                agg_samples, enrichment_columns = compute_intervals(
                    ds, 
                    enrichment,
                    window_size,
                    agg_func,
                    loc_feature
                )

                # How to group the samples
                sample_groupby_choices = list(ds.sample_metadata.values)
                sample_groupby = st.selectbox(
                    "Sample Order By",
                    ["Unordered"] + sample_groupby_choices,
                    index=0
                )


                if sample_groupby != "Unordered":
                    sample_order = []
                    for g, g_df in agg_samples.groupby(
                            [sample_groupby]
                            ):
                        sample_order.extend(g_df.index)    
                    assert len(sample_order) == len(agg_samples)
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

        with right_column:

            fig, ax = plt.subplots(figsize=[8, 8])
            cbar_kws = dict(use_gridspec=False, location="right", label=enrichment)
            sns.heatmap(
                    s_e, 
                    ax = ax, 
                    cmap=hcmap, 
                    cbar_kws=cbar_kws, 
                    norm=LogNorm() if norm=="Log" else None
            )

            if sample_groupby != "Unordered":

                base=0
                for g, g_df in agg_samples.groupby([sample_groupby]):
                    
                    height = len(g_df.index)
                
                    #rect_v = patches.Rectangle(
                    #        (-100, base),
                    #        width=70,
                    #        height=height,
                    #        clip_on=False, 
                    #        linewidth=1, 
                    #        edgecolor='black',
                    #        facecolor="None"
                    #)
                    ax.get_yaxis().set_visible(False)
                    ax.axhline(base + height, lw=2, color="black")
                    #ax.text(-1.0*float(window_size), (base+(base+height))/2, g, rotation=90, va="center", ha="right", size=14)
                    ax.text(-1.0, (base+(base+height))/2, g, va="center", ha="right", size=14)
                    #ax.text(-1.0, (base+(base+height))/2, g, rotation=90, va="center", ha="right", size=14)
                    #ax.add_patch(rect_v)
                    #axd["A"].set_xlabel("Amino acid position")
                    base = base + height

            st.write(fig)
