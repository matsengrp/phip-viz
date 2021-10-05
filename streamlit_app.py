#!/usr/bin/env python3

# https://share.streamlit.io/1edv/evolution/app/app.py

import copy
import os
import json

import altair as alt

import streamlit as st
import xarray as xr
import pandas as pd
import numpy as np
import dask
import phippery
from phippery.tidy import tidy_ds

# initialize wide view
st.set_page_config(layout='wide')


# 
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
    suppress_st_warning=True,
    max_entries=10
)
def load_data(input_file_path: str, df: pd.DataFrame, **kwargs):
    #st.write("Xarray data Cache miss")

    ds = phippery.load(input_file_path)
    sid, pid = phippery.id_coordinate_from_query(ds, df)
    return ds.loc[dict(
        sample_id=sid,
        peptide_id=pid
    )]

ph = st.sidebar.empty()

# TITLE OF THE APP
# TODO decorate
st.sidebar.title('Working Dataset Summary')
with st.sidebar:
    """
    The core feature of this app is to first select the subset of the dataset 
    which you would like to visualize. Queries are applied using the pandas df
    [query heuristic]()
    """

"""
# PhIP-Seq Interactive enrichment visualizer

### Welcome!!

*To get started, select a dataset file from the options in the sidebar to the left.*

This app is intended for viewing the many facets of an enrichment matrix resulting
from a Phage Immuno Precipitation experiment (PhIP-Seq). 
This visualization app is part of the *phippery suite* of tools, the 
[documentation]() will tell you how to create input from your own data, or obtain 
example data to play with.
"""



# SELECT WHICH FILE TO VIEW
input_file_list = [
    f"{fp}" for fp in os.listdir(".")
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
        ds.sample_table.to_pandas().reset_index().infer_objects()
)
st.session_state['peptide_table'] = copy.deepcopy(
        ds.peptide_table.to_pandas().reset_index().infer_objects()
)

unique_types = ["Heatmap"]
selected_types = st.sidebar.multiselect(
    "Types",
    unique_types,
    default=unique_types
)

sample_expand = True if qtype == 'sample' else False
with st.expander('Sample Table', expanded=sample_expand):
    st.write(st.session_state.sample_table)
    st.write('Juicy deets')

peptide_expand = True if qtype == 'peptide' else False
with st.expander('Peptide Table', expanded=peptide_expand):
    st.write(st.session_state.peptide_table)
    st.write('Juicy deets')

#if "Heatmap" in selected_types and not st.session_state.clicked_query_type_flag:
if "Heatmap" in selected_types: 


    with st.expander('Heatmap Settings', expanded=True):
        left_column, right_column = st.columns(2)
        #with left_column:
        if True:
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

                agg_func_choices = [
                        "average", 
                        "max", 
                        "mean",
                        "median",
                        "min",
                        "stderr", 
                        "stdev", 
                        "sum"
                ]
                agg_func = st.selectbox(
                    "aggregate function",
                    agg_func_choices,
                    index=0
                )

                # How to group the samples
                y_choices = list(ds.sample_metadata.values)
                y = st.selectbox(
                    "y-axis sample feature",
                    ["sample_id"] + y_choices,
                    index=0
                )

                #y = st.text_input(label=f"Y")
                # How to group the peptides
                x_choices = list(ds.peptide_metadata.values)
                x = st.selectbox(
                    "x-axis peptide feature",
                    ["peptide_id"] + x_choices,
                    index=0
                )

                facet_choices = list(ds.peptide_metadata.values)
                facet = st.selectbox(
                    "peptide facet feature",
                    ["None"] + facet_choices,
                    index=0
                )
                domain_max = st.number_input("domain max")
                submitted = st.form_submit_button("Render Heatmap")


        #with right_column:

            if submitted:
                # TODO, we'll want to check the axis they've chosen
                # are unique or throw a warning??
                sm = [y] if y != 'sample_id' else []
                pm = [x] if x != 'peptide_id' else []
                pm = pm + [facet] if facet != 'None' else pm

                # throw out all things we don't care about before 
                # creating the tall dataframe (quite memory expensive)
                subset_ds = copy.deepcopy(ds.loc[
                    dict(
                        sample_metadata = sm,
                        peptide_metadata = pm
                    )
                ])

                keep_tables = set(["sample_table", "peptide_table", enrichment])
                for dt in set(list(subset_ds.data_vars)) - keep_tables:
                    del subset_ds[dt]

                tall_subset = tidy_ds(subset_ds)

                kwargs = {}
                if domain_max:
                    kwargs["scale"] = alt.Scale(domain=[0, domain_max])
                color = alt.Color(f'{agg_func}({enrichment}):Q', **kwargs)

                c = alt.Chart(tall_subset).mark_rect().encode(
                    x=alt.X(f'{x}:O'),
                    y=alt.Y(f'{y}:O'),
                    color=color
                )
               
                c = c.facet(facet, columns=1) if facet != 'None' else c

                st.altair_chart(c, use_container_width=True)

