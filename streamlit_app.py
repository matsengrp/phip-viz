#!/usr/bin/env python3

# https://share.streamlit.io/1edv/evolution/app/app.py

import copy
import os
import json

import altair as alt
from vega_datasets import data

import streamlit as st
import xarray as xr
import pandas as pd
import numpy as np
import dask
import phippery
from phippery.tidy import tidy_ds

# initialize wide view
st.set_page_config(layout='wide')


if 'query_key_index' not in st.session_state:
    st.session_state.query_key_index = 0 

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

#ph = st.sidebar.empty()

# TITLE OF THE APP
# TODO decorate
st.sidebar.title('Dataset query')
q_help = st.sidebar.button("?", key="q_help")
if q_help:
    st.sidebar.info(f"""
        Using the sidebar you can _select_, or _remove_ subsets of the current 
        working dataset to 
    """
)

#with st.sidebar:
#    """
#    The core feature of this app is to first select the subset of the dataset 
#    which you would like to visualize. Queries are applied using the pandas df
#    [query heuristic]()
#    """
st.title('PhIP-Seq Interactive enrichment visualizer (beta)')
"""
### Welcome!!

*To get started, select a dataset file from the options in the sidebar to the left.*

:point_down: More info!
"""

g_help = st.button("?", key="g_help")
if g_help:
    st.info(f"""
        This app is intended for viewing the many facets of an enrichment matrix resulting
        from a set of Phage Immuno Precipitation (PhIP-Seq) experiments. 
        The options in the app are defined by the sample and peptides tables used to 
        produce the enrichment matrix. 

        This visualization app is part of the *phippery suite* of tools.
        For more information about the data format, the 
        [documentation](https://matsengrp.github.io/phippery/introduction.html) 
        will proivide more insight into how to create input from your own data, or obtain 
        example data to play with.

        **A few notes**

        1. The queries are not save if the page refereshes. We will soon allow users
            to upload/download query tables for convinience

        2. The application was built to be as flexible as possible with respect to the 
            types of phage libraries that may be used as input. That being the case, we
            expect that a range of inputs might not produce sensible visualizations.
        
        3. For questions and feedback: jgallowa (at) fredhutch.org 
            
        4. For bug reports, [git issues]()
        """
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

#def set_input_tables(*args, **kwargs):
#    pass

# SELECT WHICH FILE TO VIEW
input_file_list = [
    f"{fp}" for fp in os.listdir(".")
    if fp.endswith('.phip')
]
selected_input_file = st.sidebar.selectbox(
    "Input File",
    input_file_list
)

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

with st.sidebar:
    """
    ## Overlab & Matsen

    NSF, NIH, HHMI
    
    """

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

#unique_types = ["Heatmap"]
#selected_types = st.sidebar.multiselect(
#    "Types",
#    unique_types,
#    default=unique_types
#)
"""
## Current Working Dataset
"""
ds_help = st.button("?", key="ds_help")
if ds_help:
    st.info(f"""
        Each of the exanders, *Enrichments*, *Sample table*, and *Peptide*
    """
)


sample_expand = True if qtype == 'sample' else False
with st.expander('Sample Table', expanded=sample_expand):
    
    np = len(st.session_state.sample_table)
    f"""
    Total number of samples: {np}
    """

    sample_sum_choices = list(ds.sample_metadata.values)
    sample_sum = st.selectbox(
        "Feature summary",
        ["N/A"]  +sample_sum_choices,
        0
    )
    f"""
    ### {sample_sum} sample feature summary:
    """

    if sample_sum != "N/A": 
        ser = st.session_state.sample_table[sample_sum]
        dt = ser.dtype 
        des = ser.describe()
        if dt == object:
            levels = list(set(ser))
            num_levels = len(levels)
            if num_levels >= 50:
                st.info(f"""
                {sample_sum} is a categorical feature with too many factor levels (n = {num_levels}) to describe, 
                here (probably a desciptive string, and thus, not great plotting feature)
                """)
            elif num_levels < 50 and num_levels >=2:
                
                all_factors = ", ".join(levels)
                st.info(f"""
                {sample_sum} is a categorical feature with the following factor levels to investigate, 

                **{all_factors}**

                Some example query statements:

                > {sample_sum} in ['{levels[0]}', '{levels[1]}', ...]

                > {sample_sum} not in ['{levels[0]}', '{levels[-1]}', ...]

                > {sample_sum} != '{levels[-2]}'

                """)

                f"""
                ### {sample_sum} Factor level counts
                """
                st.write(ser.value_counts())
            else:
                st.info(f"""
                    There's only a single factor level, {levels[0]}, across all samples.
                """
                )

        elif dt == bool:
            st.info(f"""
                Boolean Feature:

                Some example query statements:

                > {sample_sum} == True

                > {sample_sum} == False
            """
            )
            f"""
            ### {sample_sum} Factor level counts
            """
            st.write(ser.value_counts())
        else:
            st.info(f"""
                Numerical Feature:

                Some example query statements:

                > {sample_sum} >= {des[1]}

            """
            )
            f"""
            ### {sample_sum} Summary
            """
            des = ser.describe()
            st.write(des)

    """
    ### raw data
    """

    st.write(st.session_state.sample_table)
    st.write('Juicy deets')


peptide_expand = True if qtype == 'peptide' else False
with st.expander('Peptide Table', expanded=peptide_expand):
    
    np = len(st.session_state.peptide_table)
    f"""
    Total number of peptides: {np}
    """

    sample_sum_choices = list(ds.peptide_metadata.values)
    sample_sum = st.selectbox(
        "Feature summary",
        ["N/A"]  +sample_sum_choices,
        0
    )
    f"""
    ### {sample_sum} sample feature summary:
    """

    if sample_sum != "N/A": 
        ser = st.session_state.peptide_table[sample_sum]
        dt = ser.dtype 
        des = ser.describe()
        if dt == object:
            levels = list(set(ser))
            num_levels = len(levels)
            if num_levels >= 50:
                st.info(f"""
                {sample_sum} is a categorical feature with too many factor levels (n = {num_levels}) to describe, 
                here (probably a desciptive string, and thus, not great plotting feature)
                """)
            elif num_levels < 50 and num_levels >=2:
                
                all_factors = ", ".join(levels)
                st.info(f"""
                {sample_sum} is a categorical feature with the following factor levels to investigate, 

                **{all_factors}**

                Some example query statements:

                > {sample_sum} in ['{levels[0]}', '{levels[1]}', ...]

                > {sample_sum} not in ['{levels[0]}', '{levels[-1]}', ...]

                > {sample_sum} != '{levels[-2]}'

                """)

                f"""
                ### {sample_sum} Factor level counts
                """
                st.write(ser.value_counts())
            else:
                st.info(f"""
                    There's only a single factor level, {levels[0]}, across all samples.
                """
                )

        elif dt == bool:
            st.info(f"""
                Boolean Feature:

                Some example query statements:

                > {sample_sum} == True

                > {sample_sum} == False
            """
            )
            f"""
            ### {sample_sum} Factor level counts
            """
            st.write(ser.value_counts())
        else:
            st.info(f"""
                Numerical Feature:

                Some example query statements:

                > {sample_sum} >= {des[1]}

            """
            )
            f"""
            ### {sample_sum} Summary
            """
            des = ser.describe()
            st.write(des)



    """
    ### raw data
    """

    st.write(st.session_state.peptide_table)
    st.write('Juicy deets')

#"""
### Enrichment Matrix
#"""
#et_help = st.button("et_help")
#if et_help:
#    st.info(f"""
#        below we summarize the samples in the _current working dataset_
#    """
#)
#if "Heatmap" in selected_types and not st.session_state.clicked_query_type_flag:
#if "Heatmap" in selected_types: 

enrichment_expand = True if qtype == 'enrichment' else False
with st.expander('Enrichment Matrix', expanded=enrichment_expand):

    enrichment_options = []
    for dt in set(list(ds.data_vars)) - set(["sample_table", "peptide_table"]):
        if ds[dt].values.flatten().min() != ds[dt].values.flatten().max():
            enrichment_options.append(dt)

    enrichment = st.selectbox(
        "Normalization layer",
        enrichment_options
    )
    
"""
### Heatmap
"""

with st.expander("Heatmap settings"):
    with st.form("dt"):
        st.write(f"Transform Data")
    
    
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
    
        """
        ## IP Observation groups
    
        Select how you would like to aggregate and group IP observation annotation
        groups.
        """
        
        """
        select which feature would you like to plot on the **y-axis**
        """
    
        # How to group the samples
        y_choices = list(ds.sample_metadata.values)
        index=0
        if "patient_status" in y_choices:
            index=y_choices.index("patient_status") + 1
        y = st.selectbox(
            "y-axis sample feature",
            ["sample_id"] + y_choices,
            index=index
        )
    
        #"""
        #select which sample feature would you like to use for a **column** facet
        #"""
    
        #sample_facet_choices = list(ds.sample_metadata.values)
        #sample_facet = st.selectbox(
        #    "Facet feature",
        #    ["None"] + sample_facet_choices,
        #    index=0
        #)
    
        """
        ## Peptide Observation groups
    
        Select how you would like to aggregate and group peptide annotation
        groups.
        """
    
        """
        select which feature would you like to plot on the **x-axis**
        """
    
        x_choices = list(ds.peptide_metadata.values)
        index=0
        if "Protein" in x_choices:
            index=x_choices.index("Protein") + 1
        x = st.selectbox(
            "x-axis peptide feature",
            ["peptide_id"] + x_choices,
            index=index
        )
    
        """
        select which peptide feature would you like to use for a **row** facet
        """
    
        peptide_facet_choices = list(ds.peptide_metadata.values)
        index=0
        if "Virus" in x_choices:
            index=peptide_facet_choices.index("Virus") + 1
        peptide_facet = st.selectbox(
            "Facet feature",
            ["None"] + peptide_facet_choices,
            index=index
        )
    
        domain_max = st.number_input("domain max")
        submitted = st.form_submit_button("Render Heatmap")


#with right_column:

#if submitted:
# TODO, we'll want to check the axis they've chosen
# are unique or throw a warning??
sm = [y] if y != 'sample_id' else []
pm = [x] if x != 'peptide_id' else []

pm = pm + [peptide_facet] if peptide_facet != 'None' else pm
#sm = sm + [sample_facet] if sample_facet != 'None' else sm

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
print(tall_subset)

kwargs = {}
if domain_max:
    kwargs["scale"] = alt.Scale(domain=[0, domain_max])
color = alt.Color(f'{agg_func}({enrichment}):Q', **kwargs)

c = alt.Chart(tall_subset).mark_rect().encode(
    x=alt.X(f'{x}:O'),
    y=alt.Y(f'{y}:O'),
    color=color,
    tooltip = [
        alt.Tooltip(f'{agg_func}({enrichment}):O'),
        alt.Tooltip(f'count({enrichment}):O'),
    ]
).properties(
    width=1000
)

facet_kwargs = {}
if peptide_facet != "None":
    facet_kwargs["row"] = peptide_facet
#if sample_facet != "None":
#    facet_kwargs["column"] = sample_facet
#facet_kwargs["bounds"] = "flush"

if len(facet_kwargs) >= 1:
    c = c.facet(**facet_kwargs)

st.altair_chart(c, use_container_width=False)
