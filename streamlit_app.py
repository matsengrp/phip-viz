#!/usr/bin/env python3

# https://share.streamlit.io/1edv/evolution/app/app.py

import copy
import os
import json
import io

import altair as alt
#alt.renderers.enable('altair_saver', fmts=['png', 'pdf'])
from altair_saver import save
alt.data_transformers.disable_max_rows()
#from vega_datasets import data
import phippery
from phippery.utils import *

import streamlit as st
import xarray as xr
import pandas as pd
import numpy as np
import dask
#from phippery.tidy import tidy_ds
#from phippery.string import string_feature
#from phippery.phipdata import get_annotation_table

# initialize wide view
st.set_page_config(layout='wide')
#st.write(st.__version__)

# initialize session state variables
if 'query_key_index' not in st.session_state:
    st.session_state.query_key_index = 0 

if 'drop_query_key_index' not in st.session_state:
    st.session_state.drop_query_key_index = 0 

if 'view_annotations' not in st.session_state:
    st.session_state.view_samples = False

#if 'sample_ex_switch' not in st.session_state:
#    st.session_state.sample_ex_switch = False
#
#if 'peptide_ex_switch' not in st.session_state:
#    st.session_state.peptide_ex_switch = False

req_feats = ["qkey", "expression", "dimension"]
if 'queries' not in st.session_state:
    st.session_state['queries'] = pd.DataFrame({
            feat : [] 
            for feat in req_feats
        }).set_index("qkey")

def id_coordinate_from_query(ds, query_df):

    """
    Given a df with columns 'dimension' and 
    """

    # st = ds.sample_table.to_pandas().infer_objects()
    sq = list(query_df.loc[query_df["dimension"] == "sample", "expression"].values)
    sid = sample_id_coordinate_from_query(ds, sq)

    # pt = ds.peptide_table.to_pandas().infer_objects()
    pq = list(query_df.loc[query_df["dimension"] == "peptide", "expression"].values)
    pid = peptide_id_coordinate_from_query(ds, pq)

    return sid, pid


#@st.cache(
#    hash_funcs={
#        xr.core.dataset.Dataset: dask.base.tokenize,
#    }, 
#    suppress_st_warning=True,
#    max_entries=10
#)
def load_data(input_file_path: str, df: pd.DataFrame, **kwargs):

    ds = load(input_file_path)
    sid, pid = id_coordinate_from_query(ds, df)
    return ds.loc[dict(
        sample_id=sid,
        peptide_id=pid
    )]

@st.cache
def convert_df(df):
    return df.to_csv(index=False).encode('utf-8')


def infer_dim(feature):
    if feature in st.session_state.sample_table.columns:
        return 'sample'
    elif feature in st.session_state.peptide_table.columns:
        return 'peptide'
    else:
        raise ValueError(f"{feature} not in either sample or peptide features")

def get_reasonable_features(df):
    reasonable_facet = []
    for col, data in df.items():
        l = len(set(data.values))
        if l < 50:
            reasonable_facet.append(col)
    return reasonable_facet

#with st.sidebar:
#    """
#    # Upload
#    """
#q_help = st.sidebar.button("?", key="q_help")
#if q_help:
#    st.sidebar.info(f"""
#        Looking at too much data?
#        Overlaping axis groups?
#        Using the sidebar you can _select_, or _remove_ subsets of the current 
#        working dataset.
#        
#        Use the widgets below to apply a condition which subsets the entire dataset.
#    """
#)

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
        The general idea is that you choose the subset of the dataset you're interested in,
        (using the query tools in the sidebar) and then decide how you would like to 
        aggregate (or not) the enrichments of annotation groups as determined by the 
        sample or peptide tables. Use the drop-down menus below to see sub setting options
        for any feature of interest.
        This visualization app is part of the *phippery suite* of tools.
        For more information about the data format, the 
        [documentation](https://matsengrp.github.io/phippery/introduction.html) 
        will provide more insight into how to create input from your own data, or obtain 
        example data to play with.
        **A few notes**
        1. The queries are not save if the page refreshes. We will soon allow users
            to upload/download query tables for convenience
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

    qkey = "q" + args[0].split("-")[1]
    #print("Adding: ", qkey)
    st.session_state.queries.loc[qkey] = [st.session_state[args[0]], args[1]]


def drop_query_condition(*args, **kwargs):

    st.session_state.drop_query_key_index += 1
    to_drop = st.session_state[args[0]]
    existing_keys = st.session_state.queries.index.values
    #if to_drop not in :
    if to_drop not in existing_keys:
        st.warning(f'{to_drop} does not exist with any of the condition keys, non-operation. Existing Keys include: {existing_keys}')
    else:
        st.session_state.queries.drop(to_drop, axis=0, inplace=True)

with st.sidebar:
    """
    # Upload
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

    num_queries = st.session_state.query_key_index
    num_dropped_queries = st.session_state.drop_query_key_index

    uploaded_query = st.file_uploader("Upload Query Table")
    if uploaded_query is not None:
        queries = pd.read_csv(
            uploaded_query,
            header = 0
        )
        # TODO
        #print(f"POST-FILEINPUT: {st.session_state.query_key_index}")
        ind = [f"q{i}" for i in queries.index.values]
        queries["qkey"] = ind
        queries.set_index("qkey", inplace=True)
        st.session_state.queries = queries
        st.session_state.query_key_index = len(queries) -1
        st.info("""
            Upload complete!
            Note: If you would like to edit the queries table you just loaded in,
            be sure to hit the 'x' next to the filename above before adding or
            removing individual queries.
        """)
        #print(f"POST-FILE LEN CHANGE: {st.session_state.query_key_index}")
        #uploaded_query=None

    df = copy.deepcopy(st.session_state.queries)
    #st.text(type(df))
    ds = load_data(selected_input_file, df)
    if len(ds.sample_id.values) == 0:
        raise ValueError(f'Condition file resulted in Zero-length sample table')
    if len(ds.peptide_id.values) == 0:
        raise ValueError(f'Condition file resulted in Zero-length peptide table')

    st.session_state['sample_table'] = get_annotation_table(ds)
    st.session_state['peptide_table'] = get_annotation_table(ds, dim='peptide')

    """
    # Download
    """

    # Load data (cached if no change)


    csv = convert_df(st.session_state.queries)
    st.download_button(
        label="Download queries as CSV",
        data=csv,
        file_name='queries.csv',
        mime='text/csv',
    )


    """
    ## Overlab & Matsen
    NSF, NIH, HHMI
    _Note: ^ placeholder_ 
    """


"""
***************************************
## :mag:    View & Select Working Dataset
"""
ds_help = st.button("?", key="ds_help")
if ds_help:
    st.info(f"""
        Each of the *Sample table*, and *Peptide table*), we provide a 
        summary of the data resulting form the _current working dataset_.
        This means you are provided with summary of the data '_post_-sub setting'
        each time you subset the data using the sidebar tools to the left,
    """
)

#"""
#
#"""


#sample_expand = True if qtype == 'sample' else False
left_s, right_s = st.columns(2)
#with st.expander('Sample Table', expanded=False):
with left_s:
    
    np = len(st.session_state.sample_table)
    f"""
    ### Sample Table
    Total number of samples: {np}
    """
    #def switch_s_expander():
    #    st.session_state.sample_ex_switch = not st.session_state.sample_ex_switch

    #sample_ex = st.expander("sample_expander",
    #        expanded = st.session_state.sample_ex_switch
    #)
    #clicked = sample_ex.button("View & Edit Samples", on_click=switch_s_expander)
    #with sample_ex:
    with st.expander("Working Samples"):

        buffer = io.StringIO()
        st.session_state.sample_table.info(buf=buffer, verbose=True)
        s = buffer.getvalue()
        st.text(s)

        s_q = st.session_state.queries
        st.dataframe(s_q[s_q["dimension"]=="sample"].drop("dimension",axis=1, inplace=False))

        """
        **Tell me more about** :point_down:
        """

        sample_sum_choices = list(ds.sample_metadata.values)
        sample_sum = st.selectbox(
            "Feature summary for:",
            sample_sum_choices,
            0
        )

        s_dtype = st.session_state.sample_table.dtypes[sample_sum]
        s_uniq = len(set(st.session_state.sample_table[sample_sum]))
        s_kwargs={}
        if s_dtype in [pd.Int64Dtype(), pd.Float64Dtype()] and s_uniq < 50:
            view = st.radio(
                    "Sample feature display options:", 
                    ["Quantiles", "Value Counts"])
            s_kwargs["numeric_dis"] = (view=="Quantiles")
        des = string_feature(
            ds, 
            feature= sample_sum, 
            verbosity = 0, 
            dim="sample",
            **s_kwargs
        )
        st.text(des)

        """
        **Apply a query condition to current sample table** :point_down:
        """
        # Add Query
        st.text_input(
            label=f"Sample Query Condition",
            key=f"sq-{num_queries+1}",
            on_change=add_query_condition,
            args=tuple([f"sq-{num_queries+1}", "sample"])
        )


        """
        **Remove (using row index key) a query condition from the current sample table** :point_down:
        """
        
        # Remove Query
        st.text_input(
            label=f"Remove Sample Condition (by key)",
            key=f"rm_s_key_{num_dropped_queries}",
            on_change=drop_query_condition,
            args=tuple([f"rm_s_key_{num_dropped_queries}"])
        )

        """
        **Click below to see the raw sample annotation table** :point_down:
        """
        st.write(st.session_state.sample_table)


####################################

with right_s:
    
    np = len(st.session_state.peptide_table)
    f"""
    ### Peptide Table
    Total number of peptides: {np}
    """
    with st.expander("Working Peptides"):

        buffer = io.StringIO()
        st.session_state.peptide_table.info(buf=buffer, verbose=True)
        s = buffer.getvalue()
        st.text(s)

        p_q = st.session_state.queries
        st.dataframe(p_q[p_q["dimension"]=="peptide"].drop("dimension",axis=1,inplace=False))

        """
        **Tell me more about** :point_down:
        """

        peptide_sum_choices = list(ds.peptide_metadata.values)
        peptide_sum = st.selectbox(
            "Feature summary for:",
            peptide_sum_choices,
            0
        )

        dtype = st.session_state.peptide_table.dtypes[peptide_sum]
        p_uniq = len(set(st.session_state.peptide_table[peptide_sum]))
        p_kwargs={}
        if dtype in [pd.Int64Dtype(), pd.Float64Dtype()] and p_uniq < 50:
            view = st.radio(
                    "Peptide features display options", 
                    ["Quantiles", "Value Counts"])
            p_kwargs["numeric_dis"] = (view=="Quantiles")

        des = string_feature(
            ds, 
            feature= peptide_sum, 
            verbosity = 0, 
            dim="peptide",
            **p_kwargs
        )
        st.text(des)

        """
        **Apply a query condition to current peptide table** :point_down:
        """
        # Add Query
        st.text_input(
            label=f"Peptide Query Condition",
            key=f"pq-{num_queries+1}",
            on_change=add_query_condition,
            args=tuple([f"pq-{num_queries+1}", "peptide"])
        )

        """
        **Remove (using row index key) a query condition from the current peptide table** :point_down:
        """
        
        # Remove Query
        st.text_input(
            label=f"Remove Peptide Condition (by key)",
            key=f"rm_p_key_{num_dropped_queries}",
            on_change=drop_query_condition,
            args=tuple([f"rm_p_key_{num_dropped_queries}"])
        )

        """
        **Click below to see the raw peptide annotation table** :point_down:
        """
        st.write(st.session_state.peptide_table)

"""
***************************************
## :fire:    Visualize Enrichment Heatmap
Here, you can select options for rendering a heatmap based upon the samples-peptide enrichments available 
(as described in the tables above)
"""

settings, viz = st.columns([1,3])

#with st.expander("Heatmap settings"):
with settings:
    with st.form("dt"):

        """
        **Enrichment Layer** - select which enrichment transformation layer 
        you would like to visualize
        """

        enrichment_options = []
        for dt in set(list(ds.data_vars)) - set(["sample_table", "peptide_table"]):
            if ds[dt].values.flatten().min() != ds[dt].values.flatten().max():
                enrichment_options.append(dt)

        enrichment = st.selectbox(
            "Normalization layer",
            enrichment_options
        )    

    
        """
        **Sample Groups** - select which sample annotation groups would you like on the *y-axis*
        """
        adv_help = st.sidebar.button("?", key="q_help")
        #p_kwargs["numeric_dis"] = (view=="Quantiles")

        if adv_help:
            st.info("""
                Click on the Advanceds options setting
                to be able to apply either sample
                or peptide features on *both* axis.
                Use with caution this usually
                does not produce sensible vizualizations
            """)
        adv = st.radio(
                "Advanced axis options", 
                ["On", "Off"],
                index=1
        )
        
       # How to group the samples
        index=0
        y_choices = list(ds.sample_metadata.values)
        
        if adv == "On":
            y_choices += list(ds.peptide_metadata.values)

        #if "patient_status" in y_choices:
        #    index=y_choices.index("patient_status") + 1
        y = st.selectbox(
            "y-axis sample feature",
            ["sample_id"] + y_choices,
            index=index
        )
    
        """
        **Peptide Groups** - select which peptide annotation groups you would like on the *x-axis*
        """
    
        index=0
        x_choices = list(ds.peptide_metadata.values)
        #if "Protein" in x_choices:
        #    index=x_choices.index("Protein") + 1


        if adv == "On":
            x_choices += list(ds.sample_metadata.values)

        x = st.selectbox(
            "x-axis peptide feature",
            ["peptide_id"] + x_choices,
            index=index
        )

        #if x == y: st.warning("")

        """
        **Aggregation function** - when selecting axis which may group individual sample-peptide, how would you like to
        summarize the groups within a single entry in the resulting heatmap.
        """
        
        agg_func_choices = [
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
        **Split groups** Split the vizualization into groups (subplots).
        """

        #index=0
        r_s = get_reasonable_features(st.session_state.sample_table)
        r_p = get_reasonable_features(st.session_state.peptide_table)
        facet_choices = r_s + r_p
        #print(facet_choices)
        
        facet_features = st.multiselect(
            "Facet feature",
            ["None"] + facet_choices,
        )
        #print(facet_features)

        centered = st.radio(
                "centered at 0", 
                ["On", "Off"],
                index=1
        )
    
        domain_max = st.number_input("domain max")
        domain_min = st.number_input("domain min")
        #zmid = st.button("zero centered color")

        save_dir = st.text_input(label=f"save directory")

        #domain_max = st.number_input("domain max")
        heatmap_render = st.form_submit_button("Render Heatmap")

with viz:
    if heatmap_render:
        # TODO, we'll want to check the axis they've chosen
        # are unique or throw a warning??
        yss = [y] if y != 'sample_id' else []
        xss = [x] if x != 'peptide_id' else []
        sm, pm = [], []
        for f in facet_features + xss + yss:
            
            if infer_dim(f) == 'sample':
                sm.append(f)
            else:
                pm.append(f)

        #print(f"sample metadata: {sm}")
        #print(f"peptide metadata: {pm}")
        
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
       
        tds = to_tidy(subset_ds)

        scale_args = {}
        if domain_max and not domain_min:
            scale_args["domain"] = [0, domain_max]
        elif not domain_max and domain_min:
            scale_args["domain"] = [domain_min, 0]
        elif domain_max and domain_min:
            scale_args["domain"] = [domain_min, domain_max]

        if centered == "On":
            scale_args["domainMid"] = 0

        kwargs = {}
        kwargs["scale"] = alt.Scale(**scale_args) if len(scale_args) != 0 else {}
        color = alt.Color(f'{agg_func}({enrichment}):Q', **kwargs)

        if len(facet_features) > 0:
            for group, group_df in tds.groupby(facet_features):
                if len(facet_features) == 1:
                    title = f"{facet_features[0]}: {group}"
                else:
                    title = []
                    for j in range(len(facet_features)):
                        title.append(f"{facet_features[j]}: {group[j]}")

            
                c = alt.Chart(group_df).mark_rect().encode(
                    x=alt.X(
                        f'{x}:O', 
                        axis=alt.Axis(labelOverlap=True)
                    ),
                    y=alt.Y(
                        f'{y}:O', 
                        axis=alt.Axis(labelOverlap=True)
                    ),
                    color=color,
                    tooltip = [
                        alt.Tooltip(f'{agg_func}({enrichment}):O'),
                        alt.Tooltip(f'count({enrichment}):O'),
                    ]
                ).properties(
                    width=1000,
                    title=title
                )
                #base.save("test-plot.png")
                
                #if len(facet_features) > 0:
                    #c = c.facet(row=facet_col_name).resolve_scale(y='independent')
                    #c = c.facet(row=facet_col_name, title=None)

                st.altair_chart(c, use_container_width=True)
                if save_dir:
                    group_state = [group] if type(group) != tuple else group
                    group_state = [str(gs) for gs in group_state]
                    if not os.path.exists(save_dir): os.mkdir(save_dir)
                    save(c, os.path.join(save_dir, "-".join(group_state))+".png")
            #chart = alt.hconcat()
            #for group, data in  
            #    chart |= base.transform_filter(
            #            datum.species == species
            #    )
        else:

        
            c = alt.Chart(tds).mark_rect().encode(
                x=alt.X(
                    f'{x}:O', 
                    axis=alt.Axis(labelOverlap=True)
                ),
                y=alt.Y(
                    f'{y}:O', 
                    axis=alt.Axis(labelOverlap=True)
                ),
                color=color,
                tooltip = [
                    alt.Tooltip(f'{agg_func}({enrichment}):O'),
                    alt.Tooltip(f'count({enrichment}):O'),
                ]
            ).properties(
                width=1000,
            )

            st.altair_chart(c, use_container_width=True)
