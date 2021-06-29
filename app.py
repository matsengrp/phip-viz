#!/usr/bin/env python3

from plotnine import *
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
import os
import pickle
import xarray as xr

# TITLE OF THE APP
st.title('PhIPseq Viz')

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
def load_data(input_file_path):
    ###
    # load the pickled binary dataset
    ###
    ds = pickle.load(open(input_file_path, "rb"))
    return ds.sample_table.to_pandas().infer_objects()

# The Dataset is fairly simple. a consize explanation ca be found
# in the README here: https://github.com/matsengrp/phippery
# for now, we can just look at the sample table
# Unfortunately, my xarray approach is slightly flawed in that it a data variable
# can only be of a single datatype, and this means all my metadata is stored as "Objects"
# should really find a better solution at some point. the "infer_objects" is my hack

###
# from the xarray Dataset Object, grab the sample table
###
sample_table = load_data(selected_input_file)

###
# alignment stats to visualize
###

# I want to start by visualizing a few different alignment stats in the following columns
# these all get added during the pipeline run

# "raw_total_sequenced:"
# "reads mapped:"
# "error rate:"
# "average quality:"
align_stats = [c for c in sample_table.columns if c.endswith(":")]

# Let the user select which alignment stats to display
selected_stats = st.sidebar.multiselect(
    "Alignment Stats",
    align_stats,
    default=align_stats
)

# Get the unique set of `facet` values which can be displayed
unique_facets = [
    "control_status",
    "library_batch",
    "seq_dir",
    "experiment",
]

selected_facet = st.sidebar.selectbox(
    "Facets",
    unique_facets,
    index=0
)

# Get the unique set of values to display
unique_types = list(set(sample_table[selected_facet]))
selected_types = st.sidebar.multiselect(
    "Types",
    unique_types,
    default=unique_types
)

#plotting

# generally, we want to facet the subplots by alignment stat
# and then one of a few helpful ways to split sample groups

# an example of control_status
log_scale = True

# seaborn <- preferred
fig, ax = plt.subplots(
        len(selected_types),
        len(selected_stats),
        figsize=[10,10]
)

for subplot_row, subplot_type in enumerate(selected_types):
    for subplot_col, align_stat in enumerate(selected_stats):
        ax_s = ax[subplot_row, subplot_col]
        sns.histplot(
                data=sample_table.query(
                    f"{selected_facet} == '{subplot_type}'"
                ), 
                x=align_stat, 
                log_scale=log_scale,
                ax = ax_s
        )
        if subplot_col == 0: ax_s.set_ylabel(subplot_type)
        if subplot_row != len(selected_types)-1: ax_s.set_xlabel("")

fig.subplots_adjust(hspace=0.5, wspace=0.5)
st.write(fig)
