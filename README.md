# phip-viz
Visualization of PhIPseq datasets

## Install Prerequisites

To install prerequisites, first set up a python3 virtual environment:

```#!/bin/bash
# Create the virtual environment
python3 -m venv venv

# Activate the virtual environment
source venv/bin/activate

# Install the prerequisites for running this app
# Replacing $REPO_FOLDER with the folder containing this repository
python3 -m pip install -r $REPO_FOLDER/requirements.txt

# Install phippery
python -m pip install git+https://github.com/matsengrp/phippery.git@52b8c5fcd0c4c727fe760b17a3820a60eada2bf3
```

## Launch the app

To render data from `.phip` files, run the following command
from a folder which contains >=1 `.phip` file to display:

```#!/bin/bash

# Replace $REPO_FOLDER with the folder containing this repository
streamlit run $REPO_FOLDER/app.py
```
