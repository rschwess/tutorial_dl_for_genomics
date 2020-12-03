# DL for Genomics Tutorial
A short tutorial to using DL and keras for genomics with sequence data.

## Requirements
* python 3.5 or 3.6
  * modules:
    * keras
    * numpy
  * I recommend anaconda https://docs.anaconda.com/ or miniconda to manage environments
    * to create a self contained environment after conda install:
      * conda create -n tutorial python=3.6
      * conda activate tutorial (or: source activate tutorial)
      * conda install pip
      * conda install nb_conda (to see virtual environments in jupyter notebooks)
      * pip install numpy tensorflow keras
      * to leave type: source deactivate
      
* R v3+
  * packages:
    * ggplot2 
   
## Description
This is a short tutorial running through some principles of using DL for genomics, specifically for sequence classification. We will use simulated sequences from 4 different classes that are populated with transcption factor motifs and other DNA patterns. We will build and train a small convolutional neuronal network using keras to learn to classify from the sequence only. We will also apply this network for new predictions and inspect what the network has learned to be important features.

* [Colab Notebook](https://colab.research.google.com/drive/1EeMDxj0SPRTVlZMloberSUs4BxP8e7wb) a great way to run the tutorial on google hardware :). The colab notebook actually has some updated visualizations and niceties so I strongly suggest running  it.
* [Notebook](./dl_for_genomics_tutorial_notebook.md) running through the tutorial
* python [script](./dl_intro.py) to adjust and run the entire analysis
* input [data](./data)
* [helper](./helper) scripts for simulating data and visualizing kernel weights
