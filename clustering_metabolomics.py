import os
import umap
import numpy as np
import pandas as pd
import seaborn as sns

def import_and_transform_data(file):
    data = pd.read_csv(file)
    ## TODO: transform log2fc into a wide matrix

def collect_UMAP_embeddings(data):
    reducer = umap.UMAP()
    data = data.values
    reducer.fit_transform(drug_metab_data)
    embedding.shape

def plotting_embeddings():
    plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=[sns.color_palette()[x] for x in penguins.species_short.map({"Adelie":0, "Chinstrap":1, "Gentoo":2})])
    plt.gca().set_aspect('equal', 'datalim')
    plt.title('UMAP projection of the Penguin dataset', fontsize=24)

    #plot results colored by drug sensitivity

    # HDBSCAN clustering on reduced space

def main():
    pass

if __name__ = "__main__":
    path_fig = "C:\\Users\\mauro\\Documents\\phd_results"
    path_data = "C:\\Users\\mauro\\polybox\\data_metab"
