print('script to analyze HD metabolomics data paired with drug sensitivity')
#script for high-dimentional analysis of metabolomics data
print('script to analyze HD metabolomics data paired with drug sensitivity')


#import packages
import os
import umap
import numpy as np
import pandas as pd
import seaborn as sns

path_out = "C:\\Users\\mauro\\Documents\\phd_results"
path_data = "C:\\Users\\mauro\\polybox\\data_metab"

#import data

pd.read_csv()

#fit umap

for drug in drugs:
    reducer = umap.UMAP()
    drug_metab_data = data.values
    reducer.fit_transform(drug_metab_data)
    embedding.shape

    #plot results

    plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=[sns.color_palette()[x] for x in penguins.species_short.map({"Adelie":0, "Chinstrap":1, "Gentoo":2})])
    plt.gca().set_aspect('equal', 'datalim')
    plt.title('UMAP projection of the Penguin dataset', fontsize=24)

    #plot results colored by drug sensitivity

    # HDBSCAN clustering on reduced space
