import os
import umap
import numpy as np
import pandas as pd
import seaborn as sns
from itertools import compress
import matplotlib.pyplot as plt

def import_and_transform_data(file):
    data = pd.read_csv(file,index_col=0)
    return(data)

def collect_UMAP_embeddings(data, out_path):
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(data.values)
    return(embedding)


def plot_embeddings(embedding, meta, drug):

    #plot embeddings colored by drug sensitivity
    fig = sns.relplot(x=embedding[:, 0],
                      y=embedding[:, 1],
                      hue=meta.percent_change_GR,
                      palette = sns.color_palette("rocket", as_cmap=True))
    fig.set(xlabel=None, ylabel = None)
    fig.set(title = 'UMAP projection of: '+drug + ' col=GR24')
    fig.tight_layout()
    fig.savefig(drug + "_full_col=GR24.png")

    #plot embeddings colored by cell
    fig = sns.relplot(x=embedding[:, 0],
                      y=embedding[:, 1],
                      hue=meta.cell_line)
    fig.set(xlabel=None, ylabel = None)
    fig.set(title = 'UMAP projection of: '+drug+' col=cell')
    fig.tight_layout()
    fig.savefig(drug + "_full_col=cell.png")

    # plot embeddings colored by concentration
    fig = sns.relplot(x=embedding[:, 0],
                      y=embedding[:, 1],
                      hue=meta.conc,
                      palette = sns.color_palette("viridis", as_cmap=True))
    fig.set(xlabel=None, ylabel = None)
    fig.set(title = 'UMAP projection of: '+drug + ' col=conc')
    fig.tight_layout()
    fig.savefig(drug + "_full_col=conc.png")

def parameter_selection_HDBSCAN(data,nclus):
    #param optimization HDBSCAN


def main():
    pass

if __name__ = "__main__":
    path_fig = "C:\\Users\\mauro\\Documents\\phd_results"
    path_data = "C:\\Users\\mauro\\Documents\\phd_results\\log2fc_full"

    files = os.listdir(path_data)

    files = list(compress(files, ["metadata" in x for x in files]))


    drugs = set([x.split("_")[1] for x in files])

    for drug in drugs:
        print(drug)
        #drug = 'Cladribine'
        data = import_and_transform_data(path_data+"\\data_"+drug+'_log2fc.csv')
        metadata = import_and_transform_data(path_data+"\\metadata_"+drug+'_log2fc.csv')

        embeddings = collect_UMAP_embeddings(data,"test")

        os.chdir(path_fig)

        plot_embeddings(embedding=embeddings, drug=drug, meta = metadata)





