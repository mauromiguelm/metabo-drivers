import os
import umap
import numpy as np
import pandas as pd
import seaborn as sns
from itertools import compress

def import_and_transform_data(file):
    data = pd.read_csv(file,index_col=0)
    return(data)

def collect_UMAP_embeddings(data, out_path):
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(data.values)
    return(embedding)


def plot_embeddings(embedding, drug, path):
    #plot embeddings
    plt.scatter(
    embedding[:, 0],
    embedding[:, 1])
    #c=[sns.color_palette()[x] for x in data.conc)
    plt.gca().set_aspect('equal', 'datalim')
    plt.title('UMAP projection of: '+drug , fontsize=24)
    plt.savefig(path+"_"+drug+".png")

    #plot embeddings colored by drug sensitivity
    #plot embeddings colored by cell
    # plot embeddings colored by concentration

def parameter_selection_HDBSCAN(data,nclus):
    #param optimization HDBSCAN


def main():
    pass

if __name__ = "__main__":
    path_fig = "C:\\Users\\mauro\\Documents\\phd_results"
    path_data = "C:\\Users\\mauro\\polybox\\data_metab\\log2fc"

    files = os.listdir(path_data)

    files = list(compress(files, ["metadata" in x for x in files]))


    drugs = set([x.split("_")[1] for x in files])

    for drug in drugs:
        print(drug)
        drug = 'Methotrexate'
        data = import_and_transform_data(path_data+"\\data_"+drug+'_log2fc.csv')

        embeddings = collect_UMAP_embeddings(data,"test")

        plot_embeddings(embeddings, drug, path = path_fig)





