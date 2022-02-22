import os, umap, hdbscan, pickle
import numpy as np
import pandas as pd
import seaborn as sns

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
