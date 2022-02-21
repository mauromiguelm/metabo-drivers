import os, umap, hdbscan
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as stats
from itertools import compress

import matplotlib.pyplot as plt

def import_and_transform_data(file):
    data = pd.read_csv(file,index_col=0)
    return(data)

def collect_UMAP_embeddings(data,
                            n_neighbors,
                            n_components,
                            random_state=None):
    reducer = umap.UMAP(n_neighbors=n_neighbors,
                        n_components=n_components,
                        random_state=random_state)
    embeddings = reducer.fit_transform(data.values)

    return(embeddings)

def generate_clusters(data,
                      min_cluster_size):

    clusters = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size,
                               ).fit(data.values)
    return(clusters)

def score_clusters(clusters,prob_threshold = 0.05):
    #prob scoring
    cluster_labels = clusters.labels_
    label_count = len(np.unique(cluster_labels))
    total_num = len(clusters.labels_)
    cost = (np.count_nonzero(clusters.probabilities_ < prob_threshold) / total_num)

    return(label_count, cost)


def objective(params,metadata, data, penalty):
    """
    Objective function for hyperopt to minimize, which incorporates constraints
    on the number of clusters we want to identify
    """

    embeddings = collect_UMAP_embeddings(data,
                                         n_neighbors=params['n_neighbors'],
                                         n_components=params['n_components'],
                                         random_state=params['random_state'])

    clusters = generate_clusters(data=embeddings,
                                 min_cluster_size=params['min_cluster_size'])

    label_count, cost = score_clusters(clusters, prob_threshold=0.05)

    # penalty on the cost function if no concentration effect
    samples = [condition[1] for condition in name_group[1].groupby('condition')['value']]
    f_val, p_val = f_oneway(*samples)

    if(p_val<0.05):
        penalty = penalty
    else:
        penalty = 0

    loss = cost + penalty

    return {'loss': loss, 'label_count': label_count, 'status': STATUS_OK}

def bayesian_search(data, space, num_evals):
    """
        Perform bayseian search on hyperopt hyperparameter space to minimize objective function
    """

    trials = Trials()
    fmin_objective = partial(objective, data=data, penalty)
    best = fmin(fmin_objective,
                space=space,
                algo=tpe.suggest,
                max_evals=num_evals,
                trials=trials)

    best_params = space_eval(space, best)
    print('best:')
    print(best_params)
    print(f"label count: {trials.best_trial['result']['label_count']}")

    return best_params, trials

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

def main():
    pass

if __name__ = "__main__":
    path_fig = "C:\\Users\\mauro\\Documents\\phd_results"
    path_data = "C:\\Users\\mauro\\Documents\\phd_results\\log2fc_full"

    files = os.listdir(path_data)
    files = list(compress(files, ["metadata" in x for x in files]))

    space = {'n_neighbors':range(2,20),
             'n_components':range(2,5),
             'min_cluster_size':range(4,20),
             'penalty':0.15,
             'n_evals':10,
             'random_state':13
             }

    params= {'n_neighbors': 3,
             'n_components': 2,
             'min_cluster_size': 5,
             'penalty': 0.15,
             'n_evals': 10,
             'random_state': 13
             }


    drugs = set([x.split("_")[1] for x in files])

    for drug in drugs:
        print(drug)
        #drug = 'Methotrexate'
        data = import_and_transform_data(path_data+"\\data_"+drug+'_log2fc.csv')
        metadata = import_and_transform_data(path_data+"\\metadata_"+drug+'_log2fc.csv')

        bayesian_search()

        os.chdir(path_fig)

        plot_embeddings(embedding=embeddings, drug=drug, meta = metadata)





