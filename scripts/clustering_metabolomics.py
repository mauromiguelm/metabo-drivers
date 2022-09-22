import os, umap, hdbscan, hyperopt, pickle, h5py, re
import numpy as np
import pandas
import pandas as pd
import seaborn as sns
import scipy.stats as stats
from itertools import compress
import matplotlib.pyplot as plt
import umap.plot

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
                               ).fit(data)
    return(clusters)

def score_clusters(clusters,prob_threshold = 0.05):
    #prob scoring
    cluster_labels = clusters.labels_
    label_count = len(np.unique(cluster_labels))
    total_num = len(clusters.labels_)
    cost = (np.count_nonzero(clusters.probabilities_ < prob_threshold) / total_num)

    return(label_count, cost)


def objective(params,metadata, data, sensitivity_penalty=True):
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

    label_count, cost = score_clusters(clusters = clusters, prob_threshold=params['prob_threshold'])

    if(sensitivity_penalty== True):

        # penalty on the cost function if absense of concentration effect
        samples = [value[1].dropna() for value in metadata.groupby(clusters.labels_)['percent_change_GR']]

        samples = [x  for x in samples if len(x) > 0]

        if len(samples)>1:
            f_val, p_val = stats.f_oneway(*samples)
        else:
            p_val=1

        if p_val>0.05:
            penalty = params['penalty']
        else:
            penalty = 0

        loss = cost + penalty

        return({'loss': loss, 'label_count': label_count, 'status':  hyperopt.STATUS_OK})
    else:
        loss = cost

        return ({'loss': loss, 'label_count': label_count, 'status': hyperopt.STATUS_OK})


def bayesian_search(metadata, data, space, n_evals, sensitivity_penalty):
    """
        Perform bayseian search on hyperopt hyperparameter space to minimize objective function
    """

    trials = hyperopt.Trials()
    fmin_objective = hyperopt.partial(objective, data=data, metadata = metadata, sensitivity_penalty=sensitivity_penalty)
    best = hyperopt.fmin(fmin_objective,
                space=space,
                algo= hyperopt.tpe.suggest,
                max_evals=n_evals,
                trials=trials)

    best_params = hyperopt.space_eval(space, best)
    print('best:')
    print(best_params)
    print(f"label count: {trials.best_trial['result']['label_count']}")

    return(best_params, trials)

def plot_embeddings(embedding, meta, drug):

    #plot embeddings colored by drug sensitivity
    fig = sns.relplot(x=embedding[:, 0],
                      y=embedding[:, 1],
                      hue=meta.percent_change_GR,
                      palette = sns.color_palette("rocket", as_cmap=True))
    fig.set(xlabel=None, ylabel = None)
    fig.set(title = 'UMAP projection of: '+drug + ' col=GR24')
    fig.tight_layout()
    plt.savefig(drug + "_full_col=GR24.png")
    plt.close()

    #plot embeddings colored by cell
    fig = sns.relplot(x=embedding[:, 0],
                      y=embedding[:, 1],
                      hue=meta.cell_line)
    fig.set(xlabel=None, ylabel = None)
    fig.set(title = 'UMAP projection of: '+drug+' col=cell')
    fig.tight_layout()
    plt.savefig(drug + "_full_col=cell.png")
    plt.close()

    # plot embeddings colored by drug
    fig = sns.relplot(x=embedding[:, 0],
                      y=embedding[:, 1],
                      hue=meta.drug)
    fig.set(xlabel=None, ylabel=None)
    fig.set(title='UMAP projection of: ' + drug + ' col=drug')
    fig.tight_layout()
    plt.savefig(drug + "_full_col=drug.png")
    plt.close()

    # plot embeddings colored by concentration
    fig = sns.relplot(x=embedding[:, 0],
                      y=embedding[:, 1],
                      hue=meta.conc,
                      palette = sns.color_palette("viridis", as_cmap=True))
    fig.set(xlabel=None, ylabel = None)
    fig.set(title = 'UMAP projection of: '+drug + ' col=conc')
    fig.tight_layout()
    plt.savefig(drug + "_full_col=conc.png")
    plt.close()

    # plot embeddings colored by clusters
    fig = sns.relplot(x=embedding[:, 0],
                      y=embedding[:, 1],
                      hue=meta.clusters,
                      palette="Set1")
    fig.set(xlabel=None, ylabel=None)
    fig.set(title='UMAP projection of: ' + drug + ' col=cluster')
    fig.tight_layout()
    plt.savefig(drug + "_full_col=clusters.png")
    plt.close()


def main():
    pass

if __name__ = "__main__":

    path_fig = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures_mean"
    path_metab = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\metabolomicsData_processed"
    os.chdir(path_metab)
    raw_data = h5py.File('metabolomics_raw.h5','r+')
    normalized_data = h5py.File('MM4_Mean mean_norm_DATA.h5', 'r+')

    # clustering of raw metabolomics data

    metadata = pd.DataFrame(raw_data['samples/name'])

    metadata = [str(x).split("_") for x in metadata.iloc[:,0]]

    metadata = pd.DataFrame(metadata)

    metadata = metadata.fillna("Solvent")

    metadata["ctrls"] = ["poolCtrl" if x == 'poolCtrl' else "cells" for x in metadata.iloc[:, 6]]

    data = pd.DataFrame(raw_data.get('data'))

    hspace = {'n_neighbors': hyperopt.hp.choice('n_neighbors', range(2, 500)),
              'n_components': hyperopt.hp.choice('n_components', range(2, 3)),
              'min_cluster_size': hyperopt.hp.choice('min_cluster_size', range(10, 1000)),
              "prob_threshold": 0.1,
              'penalty': 0.15,
              'n_evals': 100,
              'random_state': 13,
              'sensitivity_penalty' : False
    }



    a, b = bayesian_search(metadata= metadata, data=data, space=hspace, n_evals=hspace['n_evals'],sensitivity_penalty=hspace['sensitivity_penalty'])



    #a = {'min_cluster_size': 164, 'n_components': 2, 'n_evals': 100, 'n_neighbors': 87, 'penalty': 0.15, 'prob_threshold': 0.1, 'random_state': 13, 'sensitivity_penalty': False}

    os.chdir(path_fig + "\\clustering_figures")

    from matplotlib.axes._axes import _log as matplotlib_axes_logger

    matplotlib_axes_logger.setLevel('ERROR')

    mapper = umap.UMAP(n_neighbors=a['n_neighbors'],
                                         n_components=a['n_components']).fit(data.values)
    umap.plot.points(mapper, labels=metadata.iloc[:,5],show_legend=False)
    plt.savefig('raw_cells.png')

    umap.plot.points(mapper, labels=metadata['ctrls'],show_legend=False)
    plt.savefig('raw_ctrls.png')

    #cluster norm metabolomics  data

    metadata = pd.DataFrame(normalized_data['samples/name'])

    metadata = [str(x).split("_") for x in metadata.iloc[:,0]]

    metadata = pd.DataFrame(metadata)

    metadata = metadata.fillna("Solvent")

    metadata["ctrls"] = ["poolCtrl" if x == 'poolCtrl' else "cells" for x in metadata.iloc[:, 6]]

    data = pd.DataFrame(normalized_data.get('data'))

    hspace = {'n_neighbors': hyperopt.hp.choice('n_neighbors', range(2, 500)),
              'n_components': hyperopt.hp.choice('n_components', range(2, 3)),
              'min_cluster_size': hyperopt.hp.choice('min_cluster_size', range(10, 1000)),
              "prob_threshold": 0.1,
              'penalty': 0.15,
              'n_evals': 100,
              'random_state': 13,
              'sensitivity_penalty' : False
    }

    a, b = bayesian_search(metadata= metadata, data=data, space=hspace, n_evals=hspace['n_evals'],sensitivity_penalty=hspace['sensitivity_penalty'])

    #{'min_cluster_size': 138, 'n_components': 2, 'n_evals': 100, 'n_neighbors': 38, 'penalty': 0.15,
    # 'prob_threshold': 0.1, 'random_state': 13, 'sensitivity_penalty': False}

    os.chdir(path_fig + "\\clustering_figures")

    from matplotlib.axes._axes import _log as matplotlib_axes_logger

    matplotlib_axes_logger.setLevel('ERROR')

    mapper = umap.UMAP(n_neighbors=a['n_neighbors'],
                                         n_components=a['n_components']).fit(data.values)
    umap.plot.points(mapper, labels=metadata.iloc[:,5], show_legend=False)
    plt.savefig('norm_cells.png')

    umap.plot.points(mapper, labels=metadata['ctrls'],show_legend=False)
    plt.savefig('norm_ctrls.png')



    #cluster log2fc

    path_fig = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures_mean"
    path_data = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data_mean"

    files = os.listdir(path_data+"\\log2fc")
    files = list(compress(files, ["metadata" in x for x in files]))

    #drug-by-drug clustering

    hspace = {'n_neighbors':hyperopt.hp.choice('n_neighbors',range(2,20)),
             'n_components':hyperopt.hp.choice('n_components',range(2,7)),
             'min_cluster_size':hyperopt.hp.choice('min_cluster_size',range(4,100)),
             "prob_threshold":0.1,
             'penalty':0.15,
             'n_evals':250,
             'random_state':13
              }

    drugs = set([x.split("_")[1] for x in files])

    results = {}

    for drug in drugs:
        """
        define best params for UMAP+HDBSCAN
        """
        print(drug)
        #drug = 'Methotrexate'
        data = import_and_transform_data(path_data+"\\log2fc_full"+"\\data_"+drug+'_log2fc.csv')
        metadata = import_and_transform_data(path_data+"\\log2fc_full"+"\\metadata_"+drug+'_log2fc.csv')

        a,b = bayesian_search(data,space = hspace,n_evals=hspace['n_evals'])

        results[drug] = [a,b]


    os.chdir(path_data)

    summary_best_params = [results[x][0] for x in results.keys()]

    summary_best_params = pd.DataFrame(summary_best_params)

    summary_best_params['drug'] = results.keys()

    n_clus = [results[x][1].best_trial['result']['label_count'] for x in results.keys()]

    summary_best_params['n_clus']= n_clus

    os.chdir(path_data)

    summary_best_params.to_csv("summary_clustering.csv")

    with open('opt_clustering.pkl', 'wb') as f:
        pickle.dump(results, f)

    results = []
    with (open("opt_clustering.pkl", "rb")) as f:
        while True:
            try:
                results.append(pickle.load(f))
            except EOFError:
                break



    for drug in drugs:
        """
        get best UMAP+HDBSCAN results
        """
        #drug = 'Methotrexate'

        data = import_and_transform_data(path_data+"\\log2fc_full" + "\\data_" + drug + '_log2fc.csv')
        metadata = import_and_transform_data(path_data+"\\log2fc_full" + "\\metadata_" + drug + '_log2fc.csv')

        params = results[drug][0]

        embeddings = collect_UMAP_embeddings(data,
                                             n_neighbors=params['n_neighbors'],
                                             n_components=params['n_components'],
                                             random_state=params['random_state'])



        clusters = generate_clusters(data=embeddings,
                                     min_cluster_size=params['min_cluster_size'])

        metadata['clusters'] = clusters.labels_

        #save results

        os.chdir(path_data+'\\cluster_labels')

        metadata.to_csv('clusterLabels_'+drug+'.csv')

        embeddings = collect_UMAP_embeddings(data,
                                             n_neighbors=params['n_neighbors'],
                                             n_components=2,
                                             random_state=params['random_state'])

        os.chdir(path_fig+"\\clustering_figures")

        plot_embeddings(embedding=embeddings,
                        meta=metadata,
                        drug = drug)


    # full data clustering

    hspace = {'n_neighbors': hyperopt.hp.choice('n_neighbors', range(2, 500)),
              'n_components': hyperopt.hp.choice('n_components', range(2, 3)),
              'min_cluster_size': hyperopt.hp.choice('min_cluster_size', range(10, 1000)),
              "prob_threshold": 0.1,
              'penalty': 0.15,
              'n_evals': 100,
              'random_state': 13,
              'sensitivity_penalty': True
              }


    drugs = set([x.split("_")[1] for x in files])

    results = {}

    data = []
    metadata = []

    for drug in drugs:
        """
        define best params for UMAP+HDBSCAN
        """
        print(drug)
        # drug = 'Methotrexate'
        data.append(import_and_transform_data(path_data + "\\log2fc" + "\\data_" + drug + '_log2fc.csv'))
        metadata.append(import_and_transform_data(path_data + "\\log2fc" + "\\metadata_" + drug + '_log2fc.csv'))

    data = pd.concat(data)
    metadata = pd.concat(metadata)


    a, b = bayesian_search(data=data,metadata = metadata, space=hspace, n_evals=hspace['n_evals'],sensitivity_penalty=hspace['sensitivity_penalty'])

    #{'min_cluster_size': 54, 'n_components': 2, 'n_evals': 100, 'n_neighbors': 297, 'penalty': 0.15,
    # 'prob_threshold': 0.1, 'random_state': 13, 'sensitivity_penalty': True}

    os.chdir(path_fig+'\\clustering_figures')

    mapper = umap.UMAP(n_neighbors=a['n_neighbors'],
                       n_components=a['n_components']).fit(data.values)

    umap.plot.points(mapper, labels=metadata['cell_line'], show_legend=False)
    plt.savefig('fc_cells.png')


    umap.plot.points(mapper, labels=metadata['drug'], show_legend=False)
    plt.savefig('log2fc_drug.png')

    os.chdir(path_data)

    with open('opt_clustering_fullData.pkl', 'wb') as f:
        pickle.dump([a,b], f)

    summary_best_params = a.items()
    summary_best_params = list(summary_best_params)
    summary_best_params = pd.DataFrame(summary_best_params)

    n_clus = [results[x][1].best_trial['result']['label_count'] for x in results.keys()]

    summary_best_params['n_clus'] = b.best_trial['result']['label_count']

    os.chdir(path_data)

    summary_best_params.to_csv("summary_clustering_fullData.csv")

    """
    get best UMAP+HDBSCAN results
    """

    params = a

    embeddings = collect_UMAP_embeddings(data,
                                         n_neighbors=params['n_neighbors'],
                                         n_components=params['n_components'],
                                         random_state=params['random_state'])

    clusters = generate_clusters(data=embeddings,
                                 min_cluster_size=params['min_cluster_size'])

    metadata['clusters'] = clusters.labels_

    embeddings = collect_UMAP_embeddings(data,
                                         n_neighbors=params['n_neighbors'],
                                         n_components=2,
                                         random_state=params['random_state'])

    os.chdir(path_fig + "\\clustering_figures")

    plot_embeddings(embedding=embeddings,
                    meta=metadata,
                    drug='all-drugs')























