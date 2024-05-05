import warnings
warnings.filterwarnings("ignore")
import MENDER
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.metrics import *
import time
from sklearn.model_selection import StratifiedShuffleSplit
import anndata as ad


import networkx as nx
from pygsp import graphs, filters, plotting, utils
from sklearn.cluster import MeanShift, KMeans
import squidpy as sq
import MENDER
import scanpy as sc
import anndata as ad
import scipy.stats
import sklearn.linear_model
from sklearn.cluster import KMeans

def map_to_binary(values, case_control_labels, case_cond=1):
    
    # simple thresholding works, since the values are either very close to -1 or very close to 1
    # return np.interp(values, [np.min(values),np.max(values)], [0,1])>=0.5
    
    # but for consistency,
    # using kmeans with n_clusters=2, same as with HiDDEN p_hat binarization of the cells in the case_cond
    
    kmeans_values_res = KMeans(n_clusters=2, random_state=0).fit(pd.DataFrame(values))
    mean_values_res_kmeans_label0 = np.mean(values[kmeans_values_res.labels_==0]) 
    mean_values_res_kmeans_label1 = np.mean(values[kmeans_values_res.labels_==1])
    zero_lab_has_lower_mean = mean_values_res_kmeans_label0 < mean_values_res_kmeans_label1

    df_values_clust = pd.DataFrame(values)
    df_values_clust['kmeans'] = 0
    df_values_clust['kmeans'][(case_control_labels==case_cond).values] = [1 if x==int(zero_lab_has_lower_mean) else 0 for x in kmeans_values_res.labels_]
    
    return df_values_clust['kmeans'].values

def kmeans_clustering(y_true: np.ndarray, y_pred: np.ndarray, case_cond: int, rand_state: int) -> np.ndarray:
    """Binarize continous scores using a k-means strategy. This helps account for potential bi-modality."""
    is_case_cond = y_true == case_cond
    assert np.sum(is_case_cond), f'Found 0 examples of case condition {case_cond}!'
    kmeans_case = sklearn.cluster.KMeans(n_clusters=2, n_init='auto', random_state=rand_state)
    kmeans_case.fit(y_pred[is_case_cond].reshape(-1, 1))
    # We need to correct kmeans since the initialization can flip the labels.
    k_labels = kmeans_case.labels_
    if np.unique(k_labels).shape[0] == 1:
        print('No positive area')
    mean_p_hat_kmeans_label0 = np.mean(y_pred[is_case_cond][k_labels == 0])
    mean_p_hat_kmeans_label1 = np.mean(y_pred[is_case_cond][k_labels == 1])
    zero_lab_has_lower_mean = mean_p_hat_kmeans_label0 < mean_p_hat_kmeans_label1
    p_label = np.zeros_like(y_pred)
    p_label[is_case_cond] = np.array([1 if x == int(zero_lab_has_lower_mean) else 0 for x in k_labels])
    return p_label


def logistic_regression(x: np.ndarray, y: np.ndarray, rand_state: int) -> np.ndarray:
    """Basic logistic regression."""
    model = sklearn.linear_model.LogisticRegression(random_state=rand_state, penalty=None)
    model.fit(x, y)
    predicted_prob = model.predict_proba(x)
    return predicted_prob[:, 1]

def logistic_predictions(x: np.ndarray, y: np.ndarray, case_cond: int, rand_state: int):
    """Linear strategy for getting probabilities (continous scores) and binarized labels using kmeans strategy."""
    y_prob = logistic_regression(x, y, rand_state)
    y_labels = kmeans_clustering(y, y_prob, case_cond, rand_state)
    return y_prob, y_labels


def mms(adata, ct_obs, batch_obs, scale, radius, nn_mode):
    msm = MENDER.MENDER(
        adata,
        batch_obs = batch_obs,
        ct_obs=ct_obs,
        random_seed=42
    )

    msm.prepare()
    msm.set_MENDER_para(
        # default of n_scales is 6
        n_scales=scale,

        # for single cell data, nn_mode is set to 'radius'
        nn_mode=nn_mode,

        # default of n_scales is 15 um (see the manuscript for why).
        # MENDER also provide a function 'estimate_radius' for estimating the radius
        nn_para=radius,
        
    )
    # construct the context representation
    msm.run_representation_mp(
        100
    )

    adata.obsm['mender'] = msm.adata_MENDER.X
    train_adata = adata.copy()

    y_prob, y_labels = logistic_predictions(train_adata.obsm['mender'], train_adata.obs['condition'].values, 1, 42)
    
    train_adata.obs['new_labels'] = y_labels
    train_adata.obs['new_labels'] = train_adata.obs['new_labels'].astype('category')
    
    adata_list = []
    
    for s in train_adata.obs[batch_obs].unique():
        
        adata = train_adata[train_adata.obs[batch_obs] == s].copy()
        
        sq.gr.spatial_neighbors(adata)
        
        G = graphs.Graph(adata.obsp['spatial_connectivities'])
        
        G.estimate_lmax()
        
        if (adata.obs['new_labels'].unique().shape[0] > 1) & (adata.X.shape[0] > 1) & (adata.obs['new_labels'].values.astype(int).sum() < adata.obs['new_labels'].shape[0]) & (adata.obs['new_labels'].values.astype(int).sum() > 1):
        
            taus = [10, 25, 50]
            g = filters.Heat(G, taus)
            s = g.filter(adata.obs['new_labels'].values.astype(int), method='chebyshev')

            adata.obs['score_10'] = s[...,0]
            adata.obs['score_25'] = s[...,1]
            adata.obs['score_50'] = s[...,2]

            #sq.pl.spatial_scatter(adata, shape=None, color=['score_10', 'score_25', 'score_50'])
        
            y_seg = kmeans_clustering(adata.obs['new_labels'].values, adata.obs['score_25'].values, 1, 42)

            adata.obs['new_labels'] = y_seg 
            adata.obs['new_labels'] = adata.obs['new_labels'].astype('category')
    
        adata_list.append(adata)
    
    return ad.concat(adata_list, label='slice_id')

