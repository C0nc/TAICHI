
import anndata as ad
import squidpy as sq
import cellcharter as cc
import pandas as pd
import scanpy as sc
import scvi
import numpy as np
import matplotlib.pyplot as plt
from lightning.pytorch import seed_everything


import time
from Taichi.model import Taichi
import time

def cellcharter(adata_list):
    adata = ad.concat(adata_list, label='slice_id')
    seed_everything(12345)
    scvi.settings.seed = 12345

    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)


    scvi.model.SCVI.setup_anndata(
        adata, 
        layer="counts", 
        batch_key='slice_id',
    )

    model = scvi.model.SCVI(adata)

    model.train(early_stopping=True, enable_progress_bar=True)
    
    adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)

    sq.gr.spatial_neighbors(adata, library_key='slice_id', coord_type='generic', delaunay=True, spatial_key='spatial', percentile=99)\

    cc.gr.aggregate_neighbors(adata, n_layers=3, use_rep='X_scVI', out_key='X_cellcharter', sample_key='slice_id')

    return adata


ctrl_adata = sc.read_h5ad('merfish_control.h5ad')

full_adata = sc.read_h5ad('/home/cuiyan/mms/mms.h5ad')

res_list = []

for i in full_adata.obs['slice_id'].unique():

    cond_adata = full_adata[full_adata.obs['slice_id'] == i].copy()


    ctrl_adata.obs['condition'] = 0

    cond_adata.obs['condition'] = 1

    run_adata = cellcharter([ctrl_adata, cond_adata])


    run_adata.obs['condition'] = run_adata.obs['condition'].astype('category')


    start_time = time.time()

    model = Taichi(run_adata, ct_obs='cell_type', slice_id='slice_id')

    model.label_refinement(use_rep='X_cellcharter')

    res = model.graph_diffusion()

    res = res[res.obs['condition'] == 1]

    end_time = time.time()

    print(f'Total Running Time {end_time - start_time}')

    res_list.append(res)

ad.concat(res_list, label='slice_id').write_h5ad('taichi_cc_merfish.h5ad')


adata = sc.read_h5ad('/home/cuiyan/mms/gt_starmap.h5ad')


res_list = []

for i in adata.obs['slice_id'].unique():

    cond_adata = adata[adata.obs['slice_id'] == i].copy()

    ctrl_adata = cond_adata[cond_adata.obs['Region'].isin([2, 3])].copy()

    ctrl_adata.obs['condition'] = 0

    cond_adata.obs['condition'] = 1


    start_time = time.time()

    run_adata = cellcharter([ctrl_adata, cond_adata])

    run_adata.obs['condition'] = run_adata.obs['condition'].astype('category')

    model = Taichi(run_adata, ct_obs='ct', slice_id='slice_id')


    model.label_refinement(use_rep='X_cellcharter')

    res = model.adata

    res = model.graph_diffusion()

    res = res[res.obs['condition'] == 1].copy()

    end_time = time.time()

    res_list.append(res)

    print(f'Total Running Time {end_time - start_time}')


ad.concat(res_list, label='slice_id').write_h5ad('taichi_cc_starmap.h5ad')
