import STAGATE_pyG as STAGATE 
import scanpy as sc 
import anndata as ad
import pandas as pd
import numpy as np 

import time
from Taichi.model import Taichi



def stagate(adata_list):
    for adata in adata_list:
        STAGATE.utils.Cal_Spatial_Net(adata, rad_cutoff=150)
    adata = ad.concat(adata_list, label='slice_id')
    adata.uns['Spatial_Net'] = pd.concat([adata_list[i].uns['Spatial_Net'] for i in range(len(adata_list))])
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata = STAGATE.Train_STAGATE.train_STAGATE(adata)
    return adata

ctrl_adata = sc.read_h5ad('merfish_control.h5ad')

full_adata = sc.read_h5ad('merfish_condition.h5ad')

res_list = []

for i in full_adata.obs['slice_id'].unique():

    print(f'running {i}')
    cond_adata = full_adata[full_adata.obs['slice_id'] == i].copy()

    import anndata as ad 

    ctrl_adata.obs['condition'] = 0

    cond_adata.obs['condition'] = 1


    run_adata = stagate([ctrl_adata, cond_adata])

    run_adata.obs['condition'] = run_adata.obs['condition'].astype('category')


    start_time = time.time()

    model = Taichi(run_adata, ct_obs='cell_type', slice_id='slice_id')


    model.label_refinement(use_rep='STAGATE')


    res = model.graph_diffusion()

    res = res[res.obs['condition'] == 1]

    end_time = time.time()

    print(f'Total Running Time {end_time - start_time}')

    res_list.append(res)

ad.concat(res_list, label='slice_id').write_h5ad('taichi_stagate_merfish.h5ad')



adata = sc.read_h5ad('gt_starmap.h5ad')


res_list = []

for i in adata.obs['slice_id'].unique():

    cond_adata = adata[adata.obs['slice_id'] == i].copy()

    ctrl_adata = cond_adata[cond_adata.obs['Region'].isin([2, 3])].copy()

    ctrl_adata.obs['condition'] = 0

    cond_adata.obs['condition'] = 1


    start_time = time.time()


    run_adata = stagate([ctrl_adata, cond_adata])

    run_adata.obs['condition'] = run_adata.obs['condition'].astype('category')

    model = Taichi(run_adata, ct_obs='ct', slice_id='slice_id')


    model.label_refinement(use_rep='STAGATE')

    res = model.adata

    res = model.graph_diffusion()

    res = res[res.obs['condition'] == 1].copy()

    end_time = time.time()

    res_list.append(res)

    print(f'Total Running Time {end_time - start_time}')


ad.concat(res_list, label='slice_id').write_h5ad('taichi_stagate_starmap.h5ad')
