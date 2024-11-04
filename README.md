<div align="center">
<img src="https://github.com/C0nc/TAICHI/blob/main/fig/logo.png" width="200px">

**A Python package for the Scalale and accurate identification condition-relevant niches from spatial omics data.**

---

<p align="center">
  <a href="https://doi.org/10.1101/2024.05.30.596656" target="_blank">Preprint</a>
</p>

</div>

Taichi is able to automatically identify condition-relevant niches, and offers the downstream analysis based on obtained niches.
</p>
<p align="center">
  <img src="https://github.com/C0nc/TAICHI/blob/main/fig/pipeline.jpg" width="800px">
</p>

## Getting started




Please refer to the  
- STARmap Simulation dataset [Tutorial][link-tutorial_1] 
- MERFISH Simulation dataset [Tutorial][link-tutorial_4] (Simulation perturbated condition data  <a href="https://drive.google.com/file/d/18GGKFVeZfD1hsl17hEdRoHFLspfGd4Se/view?usp=drive_link" target="_blank">link</a> and original control data <a href="https://drive.google.com/file/d/1x5WxAU89JtnwioU4YUKAvOdTlnLj-kxq/view?usp=sharing" target="_blank">link</a>)
- Slice-seq v2 DKD mouse disease dataset [Tutorial][link-tutorial_2]. (Can be downloaded by `pysodb` package)
- CODEX proteomics CRC dataset [Tutorial][link-tutorial_3]. (Can be downloaded by `pysodb` package)
- ROC curve for spatial proteomics data [link-tutorial_5]
- Tensor Decomposition [link-tutorial_7]
- scalability Evaluation [link-tutorial_6]

## Installation

1. Create a conda environment
```bash
conda create -n taichi-env
conda activate taichi-env
```
2. Install the Taichi dependency
```bash
mamba install squidpy scanpy -c conda-forge (squidpy == 1.3.0 for reproducing CCI in manuscript)
pip insall pygsp ipykernel
```
3. Install the MENDER for batch-free niches representation:
```bash
cd MENDER
python setup.py install
```

Install the `pysodb` for efficient download processed Anndata in h5ad format (https://pysodb.readthedocs.io/en/latest/) if you want to run the DKD and CRC related analysis

We suggest using `mamba` to install the dependencies.
Installing the latest version of the dependencies may lead to dependency conflicts.


4. Additional Spatial co-embedding methods

If you want to try Taichi with other spatial co-embedding methods (CellChater, STAGATE).You should install them first and run the following code on two simulation benchmarking
```bash
python run_sta.py (STAGATE)
python run_cc.py (CellCharter)
```


## Contribution

If you found a bug or you want to propose a new feature, please use the [issue tracker][issue-tracker].

[issue-tracker]: https://github.com/C0nc/TAICHI/issues
[link-tutorial_1]: https://github.com/C0nc/TAICHI/blob/main/Tutorial.ipynb
[link-tutorial_2]: https://github.com/C0nc/TAICHI/blob/main/DKD_analysis.ipynb
[link-tutorial_3]: https://github.com/C0nc/TAICHI/blob/main/crc_analysis.ipynb
[link-tutorial_4]: https://github.com/C0nc/TAICHI/blob/main/merfish_analysis.ipynb
[link-tutorial_5]: https://github.com/C0nc/TAICHI/blob/main/roc.ipynb
[link-tutorial_6]: https://github.com/C0nc/TAICHI/blob/main/scalability_test.ipynb
[link-tutorial_7]: https://github.com/C0nc/TAICHI/blob/main/tensor_decomp.ipynb

