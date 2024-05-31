<div align="center">
<img src="https://github.com/C0nc/TAICHI/blob/main/fig/logo.png" width="200px">

**A Python package for the Scalale and accurate identification condition-relevant niches from spatial -omics data.**

---

<p align="center">
  <a href="https://www.biorxiv.org/content/10.1101/2023.01.10.523386v2" target="_blank">Preprint</a>
</p>

</div>

Taichi is able to automatically identify condition-relevant niches, and offers the downstream anaylsis based on obtained niches.
</p>
<p align="center">
  <img src="https://github.com/C0nc/TAICHI/blob/main/fig/pipeline.jpg" width="800px">
</p>

## Getting started

Please refer to the Simulation dataset [Tutorial][link-tutorial_1] and DKD mouse disease dataset [Tutorial][link-tutorial_2].

## Installation

1. Create a conda environment
```bash
conda create -n taichi-env
conda activate taichi-env
```
2. Install the Taichi dependency
```bash
mamba install squidpy scanpy -c conda-forge
pip insall pygsp ipykernel
```
3. Install the MENDER for niches embebdding:
```bash
cd MENDER
python setup.py install
```

We suggest using `mamba` to install the dependencies.
Installing the latest version of the dependencies may lead to dependency conflicts.

## Contribution

If you found a bug or you want to propose a new feature, please use the [issue tracker][issue-tracker].

[issue-tracker]: https://github.com/CSOgroup/cellcharter/issues
[link-docs]: https://cellcharter.readthedocs.io
[link-api]: https://cellcharter.readthedocs.io/en/latest/api.html
[link-tutorial_1]: https://github.com/C0nc/TAICHI/blob/main/Tutorial.ipynb
[link-tutorial_2]: https://github.com/C0nc/TAICHI/blob/main/Tutorial_real_data.ipynb
