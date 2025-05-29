<div align="center">
  <img width="1024px" src="https://raw.githubusercontent.com/ChenYangyao/dwarf_assembly_bias/master/docs/site_data/cover-github.jpg" alt="Dwarf Galaxies Adorning the Ink Wash Painting of Cosmic Web"/>
</div>


[![Last commit](https://img.shields.io/github/last-commit/ChenYangyao/dwarf_assembly_bias/master)](https://github.com/ChenYangyao/dwarf_assembly_bias/commits/master)
[![Workflow Status](https://img.shields.io/github/actions/workflow/status/ChenYangyao/dwarf_assembly_bias/run-test.yml)](https://github.com/ChenYangyao/dwarf_assembly_bias/actions/workflows/run-test.yml)
[![MIT License](https://img.shields.io/badge/License-MIT-blue)](https://github.com/ChenYangyao/dwarf_assembly_bias/blob/master/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/dwarf_assembly_bias)](https://pypi.org/project/dwarf_assembly_bias/)

The package holds the codes for the paper ***Unexpected clustering pattern in dwarf galaxies challenges formation models*** (Ziwen Zhang et al. [Nature 2025](https://www.nature.com/articles/s41586-025-08965-5), [arXiv:2504.03305](https://arxiv.org/abs/2504.03305)).

## What did we discover in this work?

A surprising pattern of spatial distribution was discovered in dwarf galaxies that diffuse ones cluster more strongly than compact ones - reversing the trends seen in massive galaxies. This challenges the structure formation models within the standard cold-dark-matter cosmology, calling for revisions for the understanding of the assembly of visible- and dark-matter structures.

## Installation

To install the package, run:
```bash
pip install dwarf-assembly-bias
```
Alternatively, you can clone the repository and install the package locally via `pip install -e /path/to/the/repo`.

## Resources for the paper

### Code samples

Code samples are organized as individual programs, scripts or jupyter-notebooks,
each for a specific procedure in producing the results in the paper.
All code samples are put under [docs/code_samples/](docs/code_samples):
- [cal_2pccf/](docs/code_samples/cal_2pccf/): program for two-point cross-correlation function (2PCCF) in observation.
- [cov_fit.py](docs/code_samples/cov_fit.py), [fit_bias.py](docs/code_samples/fit_bias.py): measuring the covariance matrix and fitting the relative bias.
- [Mhalo.py](docs/code_samples/Mhalo.py): HI-based halo mass estimation by assuming Burkert profile.
- [galaxy_web_cross.py](docs/code_samples/galaxy_web_cross.py): program for galaxy-cosmic web (knots, filaments, sheets and voids) 2PCCFs, for a sample of galaxies with given real-space locations. [galaxy_web_cross_from_obs_sample.py](docs/code_samples/galaxy_web_cross_from_obs_sample.py): is similar but takes input galaxies in redshift-space and make a correction for the redshift-space distortion before calculating the 2PCCFs.
- [theory.ipynb](docs/code_samples/theory.ipynb): theoretical interpretations (galaxy-galaxy 2PCCF in theory, abundance matching, model of self-interaction dark matter).

***Data for running the codes***: In order to run the code samples, please download the data
at [Quark cloud](https://pan.quark.cn/s/0fd2a303e03c), and put them under `docs/code_samples/data/`.

### Data for the figures

Here we provide all the data sheets, and scripts to load the data sheets and generate the figures exactly as those presented in the paper. All of these are held under [docs/figures/](docs/figures/). Specifically,
- [data/](docs/figures/data/): Excel data sheets, one for each Figure or Extended Data Figure.
- [plots_observation.py](docs/figures/plots_observation.py): scripts to generate observational figures (Fig. 1 and Extended Data Figs. 1-7).
- [plots_theory.ipynb](docs/figures/plots_theory.ipynb): scripts to generate theoretical figures (Figs. 2, 3, 4 and Extended Data Fig. 8).


## For developers

### Contribute to the project

Pull requests are welcome. For any changes, start a pull request to the ``dev`` branch and we will review it.

## For users

### Citation

We ask the users to cite the paper when using the package (code or data) in their research.

Users using specific module/code sample should also follow the copyright header in the source files.


## Details of this work

### Sample construction

The samples and subsamples of observed dwarf galaxies are described in the main text (Methods: "The sample of dwarf galaxies").

To measure the 2PCCF from observation, two additional sets of samples are constructed:
- **The reference sample**, used as field tracers, is a magnitude-limited sample constructed from the NYU-VAGC sample following the following selection criteria described in 
Sec. 2.1 of [Lixin Wang et al. 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.1452W): the $r$-band Petrosian apparent magnitude $r<17.6$; the $r$-band Petrosian absolute magnitude in $[-24, -16]$; the redshift $0.01 < z < 0.2$. 
- **The random samples**, used to account for the observational selection effects,
is obtained according to the method described in Sec. 3.1 of [Cheng Li et al. 2006](https://ui.adsabs.harvard.edu/abs/2006MNRAS.368...21L) as follows. We generated ten duplicates for each galaxy in the reference sample and randomly place them in the SDSS survey area. All other properties, including stellar mass and redshift of the duplicates, 
are the same as those of the parent galaxy. The random sample thus has the same survey geometry, the same distributions of galaxy properties and redshift, as the reference sample.


## Acknowledgements

We thank Fangzhou Jiang for his open source project of self-interaction dark matter (Fangzhou Jiang et al. 2023; [ads](https://ui.adsabs.harvard.edu/abs/2023MNRAS.521.4630J); [github](https://github.com/JiangFangzhou/SIDM)). A copy the source code can be found under [src/dwarf_assembly_bias/sidm/](src/dwarf_assembly_bias/sidm/).

We thank Hui-Jie Hu for his subroutines of HI-based halo mass estimator. A copy of the source code can be found at [docs/code_samples/Mhalo.py](docs/code_samples/Mhalo.py).