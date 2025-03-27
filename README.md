# Unexpected clustering pattern in dwarf galaxies challenges formation models

[![Last commit](https://img.shields.io/github/last-commit/ChenYangyao/dwarf_assembly_bias/master)](https://github.com/ChenYangyao/dwarf_assembly_bias/commits/master)
[![Workflow Status](https://img.shields.io/github/actions/workflow/status/ChenYangyao/dwarf_assembly_bias/run-test.yml)](https://github.com/ChenYangyao/dwarf_assembly_bias/actions/workflows/run-test.yml)
[![MIT License](https://img.shields.io/badge/License-MIT-blue)](https://github.com/ChenYangyao/dwarf_assembly_bias/blob/master/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/dwarf_assembly_bias)](https://pypi.org/project/dwarf_assembly_bias/)

The package holds the codes for the paper "Unexpected clustering pattern in dwarf galaxies challenges formation models" (Journal Name??? - Volume??? - Page???).

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
All code samples are put under ``docs/code_samples``:
- ``cal_2pccf/``: program for 2PCCF computation.
- ``cov_fit.py``, ``fit_bias.py``: measuring the covariance matrix and fitting the relative bias.
- ``cal_halomass.py``: HI-based halo mass estimation by assuming Burkert profile.
- ``theory.ipynb``: theoretical interpretations (galaxy-galaxy and galaxy-cosmic web 2PCCF, abundance matching, model of self-interaction dark matter).

### Data for the figures

Here we provide all the data sheets, and scripts to load the data sheets and generate the figures exactly as those presented in the paper. All of these are held under ``docs/figures``. Specifically,
- ``data/``: Excel data sheets, one for each Figure or Extended Data Figure.
- ``plots_observation.py``: scripts to generate observational figures (Fig. 1 and Extended Data Figs. 1-7).
- ``plots_theory.ipynb``: scripts to generate theoretical figures (Figs. 2, 3, 4 and Extended Data Fig. 8).


## For developers

### Contribute to the project

Pull requests are welcome. For any changes, start a pull request to the ``dev`` branch and we will review it.

## For users

### Citation

We ask the users to cite the paper when using the package (code or data) in their research.