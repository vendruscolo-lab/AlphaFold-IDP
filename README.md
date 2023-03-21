# AlphaFold-IDP
This repository provides with custom code and analysis scripts to generate structural ensembles of intrinsically disordered proteins by using AlphaFold generated residue distance maps as constraints using prior-structural ensembles. 

## Layout
- `BME_IDP`: This folder contains a branch of the original [Bayesian Maximum Entropy](https://github.com/KULL-Centre/BME) original code from to reweight structural ensembles using experimental information, but specially adapted in the analysis *notebooks* to incorporate pairwise aminoacid AlphaFold predicted distances as constraints.
- `AlphaFold prediction of pairwise distances`: In [this](https://colab.research.google.com/github/zshengyu14/colabfold_distmat/blob/main/AlphaFold2.ipynb) google colab notebook, the user can generate AlphaFold predictions of the means and standard deviation of pairwise aminoacid distances, from arbitrary aminoacid sequences.

## Dataset

This repository contains the full code and some small example data to reproduce our results on the kinetic ensemble of amyloid-β 42. See also the original implementation of the constrained VAMPNets.


## Reproducability information 

## Authors
[Faidon Brotzakis (@fbrotzakis)](https://github.com/fbrotzakis)

[Shengyu Zhang (@zshengyu14)](https://github.com/zshengyu14)

[Michele Vendruscolo (@vendruscolo-lab)](https://github.com/vendruscolo-lab)
## Article

https://www.biorxiv.org/content/10.1101/2023.01.19.524720v1.full
