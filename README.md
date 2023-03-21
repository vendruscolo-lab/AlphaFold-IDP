# AlphaFold-IDP
This repository provides with custom code and analysis scripts to generate structural ensembles of intrinsically disordered proteins by using AlphaFold generated residue distance maps as constraints using prior-structural ensembles. AlphaFold structural ensemble and figures generated in the publication can be found in [Zenobo](myzenobo)

## Layout
- `BME_IDP`: This folder contains a branch of the original [Bayesian Maximum Entropy](https://github.com/KULL-Centre/BME) original code to reweight structural ensembles using experimental information, but specially adapted in the analysis **notebooks** to incorporate pairwise aminoacid AlphaFold predicted distances as constraints. In  `BME_IDP/notebooks/` one is able to reproduce the results in our paper such as: 
  - Generation of the AlphaFold posterior structural ensemble by using reweighting based on AlphaFold pairwise distances predictions as constraints
  - Radius of gyration in the AlphaFold structural ensemble
  - Back calculated NMR chemical shifts from the  AlphaFold structural ensemble and comparisson to the experimental NMR chemical shifts.
  - Secondary structure prediction based on the AlphaFold structural ensemble.   
 by using as input prior trajectories of αβ and α-synuclein either MD based, CG based or FoldingDiff based, found in the datasets mentioned below.
- `AlphaFold prediction of pairwise distances`: In [this](https://colab.research.google.com/github/zshengyu14/colabfold_distmat/blob/main/AlphaFold2.ipynb) google colab notebook, the user can generate AlphaFold predictions of the means and standard deviation of pairwise aminoacid distances, from arbitrary aminoacid sequences.
- `Backmap`: This folder contains examples code on how one can generate atomistic ensembles from αβ coarse-grained ensembles, originating from [Calvados2](https://github.com/KULL-Centre/CALVADOS) or [Foldingdiff](https://github.com/microsoft/foldingdiff) packages, by using [PULCHRA]( https://cssb.biology.gatech.edu/skolnick/files/PULCHRA).
- `BME_IDP/notebooks/PDDD_comparisnon.py`,`BME_IDP/notebooks/run_comparisson_PDDF.sh` scripts perform the comparisson to PairDistanceDistribution Function derived from SAXS data (using the Raw-v2.1.4 software). The ncesessary backmapped trajectories of the Calvados2 original dataset (Pesce et. al.) have been deposited [here](myzenobo)

## Dataset

- The prior structural ensembles originating from Molecular Dynamics for αβ can be found [here](https://zenodo.org/record/4247321)
- The prior structural ensembles originating from Molecular Dynamics for α-synuclein can be found [here](https://statics.teams.cdn.office.net/evergreen-assets/safelinks/1/atp-safelinks.html)
- The prior structural ensemble originating from Calvados2 for αβ and a-synuclein can be found [here] (myzenobo)
- The prior structural ensemble originating from Foldingdiff for αβ can be found [here](myzenobo)
- The AlphaFold posterior structural ensembles for all the abovementioned priors can be found [here] (myzenobo) 
- The prior structural ensembles originating from Calvados2 for DSS1, NH6cmdd, ANAC046, Sic1, ProTa, GHR-ICD can be found [here](https://zenodo.org/record/7415039#.ZBnari0Rq1E)

## Reproducability information 

The analysis was performed on a single machine containing 12 Intel(R) Xeon(R) W-2133 CPU @ 3.60GH CPUs, with total system memory of 100 GB 
Running the notebooks that performs the AF-ensemble generation takes approximately < 10'
The google-colab notebook takes < 30' to generate the the AF pairwise aminoacid distance prediction for each protein.
The backmap analysis can take up to 3 hrs depending on the proteins tested.


## Authors
[Faidon Brotzakis (@fbrotzakis)](https://github.com/fbrotzakis)

[Shengyu Zhang (@zshengyu14)](https://github.com/zshengyu14)

[Michele Vendruscolo (@vendruscolo-lab)](https://github.com/vendruscolo-lab)
## Article

https://www.biorxiv.org/content/10.1101/2023.01.19.524720v1.full
