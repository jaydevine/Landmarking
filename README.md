# **A Registration and Deep Learning Approach to Automated Landmark Detection for GM**

This repository contains a) Python/Bash scripts for atlas construction and pairwise volumetric image registration, b) Procrustes landmark coordinate data for training and testing (with an example .dta file), c) a Julia script to implement the DNN, and d) an R script to project newly acquired shape data into your training space for optimization.  

## **Keywords**

Anatomical landmark, deep learning, geometry, image registration, micro-computed tomography, morphometrics.

## **Prerequisites**

1. Linux or macOS;
2. [Medical Imaging NetCDF (MINC) Toolkit](https://github.com/BIC-MNI/minc-toolkit-v2) with local and remote (compute cluster) installations. Note that MINC is a modality neutral imaging data format and associated set of tools and libraries developed at the Montreal Neurological Institute (MNI) and freely available online. More information can be found at the [MINC Wikibooks page](http://en.wikibooks.org/wiki/MINC);
3. Volumetric imaging data. Ideally the file sizes are 300 MB or less, but this is not a strict requirement;
4. [Python](https://www.python.org/downloads/) and associated packages;
5. [R](https://cran.r-project.org/bin/) and associated packages;
6. [Julia](https://julialang.org/downloads/) and associated packages.

## **Notes on Data**

The micro-computed tomography imaging data are currently unavailable due to file size and ongoing projects. However, if you are interested in these data, please contact the [Hallgrímsson lab](https://www.ucalgary.ca/morpho/personnel).

## **Authors**

Jay Devine, Jose D. Aponte, David C. Katz, Wei Liu, Lucas D. Lo Vercio, Nils D. Forkert, Christopher J. Percival, and Benedikt Hallgrímsson.

## **License**

This project is licensed under the GNU General Public License. See the [LICENSE file](./LICENSE.md) for details.

## **Acknowledgments**

This study was motivated by the automated landmarking code and results of [Percival et al. (2019)](https://onlinelibrary.wiley.com/doi/10.1111/joa.12973). We thank our funding sources (NIH, NSERC, CIHR, CFI) for their support. We also thank [Dmitri Rozmanov](https://www.ucalgary.ca/tieleman/people/dmitri-rozmanov) and [Doug Phillips](https://people.ucalgary.ca/~phillips/) for assisting with computations performed on the Helix and ARC compute clusters at the University of Calgary.
