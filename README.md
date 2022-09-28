# Metabolic drivers of drug sensitivity in cancer cell lines

## Introduction

To reveal the contribution of metabolism to drug sensitivity, we screened 31 drugs across 21 cancer cell lines at 5 logs of concentration, generating over 34.000 metabolic profiles coupled to 1 Mio images of over 3000 treatment groups as shown in Figure 1. From these results, we provide a map of drug-metabolite interactions that are general across many cell lines, several being independent of drug growth effect, a known confounder in metabolomics analysis. We discovered 703 drug-metabolite interactions, 488 of them being involved resistance and 215 in sensitivity. Of those, 599 are independent of drug growth effect – i.e. happens prior to growth effects in the concentration domain.

<p align="center">
<img src="https://github.com/mauromiguelm/metabo-drivers/graphical_abstract_methods.png" width="600">
 </p align="center">

<div align="center">
  Figure 1. Graphical abstract of methods for sample processing.
</div>

## Installation

Install package confluencer:
```R
 if(!require(confluencer)) {
  devtools::install_github("mauromiguelm/confluencer")
}
```

From the terminal: `start myproject.Rproj`

Metabolomics: coming soon

## About
 I am currently working towards documenting and publishing.
