# incur

A software package designed for analyzing IncuCyte-derived data, with a strong focus on generalisable and automated non-linear regression for interpretation of growth curves and dose-response relationships.

## Installation

You can install incur via:

```{r}
remote::install_github("moorezac/incur")
```

## Usage

The principle goal of incur is to provide a user-friendly pipeline to export, segment, and analyse IncuCyte-derived data. 

However, the automated curve fitting functions are able to be used on any dataset, and includes option for shared parameters across groups, outlier detection, and upper/lower limits on parameters.

It also includes a pipeline for cell line authenticiation/validation via interrogation of SNP microarray data.

## Learn

See the vignettes.
