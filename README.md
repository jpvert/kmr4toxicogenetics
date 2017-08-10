# kmr4toxicogenetics

Kernel multitask regression for toxicogenetics

## Description

The package delivers a vignette containing `R` code to reproduce the experiments in the paper "Kernel multitask regression for toxicogenetics" (bioRxiv-171298).

## Installation

This package fully relies on [`kmr`](https://github.com/jpvert/kmr) package (version >= 0.1) publicly available from github and several other packages for building the vignette locally.

```r
# install.packages("devtools")
devtools::install_github(c("jpvert/kmr", "jpvert/kmr4toxicogenetics"))
```

## Vignette

```r
vignette("kmr4toxicogenetics")
```

## References

> Elsa Bernard, Yunlong Jiao, Erwan Scornet, Veronique Stoven, Thomas Walter, and Jean-Philippe Vert. "Kernel multitask regression for toxicogenetics." Submitted. 2017. [bioRxiv-171298](https://doi.org/10.1101/171298)

## Misc

- Authors: [Jean-Philippe Vert](http://cbio.ensmp.fr/~jvert/), [Yunlong Jiao](http://cbio.ensmp.fr/~yjiao/), Elsa Bernard
- License: GPL-3
