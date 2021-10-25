
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semcloud

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/400257454.svg)](https://zenodo.org/badge/latestdoi/400257454)
<!-- badges: end -->

The goal of semcloud is to process the output of the Python Workflow of
the Nephological Semantics project (coming soon) and prepare it to be
used with [NephoVis](https://qlvl.github.io/NephoVis).

## Installation

You can install the github development version with:

``` r
devtools::install_github("montesmariana/semcloud")
```

## Example

You can find a notebook with the workflow this package is meant to
support as a vignette of the package. To access it, you’ll have to
specify that vignettes must be built:

``` r
devtools::install_github("montesmariana/semcloud", build_vignettes = T)
vignette('processClouds', 'semcloud')
```

## License

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.
