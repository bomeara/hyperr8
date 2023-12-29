
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hyperr8

<!-- badges: start -->
<!-- badges: end -->

The goal of hyperr8 is to properly handle datasets that have a rate
versus a time.

## Installation

You can install the development version of hyperr8 from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bomeara/hyperr8")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(hyperr8)
car_rates <- generate_car_simulation()
head(car_rates)
#>   distance     time      rate       dataset
#> 1 203.7439 1.684582 120.94626 simulated car
#> 2 212.6438 3.566397  59.62425 simulated car
#> 3 208.2435 2.738483  76.04336 simulated car
#> 4 210.0515 2.591704  81.04761 simulated car
#> 5 208.1002 3.246711  64.09569 simulated car
#> 6 206.9553 3.331362  62.12335 simulated car
```

This generates a dataset with 1,000 cars with an average speed of 70 mph
and an average driving time of 3 hours, so an average distance traveled
of 210 miles. The driving time and distance traveled are highly
correlated (0.99). Just looking at the raw distance versus time, we see
a strong correlation:

But if we plot the estimated rate versus time, we see a hyperbola:

And this is even more clear on a log-log plot:
