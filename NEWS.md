
# cfdecomp v0.3.0

This is the third release of the cfdecomp package. Compared to version 2, this new version has a bug fix which only showed up when cluster resampling was used and due to sample size, no cluster was resampled more than once.

# cfdecomp v0.2.0

This is a second release of the cfdecomp package. Compared to version 1, this new version now has three added functions and added functionality to existing functions.

The added functions are:
* cluster.resample: a function that allows for easy resampling by cluster ID in long-format data.
* cfd.FUN: a decomposition function that can decompose the output of any function (as long as that function returns a vector)

The added functionality is:
* cluster.sample and cluster.name have been added to cfd.mean and cfd.quantile, allowing for cluster resampling in these functions.

See DESCRIPTION for more information.
