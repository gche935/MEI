#' @keywords internal
"_PACKAGE"

#' MEI: Testing measurement Invariance Across Groups, Times, and Levels
#'
#' This package provides a collection of functions for testing configural invariance, metric invariance, scalar invariance, comparing latent means and defined parameters
#' It can be applied to multiple- groups, multiple-time in longitudinal studies, multiple-source in congruence research, and multi-level confirmatory factor analysis
#' #'
#' @section Key Features:
#' \itemize{
#'   \item Using Bootstrapping or Monte Carlo simulations to estimate confidence intervals of differences in estimated parameters. Bootstrapping 2,000 samples, or one million sets of factor loadings and intercepts are simulated based on the parameter estimates and the variance-covariance matrix of the freely estimated factor loadings and intercepts in the configural invariance model (or partial metric invariance model), and the confidence intervals of differences in factor loadings and intercepts between groups are calculated.
#'   \item When more than two groups are tested for ME/I, for each combination of reference item and argument, first compares each pair of groups using a Type I error rate of 0.01 or 0.001, then applies the list-and-delete method to identify the set of invariant groups.
#'   \item After identifying the sets of groups with invariant factor loadings or intercepts in each pair of reference and argument items, the model with the least number of freely estimated factor loadings and intercepts (the most parsimonious model with the largest number of invariant items) will be selected. When multiple models have the same number of estimated parameters, the model with the best fit should be selected.
#'   \item Missing values are handled with full-information maximum likelihood.
#'   \item All Non-converged and non-admissible bootstrapped samples are removed from the calculations.
#'   \item All The identification of invariant sets is correct up to a maximum of 25 groups. If there are more than 12 sets of invariant groups, the remaining groups will be considered as non-invariant.
#' }
#'
#' @section How to Use:
#' To get started, load the package and explore its main functions.
#' For example, use `? Full_MEI` to view individual function documentation.
#'
