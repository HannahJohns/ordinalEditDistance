#' Gets measure of cluster quality based on Markovian summaries of ordinal sequences
#'
#' @usage
#' getQuality(x)
#'
#' @param x an object of class \code{markovSummary} returned by the \code{markovSummary()} function
#'
#' @details
#'
#' This function returns several measures of cluster quality developed for
#' Ordinal Edit Distance. As the distance measure does not obey the triangle
#' inequality, measures of average between- and within-cluster distance that would typically be
#' used to measure cluster quality may be misleading. Instead, measures of
#' \emph{distinctiveness} and \emph{deviation} were proposed as alternatives.
#'
#' To achieve this, the Markovian summary for each cluster is treated as an
#' "average" value for the cluster with a known, fixed number of parameters
#' determined by the number of levels in the ordinal sequences.
#' \emph{distinctiveness} is defined as the average distance between any two
#' Markov summaries, where distance is given by the sum of squared differences
#' between each parameter value. \emph{deviation} is the average error of the
#' Markovian summary (see \link{markovSummary}) across clusters.
#'
#' @return a named vector containing the distinctiveness, deviation and minimum group size within the clusters
#'
#' @references
#'
#' Johns H, Hearne J, Bernhardt J, Churilov L. Clustering clinical and health
#' care processes using a novel measure of dissimilarity for variable-length
#' sequences of ordinal states.
#' \emph{Statistical Methods in Medical Research}.
#' 2020;29(10):3059-3075. doi:10.1177/0962280220917174
#'
#'
getQuality <- function(x){

  if(!("markovSummary" %in% class(x))) stop("getQuality should only be used by objects returned by markovSummary()")

  k <- length(x)

  maxLevel <- length(x[[1]]$start)

  deviation <- mean(sapply(x,function(y){y$error}))

  minN <- min(sapply(x,function(y){y$n}))

  if(k>1){
    distinctiveness <- do.call("rbind",lapply(1:k,function(i){
      sapply(1:k,function(j){

        sum(
          (
            c(x[[i]]$start,x[[i]]$transition[1:maxLevel,1:maxLevel]) -
              c(x[[j]]$start,x[[j]]$transition[1:maxLevel,1:maxLevel])
          )^2
        )
      })
    }))

    distinctiveness <- 2*sum(distinctiveness[upper.tri(distinctiveness)])/(k*(k-1))
  } else {
    distinctiveness <- NA
  }

  c(distinctiveness=distinctiveness,deviation=deviation,minClusterSize=minN)
}
