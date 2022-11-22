#' Produce a markovian summary of a cluster
#'
#' @usage
#' getQuality(markovSummary,by=NULL)
#'
#' @param x A list of vectors, giving each sequence of states. The target to be clustered.
#' @param by An optional vector providing how the data should be grouped.
#'
#' @details
#' This function produces a Markovian summary of sequences of vectors.
#' It estimates the proportion of observations that begin in each state,
#' and the proportion of transitions between states at each step. This summary
#' uses implicit absorbing state to represent the end of the sequence.
#' If \code{by} is specified, the sequences are grouped and a summary is produced
#' for each group.
#'
#' The \emph{error} of the markovian summary is also reported. This measure
#' is the sum of square differences in the proportion of observed sequences
#' that are in each state at time \emph{t} compared to the proportion that
#' would be expected if the sequences follow a Markov process defined according
#' to the summary. This sum is calculated from the second entry of the sequence
#' up until the maximum sequence length within the group (overall if \code{by}
#' is not specified). This sum is over all observed states and does not include
#' the implicit absorbing state used to represent the end of the sequence.
#'
#'
#' @return A list of class \code{markovSummary} containing a Markovian summary of each cluster
#'
#' @references
#'
#' Johns H, Hearne J, Bernhardt J, Churilov L. Clustering clinical and health
#' care processes using a novel measure of dissimilarity for variable-length
#' sequences of ordinal states.
#' \emph{Statistical Methods in Medical Research}.
#' 2020;29(10):3059-3075. doi:10.1177/0962280220917174
#'
#' @examples
#'
#' levelList <- by(example_data,example_data$id,function(df){
#' df$state[order(df$step)]
#' })
#' cluster <- cutree(hclust(ordinalEditDistance(levelList,0.8,5)),k = 5)
#'
#' markovSummary(levelList,cluster)
#'
markovSummary <- function(x,by=NULL){

  if(is.null(by)){
    by <- rep(1,length(x))
  }

  if(length(x)!=length(by)) stop("sequence list and by vector must be of the same length")

  maxLevel <- max(sapply(x, max))

  lapply(unique(by),function(i){
    clusterSequences <- x[by==i]

    maxLength <- max(sapply(clusterSequences,length))

    out <- ordinalEditDistance:::run_markov_profile(clusterSequences,maxLevel,maxLength)

    out$transition <- out$transition
    out$start <- out$start/sum(out$start)


    out$levelTally <- out$levelTally/length(clusterSequences)

    out$start <- out$start/sum(out$start)

    for(j in 1:nrow(out$transition)) out$transition[j,] <- out$transition[j,]/sum(out$transition[j,])
    sapply(0:(maxLength-1),function(j){out$start %*% expm::`%^%`(out$transition[1:maxLevel,1:maxLevel],j) }) -> markovDistribution

    distError <- sum((markovDistribution[,-1] - (out$levelTally)[,-1])^2)

    list(start=out$start,transition=out$transition,error=distError,n=length(clusterSequences))
  }) -> profiles

  class(profiles) <- c("markovSummary",class(profiles))
  profiles

}


