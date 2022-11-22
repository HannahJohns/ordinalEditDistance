#' Calculate ordinal edit distance
#'
#' @description Calculates Ordinal Edit Distance for sequences of ordinal data
#'
#' @usage
#' ordinalEditDistance(x,a,p)
#'
#' @param x A list containing vectors of varying length. The target to be clustered.
#' @param a A number between 0 and 1 corresponding the the relative cost of appending/deleting a level
#' @param p A positive real number corresponding to a Minkowski-like distance parameter
#'
#' @return A \code{dist} object containing the distances between each sequence
#'
#' @details
#'
#' Ordinal Edit Distance is a measure of dissimilarity for sequences of ordinal data.
#' It calculates the number of steps to convert a shorter sequence \emph{A} into a
#' longer sequence \emph{B} under the following allowed operations:
#'
#' \ennumerate{
#'     \item Append the sequence by copying the final observed sequence.
#'           This process is repeated until the sequences are of the same length
#'     \item Increment/De-increment each ordinal level within the sequence until
#'           the sequences are identical.
#' }
#'
#' Successive increments/de-increments may be performed at additional cost by
#' raising the required number of increments required to the power \code{p}.
#' The relative cost of these two steps is given by \code{a}.
#' When \code{a=1} the distance \emph{D(A,B)} is governed entirely based on
#' sequence length. When \code{a=0} the distance \emph{D(A,B)} is governed
#' entirely based on observed ordinal levels, and sequence length does not
#' matter. A value of \code{a} between these extremes forms a weighted average
#' of these two considerations.
#'
#' Ordinal Edit Distance does not neccesarily obey the triangle inequality,
#' that is for any three sequences A, B and C it is possible that
#' \emph{(A,B) + D(B,C) < D(A,C)}, for any choice of \code{a} or \code{p}.
#' Care should therefore be taken if the method is
#' to be used with methods that are sensitive to the violation of this assumption.
#'
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
#' ### Extracting ordinal sequences from a data frame and clustering
#'
#' df <- example_data
#' levelList <- by(example_data,example_data$id,function(df){
#'   df$state[order(df$step)]
#' })
#'
#' ## Producing a dendrogram
#'
#' hc <- hclust(ordinalEditDistance(levelList,a=0.8,p=5))
#'
#' plot(hc)
#'
#' ## Producing clusters, joining to a data frame, and plotting results
#'
#' library(ggplot2)
#' library(tibble)
#'
#' cluster <- cutree(hclust(ordinalEditDistance(levelList,0.8,5)),k = 5)
#' cluster <- rownames_to_column(as.data.frame(cluster),"id")
#'
#' df <- merge(df,cluster)
#'
#' ggplot(df,aes(x=step,fill=factor(state)))+
#'        geom_bar()+
#'        scale_fill_brewer(palette="RdYlGn")+
#'        facet_wrap(~cluster,scales="free_y")
#'
ordinalEditDistance <- function(x,a,p){

  if(min(sapply(x,min))!=1) stop("levels should begin at 1")

  out <- ordinalEditDistance:::run_ordinalEditDistance(x,a,1-a,p)
  as.dist(out)
}




