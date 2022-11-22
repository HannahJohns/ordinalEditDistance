#' Evaluates potential clustering settings for ordinalEditDistance
#'
#' @param x A list containing vectors of varying length. The target to be clustered.
#' @param a A vector of potential values for the append/delete cost
#' @param p A vector of potential values for the Minkowski-like parameter
#' @param k A vector of potential number of clusters
#' @param cl An optional cluster object returned by \link[parallel]{makeCluster}
#'
#' @details
#' For details on \code{a} and \code{p}, see \link{ordinalEditDist}.
#'
#' This function will perform a full-factorial evaluation across \code{a}, \code{p} and \code{k},
#' returning measures of cluster quality and the smallest observed cluster for that setting.
#' Clusters are developed using heirarchical clustering with complete linkage.
#'
#' For details on cluster quality, see \link{getQuality}.
#'
#' If the number of options for these settings is large, this function can be slow to run.
#' It may be sped up by using the optional \code{cl} argument to invoke parallelised evaluation.
#'
#' One method for selecting an optimal set of parameters based on these is to identify the
#' Pareto-optimal set of parameters based on \emph{distinctiveness} and \emph{deviation}.
#' This subset of options may then be investigated futher to determine the final set of parameters.
#'
#' Note that \emph{distinctiveness} and \emph{deviation} should not be directly compared across different numbers of
#' clusters \code{k}. The optimal values of \code{a} and \code{p} should therefore be determined for each value of
#' \code{k} that is evaluated.
#'
#'
#' @references
#'
#' Johns H, Hearne J, Bernhardt J, Churilov L. Clustering clinical and health
#' care processes using a novel measure of dissimilarity for variable-length
#' sequences of ordinal states.
#' \emph{Statistical Methods in Medical Research}.
#' 2020;29(10):3059-3075. doi:10.1177/0962280220917174
#' @examples
#'
#'
#' levelList <- by(example_data,example_data$id,function(df){
#'   df$state[order(df$step)]
#' })
#'
#' results <- evaluateClusters(levelList,
#'                             a=seq(0,1,length.out=11),
#'                             p=seq(1,5,length.out=11),
#'                             k = c(2,3,4,5))
#'
#' ## Repeating this calculation using parallelisation
#' \dontrun{
#' library(parallel)
#' cl <- makeCluster(round(0.6*parallel::detectCores()))
#' results <- evaluateClusters(levelList,
#'                             a=seq(0,1,length.out=11),
#'                             p=seq(1,5,length.out=11),
#'                             k = c(2,3,4,5),cl = cl)
#' stopCluster(cl)
#'
#'
#' ## Identifying pareto optimal solutions
#'
#'
#' library(rPref)
#' pareto <- do.call("rbind",by(results,results$k,function(df){
#'   psel(df,high(distinctiveness) * low(deviation))
#' }))
#'
#' ## Plotting results
#'
#' library(ggplot2)
#' ggplot(results,
#'        aes(x=1-deviation,
#'            y=distinctiveness,
#'            color=as.factor(k))
#' ) +
#' geom_point()+
#' geom_point(data=pareto,size=4)+
#' geom_step(data=pareto,direction = "vh")+
#' geom_label(data=pareto,size=4,hjust=0,
#'            aes(label=sprintf("a=%0.2f, p=%0.2f",a,p))
#'            )
#'}
#'
evaluateClusters <- function(x,a,p,k,cl=NULL){

  parGrid <- expand.grid(a=a,p=p,k=k)

  t(pbapply::pbsapply(1:nrow(parGrid),function(i,sequences,parGrid){
    a <- parGrid[i,"a"]
    p <- parGrid[i,"p"]
    k <- parGrid[i,"k"]

    candidate_clusters <- cutree(hclust(ordinalEditDistance::ordinalEditDistance(sequences,a,p)),k=k)
    quality <- ordinalEditDistance::getQuality(
                  ordinalEditDistance::markovSummary(sequences,candidate_clusters)
                )

    c(a=a,p=p,k=k,quality)
  },
  sequences=x,
  parGrid=parGrid,
  cl=cl)) -> results

  as.data.frame(results)

}



