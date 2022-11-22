# Generate a dataset of varying sequence length

rm(list=ls())

set.seed(23681732)

get_transition_matrix <- function(nLevels, p,b){
  shrink <- function(x,p){sign(x) * abs(x)^p}

  to_cumsum <- function(x){
    for(i in 1:nrow(x)) x[i,] <- x[i,]/sum(x[i,])
    for(i in 1:nrow(x)) x[i,] <- cumsum(x[i,])
    x[,ncol(x)] <- 1
    x
  }

  to_prob <- function(x){
    for(i in ncol(x):2) x[,i] <- x[,i]-x[,i-1]
    x
  }

  to_logit <- function(x){
    log(x/(1-x))
  }

  reverse_to_logit <- function(x){
    x <- exp(x)/(1+exp(x))
    x[,ncol(x)] <- 1
    x
  }

  X <- matrix(0,nrow=nLevels,nLevels)
  X <- 1/(abs(col(X)-row(X))+1)

  to_prob(reverse_to_logit(
    to_logit(to_cumsum(
      shrink(1/(abs(outer(1:nLevels,1:nLevels,"-"))+1),
             2)
    )) +
      -b
  ))
}

generate_one_patient <- function(p0,M,pEnd,maxSteps){

  nLevels <- length(p0)

  M <- rbind(
    cbind(M*t(sapply(1-pEnd,rep,6)),pEnd),
    c(rep(0,nLevels),1)
  )

  out <- data.frame(step=0:maxSteps,state=NA)

  currentState <- sample(1:nLevels,1,prob=p0)
  out$state[1] <- currentState

  i <- 1
  while(currentState <= nLevels &  i <= maxSteps){
    currentState <- sample(1:(nLevels+1),1,prob=M[currentState,])
    out$state[i+1] <- currentState
    i <- i+1
  }
  out[which(out$state %in% 1:nLevels),]
}

nLevels <- 6
p0List = list(even=rep(1,nLevels),
              good=1:nLevels,
              bad=nLevels:1
              )
p0List <- lapply(p0List, function(x){x/sum(x)})

pEndList = list(even=rep(1,nLevels),
              good=sqrt(1:nLevels),
              bad=sqrt(nLevels:1)
)
pEndList <- lapply(pEndList, function(x){0.2*x/sum(x)})



parGrid <- expand.grid(p0Index = 1:length(p0List),
                       pEndIndex = 1:length(pEndList),
                       b=c(-1,1)
                       )

n_per_group <- 30

shrinkPar <- 10
maxSteps <- 10

do.call("rbind",lapply(1:nrow(parGrid),function(i){
  p0 <- p0List[[parGrid[i,"p0Index"]]]
  pEnd <- pEndList[[parGrid[i,"pEndIndex"]]]
  M <- get_transition_matrix(nLevels,shrinkPar,parGrid[i,"b"])

  do.call("rbind",lapply(1:n_per_group,function(j){
    cbind(j,generate_one_patient(p0,M,pEnd,maxSteps))
  })) -> out

  out <- cbind(p0=parGrid[i,"p0Index"],
        pEnd=parGrid[i,"pEndIndex"],
        b=parGrid[i,"b"],
        out
        )

  out

})) -> example_data

example_data <- cbind(id=as.numeric(as.factor(paste(example_data$p0,example_data$pEnd,example_data$b,example_data$j))),
                      example_data)
example_data <- example_data[,c("id","step","state")]

example_data$id <- as.character(example_data$id)

usethis::use_data(example_data, overwrite = TRUE)







