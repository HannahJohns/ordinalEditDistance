#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List run_markov_profile(List X, int nLevels, int maxLength)
{

  int nObs = X.length();

  // Tally counts of markov transitions
  IntegerVector startVector(nLevels);
  IntegerMatrix transitionMatrix(nLevels+1,nLevels+1);
  IntegerMatrix levelTally(nLevels,maxLength);

  for(int i=0;i<nObs;i++)
  {
    NumericVector A = X(i);
    int Alength = A.length();

    if(Alength>0)
    {

      startVector(A(0)-1)++;
      levelTally(A(0)-1,0)++;

      for(int j=0; j<Alength-1;j++)
      {
        transitionMatrix(A(j)-1,A(j+1)-1)++;
        levelTally(A(j+1)-1,j+1)++;
      }

      transitionMatrix(A(Alength-1)-1,nLevels)++;

    }
  }

  transitionMatrix(nLevels,nLevels)=1;


  List out = List::create(Named("start")=startVector,
                          Named("transition")=transitionMatrix,
                          Named("levelTally")=levelTally
                          );

   return out;
}

