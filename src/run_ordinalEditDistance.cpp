#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix run_ordinalEditDistance(List X, double appendCost, double incrementCost, double incrementPower)
{

  int nObs = X.length();
  NumericMatrix out(nObs,nObs);

  int nVars;
  for(int i=0;i<nObs;i++)
  {
    NumericVector A = X(i);
    int Alength = A.length();

    out(i,i)=0;
    for(int j=i+1; j<nObs;j++)
    {
      NumericVector B = X(j);
      int Blength = B.length();

      int lengthDiff = std::max(Alength,Blength)-std::min(Alength,Blength);

      double d = appendCost * lengthDiff;


      //max(A,B) - min(A,B) could be removed but std::abs is integer only and I can't get
      //the double version working in Rcpp for some reason

      for (int k=0; k<std::min(Alength,Blength);k++)
      {
        d = d + incrementCost * pow(std::max(A(k),B(k))-std::min(A(k),B(k)),incrementPower);
      }

      if(Alength > Blength)
      {
        for (int k=0; k<lengthDiff;k++)
        {
          d = d + incrementCost * pow(std::max(A(Blength+k),B(Blength-1))-
                                      std::min(A(Blength+k),B(Blength-1)),
                                      incrementPower);
        }
      }
      else if (Blength > Alength)
      {
        for (int k=0; k<lengthDiff;k++)
        {
          d = d + incrementCost * pow(std::max(A(Alength-1),B(Alength+k))-
                                      std::min(A(Alength-1),B(Alength+k)),
                                      incrementPower);
        }
      }

      out(i,j)=d;
      out(j,i)=d;
    }
  }


  CharacterVector names = X.names();

  rownames(out) = names;
  colnames(out) = names;


 return out;
}

