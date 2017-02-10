/*
 * conversion.h
 *
 *  Created on: 23 juil. 2013
 *      Author: soueidat
 */



#ifndef CONVERSION_H_
#define CONVERSION_H_
#include <Rcpp.h>

template<class out,class inp>
inline out convertMatrix(const inp& matrixinput)
{
  int rows = matrixinput.rows();
  int cols = matrixinput.cols();
  out matrixOutput(rows,cols);
   for(int i=0;i<rows;i++)
   {
     for(int j=0;j<cols;j++)
     {
       matrixOutput(i,j) = matrixinput(i,j);
     }
   }
   return matrixOutput;
}


template<class out,class inp>
inline out convertvector(const inp& vectorinput)
{
    int len = vectorinput.size();
    out vectoroutput(len);
    for (int i = 0; i < len; ++i) {
      vectoroutput(i) = vectorinput(i);
    }
    return vectoroutput;
}


#endif /* CONVERSION_H_ */