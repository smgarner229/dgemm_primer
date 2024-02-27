#include <iostream>
#include <iomanip>

#include "Accelerate/Accelerate.h"

// Definitions to make editting easier
#define matsize 4

void call_dgemm();

void print_square_matrix(int n, double * mat);

int main(int argc, char * argv[])
{

  call_dgemm();

  return 0;

}

void call_dgemm()
{

  double * A = new double[matsize * matsize];
  for(int i = 0; i < matsize*matsize; i++)
    A[i]=1.0*(i+1);
  
  /*
  A = {
      1  5  9  13
      2  6  10 14
      3  7  11 15
      4  8  12 16
      }
  // NOTE: Column major as is Fortran convention
  */

  // Print the matrix
  std::cout << "A = " << std::endl;
  print_square_matrix(matsize,A);


  // For storing the result
  double * C = new double[matsize*matsize];
  // Fill in to avoid weird memory initialization things
  std::fill_n(C,matsize*matsize,0.0);

  // Flags for transposition
  char N = 'N';
  char T = 'T';

  // Dimensions of the matrices
  int m = matsize;
  int n = matsize;
  int k = matsize;

  double alpha = 1.0;
  double beta = 0.0;

  // Leading dimensions of A & C
  // tells how big of a stride to take along each memory dimension
  int LDA = matsize;
  int LDC = matsize;

  // Actual call to the matrix multiplication 
  dgemm_(&N,&N,&m,&n,&k,&alpha,A,&LDA,A,&LDA,&beta,C,&LDC);
  // Print the result
  std::cout << "A * A = " << std::endl;
  print_square_matrix(matsize,C);


  // Actual call to the matrix multiplication 
  dgemm_(&N,&T,&m,&n,&k,&alpha,A,&LDA,A,&LDA,&beta,C,&LDC);
  // Print the result
  std::cout << "A * A.T = " << std::endl;
  print_square_matrix(matsize,C);


  return;


}

void print_square_matrix(int n, double * mat)
{
  for(size_t i = 0; i < n; i++)
  {
    for(size_t j = 0; j < n; j++)
    {
      std::cout << std::setprecision(4) << std::setw(7) << mat[i+j*n];
    }
    std::cout << std::endl;
  }
}

