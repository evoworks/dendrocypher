#ifndef INT_MATRIX_H_
#define INT_MATRIX_H_


#include <cstdlib> // later versions of g++ require this
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

class int_matrix
{
protected://attributes
   int n_rows;                                  // number of rows
   int n_cols;                                  // number of columns
   int * data_;                                 // pointer to matrix entries as "one long array"

public://constructors & destructors
   int_matrix(int m=4, int n=4, int s=1);        // m by n matrix filled with s; initializer list results in a 4x4 default matrix if dimensions are not specified in the constructor
   int_matrix(const int_matrix& x);              // override default copy constructor to ensure a correct copy is made when needed
   ~int_matrix(void);                            // destructor

public://methods
   int rows(void) const {return n_rows;}          // this is a "get rows" method
   int cols(void) const {return n_cols;}          // this is a "get columns" method
   int_matrix identity();                         // construct an identity matrix of same dimensions
   int_matrix& operator = (const int_matrix& x);  // assignment operator must be coded to work from a matrix (compiler default will NOT produce a correct copy)
   void reallocate(int m, int n);                 // reallocate memory for new n x m size; all values initialized to 0

   const int& operator() (int i, int j) const {return data_[i*n_cols + j];}  // "const" subscript operator: allows access to matrix elements via familiar subscript notation M(i,j)
   int& operator() (int i, int j) {return data_[i*n_cols + j];}              // "mutable" subscript operator: accessing M(i.j) this way permits change to element at M(i,j)
   //friend istream& operator >> (istream& in, matrix& x);                      // read elements of a matrix from a text file
};

//Overloaded operators
   ostream& operator << (ostream &out, const int_matrix &x);                   // allows values of a matrix to be printed to screen by using cout
   int_matrix operator + (const int_matrix &x, const int_matrix &y);           // allows matrix addition via the + operator; returns a matrix
   int_matrix operator * (const int_matrix &x, const int_matrix &y);           // allows matrix multiplication via the * operator; returns a matrix
   int_matrix operator * (const int_matrix &x, const double &y);               // allow mult of a matrix by a scalar; returns a matrix
   int_matrix operator ^ (const int_matrix &x, const int k);                   // allow the matrix to be squared k times; i.e. M^3 = M*M*M

#endif /* INT_MATRIX_H_ */
