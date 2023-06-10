// NOTES are provided at bottom of file.

#ifndef MATRIX_H_
#define MATRIX_H_

// #include <cstdlib> // needed for the exit function
#include <cstdlib> // later versions of g++ require this
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "eigen.h"
using namespace std;

class matrix
{

public: // attributes
    static int num_Pt; // static duration (life of program) b/c it will used to count number of calls from model obj to matrix for calculation of transition probabilities

protected://attributes
   int n_rows;                                // number of rows
   int n_cols;                                // number of columns
   double * data_;                            // pointer to matrix entries as "one long array"

public://constructors & destructors
   matrix(int m=4, int n=4, double s=1.0);    // m by n matrix filled with s; initializer list results in a 4x4 default matrix if dimensions are not specified in the constructor
   matrix(const matrix& x);                   // override default copy constructor to ensure a correct copy is made when needed
   virtual ~matrix(void);                     // destructor

public://methods
   int rows(void) const {return n_rows;}      // this is a "get rows" method
   int cols(void) const {return n_cols;}      // this is a "get columns" method
   matrix identity();                         // construct an identity matrix of same dimensions
   matrix& operator = (const matrix& x);      // assignment operator must be coded to work from a matrix (compiler default will NOT produce a correct copy)
   void reallocate(int m, int n);             // reallocate memory for new n x m size; all values initialized to -1
   void set_diagonal(); // for a square matrix, set the diagonal entries so that rows sum to 0
   matrix t(); // transpose a matrix

   //NOTE: [row*n_cols + col] = (row, col)
   const double& operator() (int i, int j) const {return data_[i*n_cols + j];}  // "const" subscript operator: allows access to matrix elements via familiar subscript notation M(i,j)
   double& operator() (int i, int j) {return data_[i*n_cols + j];}              // "mutable" subscript operator: accessing M(i.j) this way permits change to element at M(i,j)
   matrix exponentiate(const double &t, int k);                                 // method 1: repeated matrix squaring
   //matrix exponentiate(const matrix &x, const double &t);                     // method 2: does not work!!!!
   //friend istream& operator >> (istream& in, matrix& x);                      // read elements of a matrix from a text file
   void diag(double x); // replace the diagonal of a matrix with the supplied value

   friend matrix operator&(const matrix &A, const matrix &B); // returns the Kronecker product of two matrices

public://speed issues have cropped up: seem to have been resolved!
   const double& get_entry(int i, int j) const {return data_[i*n_cols + j];}    // no gains with this function as compared to (i,j) overload above; delete it when safe
   void set_transition_prob(int type, double pi[], double t);                   // moved here to allow computation of TPs in tree; the reduced deep copying of a matrix by 1 in each node o tree and yielded an increase in speed
   void set_transition_prob(double pi[], double t, int Q_dim); // JRM
   //void set_transition_prob(double pi[], double t, int dim, double sqrt_pi[], double lambda[], double QBV[], double W[]); // JRM
   //void set_transition_prob(double pi[], double t, int Q_dim, double lambda[], double QBV[], double W[]); // JRM
};

//Overloaded operators
ostream& operator << (ostream &out, const matrix &x); // allows values of a matrix to be printed to screen by using cout
matrix operator + (const matrix &x, const matrix &y); // allows matrix addition via the + operator; returns a matrix
matrix operator * (const matrix &x, const matrix &y); // allows matrix multiplication via the * operator; returns a matrix
matrix operator * (const matrix &x, const double &y); // allow mult of a matrix by a scalar; returns a matrix
matrix operator ^ (const matrix &x, const int k);     // allow the matrix to be squared k times; i.e. M^3 = M*M*M
matrix operator&(const matrix &A, const matrix &B); // returns the Kronecker product of two matrices


#endif /* MATRIX_H_ */

/* NOTES:

1. The declaration matrix x; gives a 4 by 4 matrix with entry 1; m = n_rows and n = n_cols;

2. The matrix is represented (internally as a linear array for a faster implementation.
   Access to the matrix is via a paired operator (i,j) for ease of use. See also note 3 below.

3. The reason for choosing double& operator()(int i, int j) rather than
   double*& operator[](int i) is spelled out in the C++ faq Lite:
   http://www.parashift.com/c++-faq-lite/operator-overloading.html#faq-13.11 and
   http://www.parashift.com/c++-faq-lite/operator-overloading.html#faq-18.12

4. When the subscript operator is applied to an object that is non-const
   the compiler will call the non-const version.  When the
   subscript operator is applied to a const object the compiler will call
   the const version.

5. All "float" type variables were changed to "doubles" on 13 March 2009

6. Base class matrices can have any dimension; subclassed matrices (DNA_matrix,
   CODON_matrix, & AA_matrix) have fixed dimensions and cannot exist otherwise.

*/


