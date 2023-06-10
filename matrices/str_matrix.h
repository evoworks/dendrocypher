/*
 * Seems to be UNTESTED!!!!!
 */

#ifndef STR_MATRIX_H_
#define STR_MATRIX_H_

#include <string>
#include <iostream>
#include <fstream>
using namespace std;

class str_matrix
{
protected://attributes
   int n_rows;                                  // number of rows
   int n_cols;                                  // number of columns
   string **matrix;                             // double pointer(pp) to character state matrix that works with strings

public://constructors & destructors
   str_matrix();                                 // default constructor sets matrix to 4x4 (entries = "NULL")
   str_matrix(int m, int n);                     // m by n matrix; initialized with str = "NULL" at each element of the matrix
   str_matrix(int m, int n, string str);         // m by n matrix; initialized by user with value passed with str
   str_matrix(const str_matrix& x);              // override default copy constructor to ensure a correct copy is made when needed
   ~str_matrix(void);                            // destructor

public://methods
   int rows(void) const {return n_rows;}          // IN-LINED: this is a "get rows" method
   int cols(void) const {return n_cols;}          // IN-LINED: this is a "get columns" method (in-lined)
   str_matrix& operator = (const str_matrix& x);  // assignment operator must be coded to work from a matrix (compiler default will NOT produce a correct copy)
   void reallocate(int m, int n);                 // reallocate memory for new n x m size; Initialization is different here; reallocation yields = "*"
   void show_row(int m);                          // print just row indexed by m

   const string& operator() (int i, int j) const {return matrix[i][j];}   // IN-LINED: "const" subscript operator: allows access to matrix elements via familiar subscript notation M(i,j)
   string& operator() (int i, int j) {return matrix[i][j];}               // IN-LINED: "mutable" subscript operator: accessing M(i.j) this way permits change to element at M(i,j)
   //friend istream& operator >> (istream& in, str_matrix& x);            // read elements of a matrix from a text file
};

//Overloaded operators
   ostream& operator << (ostream &out, const str_matrix &x);              // allows values of a matrix to be printed to screen by using cout; places space between each element (could easily be changes to tab, etc)


#endif /* STR_MATRIX_H_ */


