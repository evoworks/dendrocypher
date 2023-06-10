#include "str_matrix.h"

/******************* Implementation of matrix ********************/


//-----------------------------------------------------------------
// CONSTRUCTOR_1: default constructor sets matrix to 4x4 (entries = "NULL")
//-----------------------------------------------------------------
str_matrix::str_matrix()
{
	int i,j;
   n_rows = 4;
   n_cols = 4;

   matrix = new string *[n_rows];                              // "new" allocates memory for data_matrix rows
   matrix[0] = new string[n_rows*n_cols];                      // "new" allocates memory for data_matrix rows * cols
   for(i=1; i<n_rows; i++) matrix[i] = matrix[0] + i*n_cols;   // ptr math: sets n_row pts to elements allocated for the rows

   // initialize all elements to "NULL"
   for(i=0; i<n_rows; i++)
      for(j=0; j<n_cols; j++)
         matrix[i][j] = "NULL";
}


//-----------------------------------------------------------------
// CONSTRUCTOR_2:
//-----------------------------------------------------------------
str_matrix::str_matrix(int m, int n)
{
	int i,j;
   n_rows = m;
   n_cols = n;

   matrix = new string *[n_rows];                              // "new" allocates memory for data_matrix rows
   matrix[0] = new string[n_rows*n_cols];                      // "new" allocates memory for data_matrix rows * cols
   for(i=1; i<n_rows; i++) matrix[i] = matrix[0] + i*n_cols;   // ptr math: sets n_row pts to elements allocated for the rows

   // initialize all elements to "NULL"
   for(i=0; i<n_rows; i++)
      for(j=0; j<n_cols; j++)
         matrix[i][j] = "NULL";
}

//-----------------------------------------------------------------
// CONSTRUCTOR_3:
//-----------------------------------------------------------------
str_matrix::str_matrix(int m, int n, string str)
{
	int i,j;
   n_rows = m;
   n_cols = n;

   matrix = new string *[n_rows];                               // "new" allocates memory for data_matrix rows
   matrix[0] = new string[n_rows*n_cols];                       // "new" allocates memory for data_matrix rows * cols
   for(i=1; i<n_rows; i++) matrix[i] = matrix[0] + i*n_cols;    // ptr math: sets n_row pts to elements allocated for the rows

   // initialize all elements to value supplied in str
   for(i=0; i<n_rows; i++)
      for(j=0; j<n_cols; j++)
         matrix[i][j] = str;
}


//-----------------------------------------------------------------
// CONSTRUCTOR: override default copy constructor
//-----------------------------------------------------------------
str_matrix::str_matrix(const str_matrix& x)
{
   int i,j;
   n_rows = x.n_rows;
   n_cols = x.n_cols;

   matrix = new string *[n_rows];                               // "new" allocates memory for data_matrix rows
   matrix[0] = new string[n_rows*n_cols];                       // "new" allocates memory for data_matrix rows * cols
   for(i=1; i<n_rows; i++) matrix[i] = matrix[0] + i*n_cols;    // ptr math: sets n_row pts to elements allocated for the rows

   for(i=0; i<n_rows; i++)
      for(j=0; j<n_cols; j++)
      	matrix[i][j] = x.matrix[i][j];
}


//-----------------------------------------------------------------
// DESTRUCTOR: release matrix memory
//-----------------------------------------------------------------
str_matrix::~str_matrix(void)
{
   delete [] matrix[0];
   delete [] matrix;
}


//-----------------------------------------------------------------
// ASSIGNEMT OPERATOR: default assignment operator will not produce
// a correct copy.
//-----------------------------------------------------------------
str_matrix& str_matrix::operator = (const str_matrix& x)
{
   int i,j;
   if (this != &x)   // watch out for x = x
   {
      if( !((n_rows == x.n_rows)&&(n_cols == x.n_cols)) ) // must call reallocate
      {
      	cout << "\nWARNING::str_matrix::operator=:: assignment of different sized matrices!\n";
      	cout << "Reallocating the target matrix; you have been warned!\n\n";
         delete [] matrix[0];
         delete [] matrix;
         n_rows = x.n_rows;
         n_cols = x.n_cols;

         matrix = new string *[n_rows];                               // "new" allocates memory for data_matrix rows
         matrix[0] = new string[n_rows*n_cols];                       // "new" allocates memory for data_matrix rows * cols
         for(i=1; i<n_rows; i++) matrix[i] = matrix[0] + i*n_cols;    // ptr math: sets n_row pts to elements allocated for the rows
      }
      for(i=0; i<n_rows; i++)
         for(j=0; j<n_cols; j++)
         	matrix[i][j] = x.matrix[i][j];
   }
   return *this;
}


//-----------------------------------------------------------------
// PUBLIC: delete old memory allocation and allocate new matrix
// Initialization is different here; reallocation yields = "*"
//-----------------------------------------------------------------
void str_matrix::reallocate(int m, int n)
{
   delete [] matrix[0];
   delete [] matrix;

   int i,j;
   n_rows = m;                           // set new number of rows
   n_cols = n;                           // set new number of columns

   matrix = new string *[n_rows];                              // "new" allocates memory for data_matrix rows
   matrix[0] = new string[n_rows*n_cols];                      // "new" allocates memory for data_matrix rows * cols
   for(i=1; i<n_rows; i++) matrix[i] = matrix[0] + i*n_cols;   // ptr math: sets n_row pts to elements allocated for the rows

   // a reallocation in initialized so all elements = "*"
   for(i=0; i<n_rows; i++)
      for(j=0; j<n_cols; j++)
         matrix[i][j] = "*";
}


//-----------------------------------------------------------------
// PUBLIC: print just row indexed by m
//-----------------------------------------------------------------
void str_matrix::show_row(int m)
{
   for(int j=0; j<n_cols; j++) cout << matrix[m][j] << "\t";
}


/*
//-----------------------------------------------------------------
// OVERLOADED OPERATOR: to read matrix from a file
//-----------------------------------------------------------------
istream& operator >> (istream &in, str_matrix &x)
{
   int m, n, i, j;

   in >> m >> n;

   // check whether the sizes are compatible
   if(!( (m == x.n_rows) && (n == x.n_cols) ) )
   {
   	delete [] matrix[0];
   	delete [] matrix;

   	x.matrix = new string *[n_rows];                              // "new" allocates memory for data_matrix rows
   	x.matrix[0] = new string[n_rows*n_cols];                      // "new" allocates memory for data_matrix rows * cols
      for(i=1; i<n_rows; i++) x.matrix[i] = x.matrix[0] + i*n_cols;   // ptr math: sets n_row pts to elements allocated for the rows

      x.n_rows = m;
      x.n_cols = n;
   }

   for(i=0; i<n_rows; i++)
      for(j=0; j<n_cols; j++)
      	in >> x.matrix[i][j]
         matrix[i][j] = "*";

   return in;
}
*/

//-----------------------------------------------------------------
// OVERLOADED OPERATOR:
//-----------------------------------------------------------------
ostream& operator << (ostream &out, const str_matrix &x)
{
   int i, j;
   int m = x.rows(), n = x.cols();

   for(i=0; i < m; i++)
   {
      for(j=0; j<n; j++)
         out << x(i,j) << '\t';

      out << endl;
   }

   return out;
}


