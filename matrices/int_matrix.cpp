#include "int_matrix.h"

/******************* Implementation of matrix ********************/

//-----------------------------------------------------------------
// CONSTRUCTOR:
//-----------------------------------------------------------------
int_matrix::int_matrix(int m, int n, int s)
{
   int i;
   n_rows = m;
   n_cols = n;
   data_ = new int[m*n];
   for(i=0; i<m*n ; i++) data_[i] = s;
}


//-----------------------------------------------------------------
// CONSTRUCTOR: override default copy constructor
//-----------------------------------------------------------------
int_matrix::int_matrix(const int_matrix& x)
{
   int i;
   n_rows = x.n_rows;
   n_cols = x.n_cols;
   data_ = new int[n_rows * n_cols];
   for(i=0; i< n_rows * n_cols ; i++)
      data_[i] = x.data_[i];
}

//-----------------------------------------------------------------
// ASSIGNEMT OPERATOR: default assignment operator will not produce
// a correct copy.
//-----------------------------------------------------------------
int_matrix& int_matrix::operator = (const int_matrix& x)
{
   int i;
   if (this != &x)   // watch out for x = x
   {
      if (!((n_rows == x.n_rows)&&(n_cols == x.n_cols)) ) // must call reallocate
      {
         delete [] data_;
         n_rows = x.n_rows;
         n_cols = x.n_cols;
         data_ = new int[n_rows * n_cols];
      }
      for(i=0; i<n_rows * n_cols; i++)
            data_[i] = x.data_[i];
   }
   return *this;
}


//-----------------------------------------------------------------
// PUBLIC: delete old memory allocation and allocate new matrix
//-----------------------------------------------------------------
void int_matrix::reallocate(int m, int n)
{
   delete [] data_;                      // delete the old allocation
   int i;
   n_rows = m;                           // set new number of rows
   n_cols = n;                           // set new number of columns
   data_ = new int[m*n];                 // reallocate memory for new matrix sized m x n
   for(i=0; i<m*n ; i++) data_[i] = 0;   // set values of matrix to 0
}


//-----------------------------------------------------------------
// DESTRUCTOR:
//-----------------------------------------------------------------
int_matrix::~int_matrix(void)
{
   delete [] data_;
}


/**************** matrix functions and operations *****************/

//-----------------------------------------------------------------
// OVERLOADED OPERATOR:
//-----------------------------------------------------------------
ostream& operator << (ostream &out, const int_matrix &x)
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

//-----------------------------------------------------------------
// OVERLOADED OPERATOR:
//-----------------------------------------------------------------
int_matrix operator + (const int_matrix &x, const int_matrix &y)
{
   int i, j;
   int m = x.rows(), n = x.cols();

   if(! ( (x.rows() == y.rows()) && (x.cols() == y.cols() ) ) )
   {
      cout << "matrices incompatible\n";
      exit(1);
   }

   int_matrix c(m,n);   // new matrix constructed here

   for(i=0; i < m; i++)
      for(j=0; j < n; j++)
         c(i,j) = x(i,j) + y(i,j);

   return c;
}

//-----------------------------------------------------------------
// OVERLOADED OPERATOR:
//-----------------------------------------------------------------
int_matrix operator * (const int_matrix &x, const int_matrix &y)
{
   int i, j, k;
   int m = x.rows(), n = x.cols(), p = y.cols();
   int sum;

   if( y.rows() != n )
   {
      cout << "matrices incompatible\n";
      exit(1);
   }

   int_matrix c(m,p);   // new matrix constructed here

   for(i=0; i < m; i++)
      for(j=0; j < p; j++)
      {
         sum = 0.0f;
         for(k=0; k<n; k++)
            sum = sum + x(i,k) * y(k,j);

         c(i,j) = sum;
      }

   return c;
}

//-----------------------------------------------------------------
// OVERLOADED OPERATOR:
//-----------------------------------------------------------------
int_matrix operator * (const int_matrix &x, const double &y)
{
   int i, j;
   int m = x.rows(), n = x.cols();
   double product;

   int_matrix c(m,n);  // new matrix constructed here

   for(i=0; i < m; i++)
      for(j=0; j < n; j++)
      {
         product = x(i,j) * y;
         c(i,j) = product;
      }

   return c;
}

//-----------------------------------------------------------------
// OVERLOADED OPERATOR:
//-----------------------------------------------------------------
int_matrix operator ^ (const int_matrix &x, const int k)
{
	int m = x.rows(), n = x.cols();
	int_matrix c(m,n);   // new matrix constructed here
	c = x;

	for(int i=0; i < k; i++) {c = c * c;}

	return c;
}

//-----------------------------------------------------------------
// PUBLIC: construct an Identity matrix
//-----------------------------------------------------------------
int_matrix int_matrix::identity()
{
	int m = n_rows, n = n_cols;
    int_matrix I(m,n);   // new matrix constructed here

	for(int i=0; i < m; i++)
	   for(int j=0; j < n; j++)
	   {
         if(i==j) I(i,j) = 1.0;
         else I(i,j) = 0;
	   }
    //cout << I;
    return I;
}


/*
//-----------------------------------------------------------------
// OVERLOADED OPERATOR: to read matrix from a file
//-----------------------------------------------------------------
istream& operator >> (istream &in, int_matrix &x)
{
   int m, n, i, j;

   in >> m >> n;

   // check whether the sizes are compatible
   if(!( (m == x.n_rows) && (n == x.n_cols) ) )
   {
      delete [] x.data_;
      x.data_ = new int[m*n];
      x.n_rows = m;
      x.n_cols = n;
   }

   for(i=0; i < m*n; i++)
         in >> x.data_[i];  // read by rows

   return in;
}
*/

