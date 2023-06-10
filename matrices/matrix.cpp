#include "matrix.h"

//------------------------------------------------------------------------------
// STATIC: static variable retains its state between different function calls and is
// shared by ALL instances of this class. This make it convenient to keep count
// of the number of calls made by model obj to a matrix obj for calculation of
// transition probabilities. This one is PUBLIC so that it can be reset by
// optimizer obj via a method in opt_tools
//------------------------------------------------------------------------------
int matrix::num_Pt = 0;

/******************* Implementation of matrix ********************/

//-----------------------------------------------------------------
// CONSTRUCTOR:
//-----------------------------------------------------------------
matrix::matrix(int m, int n, double s)
{
	int i;
	n_rows = m;
	n_cols = n;
	data_ = new double[m*n];
	for(i=0; i<m*n ; i++) data_[i] = s;
}


//-----------------------------------------------------------------
// CONSTRUCTOR: override default copy constructor
//-----------------------------------------------------------------
matrix::matrix(const matrix& x)
{
	n_rows = x.n_rows;
	n_cols = x.n_cols;
	data_ = new double[n_rows * n_cols];
	for(int i=0; i< n_rows * n_cols ; i++)
		data_[i] = x.data_[i];
}

//-----------------------------------------------------------------
// ASSIGNEMT OPERATOR: default assignment operator will not produce
// a correct copy.
//-----------------------------------------------------------------
matrix& matrix::operator = (const matrix& x)
{
	int i;
	if (this != &x)   // watch out for x = x
	{
		if (!((n_rows == x.n_rows)&&(n_cols == x.n_cols)) ) // must call reallocate
		{
			delete [] data_;
			n_rows = x.n_rows;
			n_cols = x.n_cols;
			data_ = new double[n_rows * n_cols];
		}
		for(i=0; i<n_rows * n_cols; i++)
		{
			data_[i] = x.data_[i];
		}
	}
	// cout << "matrix::operator=:: rows = " << n_rows << " and cols = " << n_cols << "\n";
	return *this;
}

//-----------------------------------------------------------------
// DESTRUCTOR:
//-----------------------------------------------------------------
matrix::~matrix()
{
	delete [] data_;
}


/**************** matrix functions and operations *****************/

//-----------------------------------------------------------------
// OVERLOADED OPERATOR:
//-----------------------------------------------------------------
ostream& operator << (ostream &out, const matrix &x)
{
	int i, j;
	int m = x.rows(), n = x.cols();

	for(i=0; i < m; i++)
	{
		if (i == 0)
		{
			out << '\t' << '\t';
			for(j=0; j<n; j++) out << setw(11) << setiosflags(ios::fixed) << setprecision(6) << j+1 << '\t';
			out << endl;
		}
		out << setw(11) << setiosflags(ios::fixed) << setprecision(6) << i+1 << '\t';
		for(j=0; j<n; j++)
			out << setw(11) << setiosflags(ios::fixed) << setprecision(6) << x(i,j) << '\t';
		out << endl;
	}

	return out;
}

//-----------------------------------------------------------------
// OVERLOADED OPERATOR:
//-----------------------------------------------------------------
matrix operator + (const matrix &x, const matrix &y)
{
	int i, j;
	int m = x.rows(), n = x.cols();

	if(! ( (x.rows() == y.rows()) && (x.cols() == y.cols() ) ) )
	{
		cout << "matrices incompatible\n";
		exit(1);
	}

	matrix c(m,n);   // new matrix constructed here

	for(i=0; i < m; i++)
		for(j=0; j < n; j++)
			c(i,j) = x(i,j) + y(i,j);

	return c;
}

//-----------------------------------------------------------------
// OVERLOADED OPERATOR:
//-----------------------------------------------------------------
matrix operator * (const matrix &x, const matrix &y)
{
	int i, j, k;
	int m = x.rows(), n = x.cols(), p = y.cols();
	double sum;

	if( y.rows() != n )
	{
		cout << "matrices incompatible\n";
		exit(1);
	}

	matrix c(m,p);   // new matrix constructed here

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
matrix operator * (const matrix &x, const double &y)
{
	int i, j;
	int m = x.rows(), n = x.cols();
	double product;

	matrix c(m,n);  // new matrix constructed here

	// NOTES:
	// this could be source of slow peeling algorithm
	// above construction could be causing a delay!
	// also this is the most simple-minded multiplication routine!  Look up more efficient algorithms for matrix multiplication?
	// also check if checking for zero and skipping the multiply is step

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
matrix operator ^ (const matrix &x, const int k)
{
	int m = x.rows(), n = x.cols();
	matrix c(m,n);   // new matrix constructed here
	c = x;

	for(int i=0; i < k; i++) {c = c * c;}

	return c;
}

//-------------------------------------------------------------------
// OVERLOADED OPERATOR: returns the Kronecker product of two matrices
// Note that (A & B) = A*I_cA & I_rB*B = (A & I_rB)(I_cA & B).
//-------------------------------------------------------------------
matrix operator & (const matrix &A, const matrix &B)
{
	int rA = A.rows(); int cA = A.cols(); int rB = B.rows(); int cB = B.cols();

	matrix I_cA(cA, cA, 0.0); I_cA.diag(1.0);
	matrix I_rB(rB, cB, 0.0); I_rB.diag(1.0);
	matrix K1(rA*rB, cA*cB, 0.0);
	matrix K2(rA*rB, cA*cB, 0.0);

	/* Calculate A XX I_rB and store it in K1. */
	for(int i=0; i<rA; i++)
		for(int j=0; j<rB; j++)
			for(int k=0; k<cA; k++)
				K1.data_[i*cA*rB*rB + j*(cA*rB) + k*rB+j] = A.data_[i*cA + k];

	/* Calculate I_cA XX B and store it in K2. */
	for(int i=0; i<cA; i++)
		for(int j=0; j<rB; j++)
			for(int k=0; k<cB; k++)
				K2.data_[ j*(cA*cB) + i*cB+k + i*cA*cB*rB] = B.data_[j*cB + k];

	return K1*K2;
}

//-----------------------------------------------------------------
// PUBLIC: replace the diagonal of a matrix with the supplied value
//-----------------------------------------------------------------
void matrix::diag(double x)
{
	for(int i=0; i<n_rows; i++)
		for(int j=0; j<n_cols; j++)
			if(i == j)
				data_[i*n_cols + j] = x;
}

//-----------------------------------------------------------------
// PUBLIC: construct an Identity matrix
//-----------------------------------------------------------------
matrix matrix::identity()
{
	int m = n_rows, n = n_cols;
    matrix I(m,n);   // new matrix constructed here

	for(int i=0; i < m; i++)
		for(int j=0; j < n; j++)
		{
			if(i==j) I(i,j) = 1.0;
			else I(i,j) = 0;
		}
    //cout << I;
    return I;
}

void matrix::set_diagonal() // for a square matrix, set the diagonal entries so that rowd sum to 0
{
	if(n_rows != n_cols)
	{
		cout << "Trying to call matrix::set_diagonal() with a non-square matrix\n";
		exit(1);
	}

	double row_sum = 0.0;
	for(int i=0; i<n_rows; i++)
	{
		for(int j=0; j<n_rows; j++)
		{
			if (i != j) row_sum = row_sum + data_[i*n_rows+j];
		}
		row_sum = -1.0 * row_sum;
		data_[i*n_rows+i] = row_sum;
		row_sum = 0.0;
	}
}

//-----------------------------------------------------------------
// PUBLIC: transpose a matrix
//-----------------------------------------------------------------
matrix matrix::t()
{
	matrix tmp(n_cols, n_rows);
	for (int i=0; i<n_rows; i++)
		for (int j=0; j<n_cols; j++)
			tmp(j,i) = data_[i*n_rows+j];

	return tmp;
}

//-----------------------------------------------------------------
// PUBLIC: delete old memory allocation and allocate new matrix
//-----------------------------------------------------------------
void matrix::reallocate(int m, int n)
{
	// Note: This is dangerous for derived types "DNA_matrix", etc.
	//       Write some code to check the matrix type and stop
	//       reallocation if this obj is NOT a base type matrix!!!

	//cout << "\n\n *************** void matrix::reallocate(int m, int n) ****************\n\n";

	// Delete the old allocation; if you do not the destructor will
        // only delete one "data_" and the other new one below represents
 	// a memory leak
	delete [] data_;
	int i;
	n_rows = m;                           // set new number of rows
	n_cols = n;                           // set new number of columns
	data_ = new double[m*n];              // reallocate memory for new matrix sized m x n
	// JRM: There is probably a good reason to set this to 0.0 and not -1??
	for(i=0; i<m*n ; i++) data_[i] = 0.0;  // set values of matrix to 0

	//cout << "\n\n******  m=" << m << " n=" << n << " *****\n\n";
}


//-----------------------------------------------------------------
// PUBLIC:
//-----------------------------------------------------------------
matrix matrix::exponentiate(const double &t, int k)
{
// Compute the matrix exponential by repeated matrix squaring.
// We want P(t) = e^Qt
// We approximate e^QT with P(t) ~ (I + Qt/m + (Qt/m)^2/2)^k,  with m=2^k
// Q is the inst rate matrix
// P is the transition probability matrix for time t
// t is the branch length
// k is the number of times the we square the first three terms of the Taylor expansion
//
// This approach is not as accurate, or as fast, as matrix exponentiation by
// using eigenvectros and eigenvalues. It is, however, applicable to any rate
// matrix, including non-reversible models.

// NOTE: A more accurate method is available by using the "eigen_GTR" method
// in the eigan namespace; the method "set_GTRb" in DNA_matrix uses the eigan
// the eigen solution of the rate matrix Q; the method "set_GTRa" in DNA_matrix
// uses the "exponentiate" method coded here.

	double m = pow(2,k);
	double scalar = t/m;
	matrix Q = *this;         // The matrix object that calls this function on will be Q!
	matrix P(n_rows,n_cols);  // The new matrix P will be the transition prob matrix
	matrix I = P.identity();  // Identity matrix of same dimension as P (and Q for that matter)

	P = I + (Q * scalar) + ((Q * scalar)^2)*0.5;  // This expression relies on code for overloaded operators
	//P = I + (Q * scalar);
	P = P ^ k;

	return P;
}


//-----------------------------------------------------------------
// PUBLIC: moved here to try minimize matrix deep copies and speed
//         up the likelihood calculation
//-----------------------------------------------------------------
void matrix::set_transition_prob(int type, double pi[], double t)
{
    num_Pt++;

	switch(type)
	{
	case 1:/*DNA*/
		if (!((n_rows == 4)&&(n_cols == 4)) )
		{
			cout << "matrix::set_transition_prob:This matrix is not a 4x4 codon matrix!\n";
			exit(1);
		}
		else
		{
			eigen::eigenQ04(data_, pi, t);  //data_ is an array and is always passed as a ptr (hence, data_ gets treated as a reference)
		}
		break;

	case 2:/*CODONS*/
		if (!((n_rows == 60)&&(n_cols == 60)) && !((n_rows == 61)&&(n_cols == 61)) )
		{
			cout << "ERROR::matrix::set_transition_prob:This matrix is not a 60x60 or 61x61 codon matrix!\n";
			exit(1);
		}
		else
		{
			//cout << "\n\n*** NUMBER OF ROWS: " << n_rows << "****\n";
			if (n_rows == 60)
				eigen::eigenQ60(data_, pi, t);  //data_ is an array and is always passed as a ptr (hence, data_ gets treated as a reference)
			if (n_rows == 61)
				eigen::eigenQ61(data_, pi, t);  //data_ is an array and is always passed as a ptr (hence, data_ gets treated as a reference)
		}
		break;

	case 3:/*AA*/
		cout << "Base class matrix does not support AA TP at this time.\n";
		break;
	}
}

// JRM
void matrix::set_transition_prob(double pi[], double t, int dim)
{
	eigen::eigenQn(data_, pi, t, dim);
}

// JRM
// void matrix::set_transition_prob(double pi[], double t, int dim, double sqrt_pi[], double lambda[], double QBV[], double W[])
// {
// 	int status;

// 	// Calculate B=pi^(1/2)*Q*pi^(-1/2) and store it in QBV
// 	for(int i=0; i<dim; i++)
// 		for(int j=0; j<i; j++)
// 			QBV[i*dim+j] = QBV[j*dim+i] = data_[i*dim+j]*sqrt_pi[i]/sqrt_pi[j]; // JRM TODO THIS LINE IS THE PROBLEM

// 	// Get B's (same as Q's) eigenvectors and store them in columns of B (Turn B into V) and put the eigenvalues in lambda.  W is just working space.
// 	status = eigen::eigenRealSym(QBV, dim, lambda, W);

// 	for(int i=0; i<dim; i++)
// 		for(int j=0; j<dim; j++)
// 			W[i*dim+j] = sqrt_pi[i] * QBV[i*dim+j];

// 	/* Calculate P. */
// 	eigen::compute_Pt(t, dim, pi, W, lambda, 0, data_, NULL, NULL);

// }

/*
//-----------------------------------------------------------------
// OVERLOADED OPERATOR: to read matrix from a file
//-----------------------------------------------------------------
istream& operator >> (istream &in, matrix &x)
{
int m, n, i, j;

in >> m >> n;

// check whether the sizes are compatible
if(!( (m == x.n_rows) && (n == x.n_cols) ) )
{
delete [] x.data_;
x.data_ = new double[m*n];
x.n_rows = m;
x.n_cols = n;
}

for(i=0; i < m*n; i++)
in >> x.data_[i];  // read by rows

return in;
}
*/

