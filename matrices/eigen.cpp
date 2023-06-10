//NOTE: for branch length optimization, I need to use U[], root[], and V[] (this is U, lambda, and U^-1)
//      -- 1st and 2nd derivatives of P(t) computed in method called "compute_Pt2"
//      -- arguments to "compute_Pt2" include P_dp[] and P_d2p[] these will contain the 1st and 2nd derivatives of P(t)
//      -- compute_likelihood_updates(t) is interface for an optimizer object that will carry out the newton branch length estimation; in CODON_model, DNA_model, and AA_model only
//      -- update_lnl_df(double P_dp[], double P_d2p[]) does the work and has model specific implementation!


#include <iostream>
#include <math.h>
#include "eigen.h"
#include <cstdlib>

using namespace std;

namespace eigen
{
/* Notes on "space" array:
 * -- size of U = n*n
 * -- size of V = n*n
 * -- size of Root = n
 * -- size of sqrt_pi = n
 * -- size of exp_lam_t = n
 * -- "space" holds U, V, Root, sqrt_pi and exp_lam_t
 * -- hence, size of "space" = n*n*2+n*3
 */

// Might add something like #define space_eigen(n) ( ((n*n*2)+(n*3))*sizeof(double) ) ???

//--------functions--------//
	void eigenQn(double Q[], double pi[], double t, int n)
	{
		double* space = new double[n*n*2+n*3];
		expQt(Q, pi, t, n, space);
		delete [] space;
	}


	void eigenQ04(double Q[], double pi[], double t)
	{
		// do not yet know if this works; TEST ME!!!
		int n = 4;
		double space[44];
		expQt(Q, pi, t, n, space);
	}


	void eigenQ20(double Q[], double pi[], double t) // UN-TESED!!!
	{
		int n = 20;
		double space[860];
		expQt(Q, pi, t, n, space);
	}

	void eigenQ60(double Q[], double pi[], double t)
	{
		//cout << "\n\n ****** INSIDE eigenQ60() *******\n\n";
		int n = 60;
		double space[7380];
		expQt(Q, pi, t, n, space);
	}


	void eigenQ61(double Q[], double pi[], double t)
	{
		//cout << "\n\n ****** INSIDE eigenQ61() *******\n\n";
		int n = 61;
		double space[7625];
		expQt(Q, pi, t, n, space);
	}

//-----------------------------------------------------------------
// exponentiate Q * t
//-----------------------------------------------------------------
	void expQt(double Q[], double pi[], double t, int n, double space[])
	{
/*  Calculates P(t) = exp(Q*t)
 *  -- Q is the rate matrix for a time-reversible Markov process.
 *  -- Q[] has the rate matrix as input
 *  -- P(t) is within Q[] in return.
 *  -- space[n*n*2+n*2]; see notes at top
 */

		//partition "space" by using pointers and pass these to the function eigen_Q_decomp
		double *U=space, *V=U+n*n, *Lambda=V+n*n, *sqrt_pi=Lambda+n, *exp_lam_t=sqrt_pi+n;   // NOTE: this is done b/c the matrix dimensions are variable; this way no need to deal with size in the function that does the decomp

		//---Find the eigen solution of the rate matrix Q for a time-reversible Markov process
		eigen_Q_decomp(Q, pi, n, U, Lambda, V, sqrt_pi);                                     // NOTE: Q-matrix = U * diag{Root} * V; see Adachi and Hasegawa 1996b


		//--- P(t) = U * exp{Root*t} * V ---//
		compute_Pt(Q, pi, t, n, U, Lambda, V, exp_lam_t);                                    // NOTE: Q will get zeroed out and set to P(t)


		/*
		 * Temp....
		 double P_dp[61*61];   //temp
		 double P_d2p[61*61];  //temp
		 compute_Pt2(Q, P_dp, P_d2p, pi, t, n, U, Lambda, V, exp_lam_t);
		*/
	}


//-----------------------------------------------------------------
//
//-----------------------------------------------------------------
    void eigen_Q_decomp(double Q[], double pi[], int n, double U[], double Lambda[], double V[], double sqrt_pi[])
    {
	//  Q = R * diag{pi} = U * diag{Lambda} * V

	//int status;
	double small=1e-100;

	//--------Check pi's, take square root, store in sqrt_pi--------//
	int counter = 0;
	for(int i=0; i<n; i++)
	{
	    // if pi is big enough, take square root and store it
	    if(pi[i]>small) sqrt_pi[counter++]=sqrt(pi[i]);  // TODO JRM, this can cause a problem if seq. data is short
	    for(int j=0; j<n; j++)  U[i*n+j] = Q[i*n+j]; // Copy Q to U
	}

	//-----Eigen decomposition of U = Pi^(1/2) * Q * Pi^(-1/2)-------//
	if(counter == n)
	{
	    for(int i=0; i<n; i++)
		for(int j=0; j<i; j++)
		    // Set U = Pi^(1/2) * Q * Pi^(-1/2)
		    U[i*n+j] = U[j*n+i] = (Q[i*n+j] * sqrt_pi[i]/sqrt_pi[j]);

	    // Find the eigen solution for real symmetrical matrix U = Pi^(1/2) * Q * Pi^(-1/2)
	    // Result: (1) Rows of U will = eigenvectors and (2) Lambda is diagonal matrix having the eigenvalues
	    //status=eigenRealSym(U, n, Lambda, V);
	    eigenRealSym(U, n, Lambda, V);

	    // construct V = U * Pi^(1/2)
	    for(int i=0;i<n;i++) for(int j=0;j<n;j++)  V[i*n+j] = U[j*n+i] * sqrt_pi[j];

	    // construct U = U / Pi^(1/2)
	    for(int i=0;i<n;i++) for(int j=0;j<n;j++)  U[i*n+j] /= sqrt_pi[i];
	}
	else
	{
	    cout << "\nERROR::eigen_decomp_Q::counter != n\n";  // deal with this later if it ever comes up
	    cout << "counter = " << counter << "  and  n = " << n << "\n";
	    cout << "Some equilibrium frequencies may be too small. Try using the F3x4 model.\n";
	    cout << "Sorry, Proteus must terminate.\n";
	    exit(1);

      /* Note: Below is Ziheng Yang's code; this might offer a solution!
      int i,j, inew, jnew, nnew, status;
      double small=1e-100;

      for(i=0, inew=0; i<n; i++)
      {
         if(pi[i]>small)
         {
            for(j=0,jnew=0; j<i; j++)
            {
               if(pi[j]>small)
               {
                  U[inew*nnew+jnew] = U[jnew*nnew+inew]
                                    = Q[i*n+j] * pi_sqrt[inew]/pi_sqrt[jnew];
                  jnew++;
               }
            }
            U[inew*nnew+inew] = Q[i*n+i];
            inew++;
         }
      }

      status=eigenRealSym(U, nnew, Root, V);


      // construct Root //
      for(i=n-1,inew=nnew-1; i>=0; i--)
      {
         Root[i] = (pi[i]>small ? Root[inew--] : 0);
      }


      // construct V //
      for(i=n-1,inew=nnew-1; i>=0; i--)
      {
         if(pi[i]>small)
         {
            for(j=n-1,jnew=nnew-1; j>=0; j--)
               if(pi[j]>small)
               {
                  V[i*n+j] = U[jnew*nnew+inew]*pi_sqrt[jnew];
                  jnew--;
               }
               else
                  V[i*n+j] = (i==j);
            inew--;
         }
         else
         {
            for(j=0; j<n; j++)  V[i*n+j] = (i==j);
         }
      }


      // construct U //
      for(i=n-1,inew=nnew-1; i>=0; i--)
      {
         if(pi[i]>small) {
            for(j=n-1,jnew=nnew-1;j>=0;j--)
               if(pi[j]>small) {
                  U[i*n+j] = U[inew*nnew+jnew]/pi_sqrt[inew];
                  jnew--;
               }
               else
                  U[i*n+j] = (i==j);
            inew--;
         }
         else
            for(j=0;j<n;j++)
               U[i*n+j] = (i==j);
      }
      */		    
	}
    }


//-----------------------------------------------------------------
// Note: There are two compute_Pt methods with different signatures
// JPB
//-----------------------------------------------------------------
	void compute_Pt(double P[], double pi[], double t, int n, double U[], double Lambda[], double V[], double exp_lam_t[])
	{
// P(t) = U * exp{Lamda*t} * V
// Lambda is passed as argument and has eigenvalues (lambda) of Q
// Lambda = diag{lambda1, lambda2, ... , lamda_n}
// exp(Lambda,t) = diag{exp(lambda1,t), exp(lambda2.t), ... , exp(lambda,t)}
// columns of U = right eigenvectors of Q
// rows of V = left eigenvectors of Q
// P(t) = U * exp(Root,t) * U^-1
// P(t) = U * exp(Root,t) * V
// n = 4, 20 or 61, or multiples of those for covarion-like models

		int i,j,k;
		double sum;
		double lower_bound = 0;

		if (t<-0.1) cout << "\nCAUTION::eigen::compute_Pt: t = " << t << "\n";
		if (t<1e-100) { identity (P, n); return; }  // if t is very very small, set P = identity matrix and exit function

		for(k = 0; k < n; k++) exp_lam_t[k] = exp(t * Lambda[k]);

		zero(P,n*n);
		i=0; j=0; k=0;

		// This computes P(t) from diagonalized Q

		//---  U * exp(Root,t)  ---//
		for(i=0; i<n; i++)
		{
			for(j=0; j<n; j++)
			{
				U[j*n+i] = U[j*n+i] * exp_lam_t[i];
			}
		}

		//--- P[] = U * exp(Root,t) * V  ---//
		for(i=0; i<n; i++)
		{
			for(j=0; j<n; j++)
			{
				sum = 0.0;
				for(k=0; k<n; k++ )
				{
					sum += U[i*n+k] * V[k*n+j]; // accumulate element-wise mult of U(i,k) * V(k,j)
				}
				P[i*n+j] = sum;                // sum = P(i,j)
			}
		}

		/*
		  for(int l=0; l<n*n; l++)
		  {
		  cout << "\n" << P[l] << " ";
		  }
		  cout << "\n\n";
		*/

		for(i=0;i<n*n;i++)  if(P[i]<lower_bound)  P[i]=0;
	}


//-----------------------------------------------------------------
//
//-----------------------------------------------------------------
	void compute_Pt2(double P[], double P_dp[], double P_d2p[], double pi[], double t, int n, double U[], double Lambda[], double V[], double exp_lam_t[])
	{
//NOTE: This is a potentially faster version of compute_Pt that also
//      computes the fist and second derivatives of Pt (P_dp and P_d2p).
//      Code is based on algorithm of Ed Susko.

		int i,j,k;
		double lower_bound = 0;
		double w,pt,dpt,d2pt;

		if (t<-0.1) cout << "\nCAUTION::eigen::compute_Pt: t = " << t << "\n";
		if (t<1e-100) { identity (P, n); return; }

		i=0; j=0; k=0; zero(P,n*n);
		for(k = 0; k < n; k++) exp_lam_t[k] = exp(t * Lambda[k]);

		// This computes P(t) from diagonalized Q
		for(i = 0; i < n; i++)
		{
			for(j = 0; j < n; j++)
			{
				pt = 0.0; dpt = 0.0; d2pt = 0.0;
				for(k = 0; k < n; k++)
				{
					w = V[i+k*n] * V[j+k*n] * exp_lam_t[k];
					//w = U[j*n+k] * exp_lam_t[k] * V[k*n+j];   BROKE
					pt += w;
					dpt += Lambda[k]*w;
					d2pt += Lambda[k]*Lambda[k]*w;
				}
				P[i*n+j] = pt/pi[i];       // This works
				P_dp[i*n+j] = dpt/pi[i];   // Bug here? Check in Mathmatica, or against Joey's code.
				P_d2p[i*n+j] = d2pt/pi[i]; // Bug here? Check in Mathmatica, or against Joey's code.
			}
		}

		/*
		  for(int l=0; l<n*n; l++)
		  {
		  cout << P[l] << " ";
		  }
		  cout << "\n\n";

		  for(int l=0; l<n*n; l++)
		  {
		  cout << P_dp[l] << " ";
		  }
		  cout << "\n\n";

		  for(int l=0; l<n*n; l++)
		  {
		  cout << P_d2p[l] << " ";
		  }
		  cout << "\n\n";
		*/


		for(i=0;i<n*n;i++)  if(P[i]<lower_bound)  P[i]=0;
	}

/* Complete transition matrix and optionally, 1st and 2nd derivatives. */
/* Adopted from code from Dr. Susko.  Added by JRM for covarion models. */
// Note: There are two compute_Pt methods with different signatures
    void compute_Pt(double t, int Q_dim, double* PiVec, double* W, double* Lam,
	    int deriv, double* P, double* Pp, double* Ppp)
    {
	/* IN:
	 *  t - edge length
	 *  m - number of possible character states
	 *  g - number of rate classes for the covarion-like model
	 *  PiVec - m dimensional. stationary frequencies for the model
	 *  W, Lam - m*2*g x m*2*g and m*2*g dimensional. Giving W and diagonal
	 *           entries of Lam in the conventional decompostion:
	 *                    P(t) = Pi**(-1) W Lam W'
	 *  deriv - 0, 1 or 2. The number of derivatives required
	 *
	 * OUT:
	 *  P - 2*g*m x 2*g*m. P[i+j*m] gives P(t)_ij
	 *  Pp - (not used if deriv = 0). m dimensional. Pp[i+j*m] gives P'(t)_ij
	 *  Ppp - (not used if deriv <= 1). m dimensional.
	 *       Ppp[i+j*m] gives P''(t)_ij
	 */
	int i,j,k;
	double w;

	for(i = 0; i < Q_dim; i++)
	{
	    for(j = 0; j < Q_dim; j++)
	    {
		P[i*Q_dim + j] = 0.0;
		if(deriv==1 || deriv==2) Pp[i*Q_dim+ j] = 0.0;
		if(deriv==2) Ppp[i*Q_dim + j] = 0.0;

		for(k = 0; k < Q_dim; k++)
		{
		    /* W %*% e^(lt) %*% W^T */
		    w = W[i*Q_dim + k] * W[j*Q_dim + k] * exp(Lam[k]*t);
		    P[i*Q_dim + j] += w;

		    if(deriv==1 || deriv==2) Pp[i*Q_dim + j] += Lam[k]*w;
		    if(deriv==2) Ppp[i*Q_dim + j] += Lam[k]*Lam[k]*w;
		}
		P[i*Q_dim + j] /= PiVec[i];
		if(deriv==1 || deriv==2) Pp[i*Q_dim + j] /= PiVec[i];
		if(deriv==2) Ppp[i*Q_dim + j] /= PiVec[i];
	    }
	}
    }

//-----------------------------------------------------------------
// create an identity matrix
//-----------------------------------------------------------------
	int identity (double x[], int n)
	{ int i,j;  FOR (i,n)  { FOR(j,n)  x[i*n+j]=0;  x[i*n+i]=1; }  return (0); }

//-----------------------------------------------------------------
// zero out the elements of a matrix
//-----------------------------------------------------------------
	int zero (double x[], int n)
	{ int i; FOR (i,n) x[i]=0; return (0);}


//-----------------------------------------------------------------
//
//-----------------------------------------------------------------
	int eigenRealSym(double A[], int n, double Root[], double work[])
	{
/* This finds the eigen solution of a real symmetrical matrix A[n*n].  In return,
   A has the right vectors and Root has the eigenvalues.
   work[n] is the working space.
   The matrix is first reduced to a tridiagonal matrix using HouseholderRealSym(),
   and then using the QL algorithm with implicit shifts.

   Adapted from routine tqli in Numerical Recipes in C, with reference to LAPACK
   Ziheng Yang, 23 May 2001
*/

		int status=0;
		HouseholderRealSym(A, n, Root, work);                // reduce a real symmetric matrix to tri-diagonal matrix
		status=EigenTridagQLImplicit(Root, work, n, A);      // find the eigen solution of a tridiagonal matrix represented by root and work
		EigenSort(Root, A, n);                               // put values in proper order
		return(status);                                      // if went well: "A" has rt eigen vectors & "Root" has eigen values
	}

//-----------------------------------------------------------------
//
//-----------------------------------------------------------------
	void EigenSort(double d[], double U[], int n)
	{
/* this sorts the eigen values d[] and rearrange the (right) eigen vectors U[]
 */
		int k,j,i;
		double p;

		for (i=0;i<n-1;i++) {
			p=d[k=i];
			for (j=i+1;j<n;j++)
				if (d[j] >= p) p=d[k=j];
			if (k != i) {
				d[k]=d[i];
				d[i]=p;
				for (j=0;j<n;j++) {
					p=U[j*n+i];
					U[j*n+i]=U[j*n+k];
					U[j*n+k]=p;
				}
			}
		}
	}

//-----------------------------------------------------------------
//
//-----------------------------------------------------------------
	void HouseholderRealSym(double a[], int n, double d[], double e[])
	{
/* This uses HouseholderRealSym transformation to reduce a real symmetrical matrix
   a[n*n] into a tridiagonal matrix represented by d and e.
   d[] is the diagonal (eigends), and e[] the off-diagonal.
*/
		int m,k,j,i;
		double scale,hh,h,g,f;

		for (i=n-1;i>=1;i--) {
			m=i-1;
			h=scale=0;
			if (m > 0) {
				for (k=0;k<=m;k++)
					scale += fabs(a[i*n+k]);
				if (scale == 0)
					e[i]=a[i*n+m];
				else {
					for (k=0;k<=m;k++) {
						a[i*n+k] /= scale;
						h += a[i*n+k]*a[i*n+k];
					}
					f=a[i*n+m];
					g=(f >= 0 ? -sqrt(h) : sqrt(h));
					e[i]=scale*g;
					h -= f*g;
					a[i*n+m]=f-g;
					f=0;
					for (j=0;j<=m;j++) {
						a[j*n+i]=a[i*n+j]/h;
						g=0;
						for (k=0;k<=j;k++)
							g += a[j*n+k]*a[i*n+k];
						for (k=j+1;k<=m;k++)
							g += a[k*n+j]*a[i*n+k];
						e[j]=g/h;
						f += e[j]*a[i*n+j];
					}
					hh=f/(h*2);
					for (j=0;j<=m;j++) {
						f=a[i*n+j];
						e[j]=g=e[j]-hh*f;
						for (k=0;k<=j;k++)
							a[j*n+k] -= (f*e[k]+g*a[i*n+k]);
					}
				}
			}
			else
				e[i]=a[i*n+m];
			d[i]=h;
		}
		d[0]=e[0]=0;

		/* Get eigenvectors */
		for (i=0;i<n;i++) {
			m=i-1;
			if (d[i]) {
				for (j=0;j<=m;j++) {
					g=0;
					for (k=0;k<=m;k++)
						g += a[i*n+k]*a[k*n+j];
					for (k=0;k<=m;k++)
						a[k*n+j] -= g*a[k*n+i];
				}
			}
			d[i]=a[i*n+i];
			a[i*n+i]=1;
			for (j=0;j<=m;j++) a[j*n+i]=a[i*n+j]=0;
		}
	}

//-----------------------------------------------------------------
//
//-----------------------------------------------------------------
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
	int EigenTridagQLImplicit(double d[], double e[], int n, double z[])
	{
/* This finds the eigen solution of a tridiagonal matrix represented by d and e.
   d[] is the diagonal (eigenvalues), e[] is the off-diagonal
   z[n*n]: as input should have the identity matrix to get the eigen solution of the
   tridiagonal matrix, or the output from HouseholderRealSym() to get the
   eigen solution to the original real symmetric matrix.
   z[n*n]: has the orthogonal matrix as output

   Adapted from routine tqli in Numerical Recipes in C, with reference to
   LAPACK fortran code.
   Ziheng Yang, May 2001
*/
		int m,j,iter,niter=30, status=0, i,k;
		double s,r,p,g,f,dd,c,b, aa,bb;

		for (i=1;i<n;i++) e[i-1]=e[i];  e[n-1]=0;
		for (j=0;j<n;j++) {
			iter=0;
			do {
				for (m=j;m<n-1;m++) {
					dd=fabs(d[m])+fabs(d[m+1]);
					if (fabs(e[m])+dd == dd) break;  /* ??? */
				}
				if (m != j) {
					if (iter++ == niter) {
						status=-1;
						break;
					}
					g=(d[j+1]-d[j])/(2*e[j]);

					/* r=pythag(g,1); */

					if((aa=fabs(g))>1)  r=aa*sqrt(1+1/(g*g));
					else                r=sqrt(1+g*g);

					g=d[m]-d[j]+e[j]/(g+SIGN(r,g));
					s=c=1;
					p=0;
					for (i=m-1;i>=j;i--) {
						f=s*e[i];
						b=c*e[i];

						/*  r=pythag(f,g);  */
						aa=fabs(f); bb=fabs(g);
						if(aa>bb)       { bb/=aa;  r=aa*sqrt(1+bb*bb); }
						else if(bb==0)             r=0;
						else            { aa/=bb;  r=bb*sqrt(1+aa*aa); }

						e[i+1]=r;
						if (r == 0) {
							d[i+1] -= p;
							e[m]=0;
							break;
						}
						s=f/r;
						c=g/r;
						g=d[i+1]-p;
						r=(d[i]-g)*s+2*c*b;
						d[i+1]=g+(p=s*r);
						g=c*r-b;
						for (k=0;k<n;k++) {
							f=z[k*n+i+1];
							z[k*n+i+1]=s*z[k*n+i]+c*f;
							z[k*n+i]=c*z[k*n+i]-s*f;
						}
					}
					if (r == 0 && i >= j) continue;
					d[j]-=p; e[j]=g; e[m]=0;
				}
			} while (m != j);
		}
		return(status);
	}
#undef SIGN

//end namespace eigen//
}
