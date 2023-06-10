#ifndef EIGEN_H_
#define EIGEN_H_

// NOTE: Don't forget that arrays get implicitly converted to a pointer, so ALL arrays
//       can/will be modified by function.  Only variables t and n get passed by value.

namespace eigen
{
   /// interface (legacy of implementations at different times) ///
   void eigenQn(double Q[], double pi[], double t, int n);
   void eigenQ04(double Q[], double pi[], double t);
   void eigenQ20(double Q[], double pi[], double t);
   void eigenQ60(double Q[], double pi[], double t);
   void eigenQ61(double Q[], double pi[], double t);

   /// wrappers ///
   void expQt(double Q[], double pi[], double t, int n, double space[]);
   void eigen_Q_decomp(double Q[], double pi[], int n, double U[], double Lambda[], double V[], double sqrt_pi[]);
   void compute_Pt(double P[], double pi[], double t, int n, double U[], double Lambda[], double V[], double exp_lam_t[]);
   void compute_Pt2(double P[], double P_dp[], double P_d2p[], double pi[], double t, int n, double U[], double Lambda[], double V[], double exp_lam_t[]);

   /// eigen solution for real symmetric matrix ///
   int eigenRealSym(double A[], int n, double Root[], double work[]);
   void HouseholderRealSym(double a[], int n, double d[], double e[]);
   int EigenTridagQLImplicit(double d[], double e[], int n, double z[]);
   void EigenSort(double d[], double U[], int n);

   /// tools ///
   #define FOR(i,n) for(i=0; i<n; i++)
   int identity (double x[], int n);
   int zero (double x[], int n);

   /* Complete transition matrix and optionally, 1st and 2nd derivatives */
   /* Adopted from code from Dr. Ed. Susko.  Added by JRM for covarion models. */
   void compute_Pt(double t, int Q_dim, double *PiVec, double *W, double *Lam,
      int deriv, double *P, double *Pp, double *Ppp);
}
#endif /* EIGEN_H_ */
