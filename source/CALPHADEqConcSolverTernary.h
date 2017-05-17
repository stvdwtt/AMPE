#ifndef included_CALPHADEqConcSolverTernary
#define included_CALPHADEqConcSolverTernary

#include "DampedNewtonSolver.h"

class CALPHADEqConcentrationSolverTernary :
   public DampedNewtonSolver
{
public :

   CALPHADEqConcentrationSolverTernary(){};
      
   virtual ~CALPHADEqConcentrationSolverTernary() {};
      
   int ComputeConcentration(
      double* const conc,
      const double RTinv,
      const double* const L_AB_L,
      const double* const L_AC_L,
      const double* const L_BC_L,
      const double* const L_AB_A,
      const double* const L_AC_A,
      const double* const L_BC_A,
      const double* const fA,
      const double* const fB,
      const double* const fC );
   
protected :

   virtual void RHS(
      const double* const x,
      double* const fvec );

   virtual void Jacobian(
      const double* const x,
      double** const fjac );

   double d_RTinv;
   double d_RT;
   double d_c0;
   double d_hphi;

   //energies of 3 species, in two phase each
   double d_fA[2];
   double d_fB[2];
   double d_fC[2];
   
   // L coefficients for 2 possible phases (L and A)
   double d_L_AB_L[4];
   double d_L_AC_L[4];
   double d_L_BC_L[4];

   double d_L_AB_A[4];
   double d_L_AC_A[4];
   double d_L_BC_A[4];

};

#endif