#include "PhaseTemperatureFACOps.h"

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/tbox/Utilities.h"
#include "FuncFort.h"

#include <cassert>
using namespace std;

//======================================================================

PhaseTemperatureFACOps::PhaseTemperatureFACOps(
   const std::string &object_name,
   boost::shared_ptr<tbox::Database> database )
   :
   EllipticFACOps( object_name, database )
{
   d_C_is_set = false;
   d_D_is_set = false;
   d_M_is_set = false;

   return;
}

//======================================================================
// expects M/cp as mobility
void PhaseTemperatureFACOps::setOperatorCoefficients(
   const int phase_id,
   const int mobility_id,
   const double epsilon_phase, 
   const double latent_heat,
   const string phase_interp_func_type,
   const double phase_well_scale,
   const string phase_well_func_type)
{
   assert( mobility_id>=0 );
   assert( latent_heat>0. );

   setM( mobility_id );

   // C to be set after M since it uses M
   setC(
      phase_id,
      latent_heat,
      phase_interp_func_type,
      phase_well_scale,
      phase_well_func_type);
   
   setDConstant( -latent_heat * epsilon_phase * epsilon_phase );
}

//======================================================================

// C = (L/cp) * phi_mobility * (
//       phi_well_scale * phi_well_func'' )

void PhaseTemperatureFACOps::setC(
   const int phi_id,
   const double factor,
   const string phi_interp_func_type,
   const double phi_well_scale,
   const string phi_well_func_type)
{
   assert( phi_id >= 0 );
   assert( d_m_id >= 0 );
   assert( d_c_id >= 0 );
   assert( d_M_is_set );
  
   //tbox::pout<<"PhaseTemperatureFACOps::setC()..."<<endl;
 
   for ( int ln = d_ln_min; ln <= d_ln_max; ++ln ) {
      boost::shared_ptr< hier::PatchLevel > level =
         d_hierarchy->getPatchLevel( ln );

      for (hier::PatchLevel::iterator pi(level->begin());
           pi != level->end(); ++pi) {

         boost::shared_ptr< hier::Patch > patch = *pi;
         const hier::Box& patch_box = patch->getBox();
      
         boost::shared_ptr< pdat::CellData<double> > phi_data (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( phi_id) ) );
         
         boost::shared_ptr< pdat::CellData<double> > local_m_data (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_m_id) ) );

         boost::shared_ptr< pdat::CellData<double> > cdata (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_c_id) ) );

         setCOnPatchForPreconditionODE(
            phi_data,
            local_m_data,
            cdata,
            factor,
            phi_interp_func_type.c_str(),
            phi_well_scale,
            phi_well_func_type.c_str(),
            patch_box );
      }
   }
   
   d_poisson_spec.setCPatchDataId( d_c_id );
   d_C_is_set = true;

   return;
}

void PhaseTemperatureFACOps::setCOnPatchForPreconditionODE(
   boost::shared_ptr< pdat::CellData<double> > cd_phi,
   boost::shared_ptr< pdat::CellData<double> > cd_m,
   boost::shared_ptr< pdat::CellData<double> > cd_c,
   const double latent_heat,
   const char* phi_interp_func_type,
   const double phi_well_scale,
   const char* phi_well_func_type,
   const hier::Box& pbox )
{
   double* ptr_phi = cd_phi->getPointer();
   double* ptr_m = cd_m->getPointer();
   double* ptr_c = cd_c->getPointer();
   double* ptr_eta = NULL;
   
   const hier::Box& c_gbox = cd_c->getGhostBox();
   int imin_c = c_gbox.lower(0);
   int jmin_c = c_gbox.lower(1);
   int jp_c = c_gbox.numberCells(0);
   int kmin_c = 0;
   int kp_c = 0;
#if (NDIM == 3)
   kmin_c = c_gbox.lower(2);
   kp_c = jp_c * c_gbox.numberCells(1);
#endif

   const hier::Box& m_gbox = cd_m->getGhostBox();
   int imin_m = m_gbox.lower(0);
   int jmin_m = m_gbox.lower(1);
   int jp_m = m_gbox.numberCells(0);
   int kmin_m = 0;
   int kp_m = 0;
#if (NDIM == 3)
   kmin_m = m_gbox.lower(2);
   kp_m = jp_m * m_gbox.numberCells(1);
#endif

   const hier::Box& pf_gbox = cd_phi->getGhostBox();
   int imin_pf = pf_gbox.lower(0);
   int jmin_pf = pf_gbox.lower(1);
   int jp_pf = pf_gbox.numberCells(0);
   int kmin_pf = 0;
   int kp_pf = 0;
#if (NDIM == 3)
   kmin_pf = pf_gbox.lower(2);
   kp_pf = jp_pf * pf_gbox.numberCells(1);
#endif

   int imin = pbox.lower(0);
   int imax = pbox.upper(0);
   int jmin = pbox.lower(1);
   int jmax = pbox.upper(1);
   int kmin = 0;
   int kmax = 0;
#if (NDIM == 3)
   kmin = pbox.lower(2);
   kmax = pbox.upper(2);
#endif
         
   for ( int kk = kmin; kk <= kmax; kk++ ) {
      for ( int jj = jmin; jj <= jmax; jj++ ) {
         for ( int ii = imin; ii <= imax; ii++ ) {

            const int idx_c = (ii - imin_c) +
               (jj - jmin_c) * jp_c + (kk - kmin_c) * kp_c;

            const int idx_m = (ii - imin_m) +
               (jj - jmin_m) * jp_m + (kk - kmin_m) * kp_m;

            const int idx_pf = (ii - imin_pf) +
               (jj - jmin_pf) * jp_pf + (kk - kmin_pf) * kp_pf;

            const double m = ptr_m[idx_m];
            const double phi = ptr_phi[idx_pf];

            const double g_phi_dbl_prime =
               FORT_SECOND_DERIV_WELL_FUNC(
                  phi,
                  phi_well_func_type );

            ptr_c[idx_c] =
               latent_heat * m * phi_well_scale * g_phi_dbl_prime;

         }
      }
   }
}

//======================================================================
void PhaseTemperatureFACOps::multiplyDTDPhiBlock(
   const int phase_id,
   const int out_id)
{
   evaluateRHS(phase_id,out_id);
}
