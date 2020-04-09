// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#ifndef included_SteadyStateTemperatureStrategy
#define included_SteadyStateTemperatureStrategy

#include "TemperatureStrategy.h"
#include "HeatCapacityStrategy.h"

#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/tbox/Database.h"

#include <cassert>
#include <vector>

// This strategy fills the temperature field from the solution of
// diffusion equation
#include "TemperatureFACSolver.h"

class SteadyStateTemperatureStrategy : public TemperatureStrategy
{
 public:
   SteadyStateTemperatureStrategy(
       const int temperature_scratch_id,
       const int rhs_id,  // used internally only, but allocated outside class
       const int weight_id,
       boost::shared_ptr<tbox::Database> temperature_sys_solver_database,
       solv::LocationIndexRobinBcCoefs* bc_coefs);

   virtual ~SteadyStateTemperatureStrategy();

   void initialize(
       const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy);

   virtual double getCurrentMinTemperature(
       boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time)
   {
      (void)time;
      math::HierarchyCellDataOpsReal<double> cellops(patch_hierarchy);

      return cellops.min(d_temperature_scratch_id);
   }

   virtual double getCurrentMaxTemperature(
       boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time)
   {
      (void)time;
      math::HierarchyCellDataOpsReal<double> cellops(patch_hierarchy);

      return cellops.max(d_temperature_scratch_id);
   }

   virtual double getCurrentAverageTemperature(
       boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time)
   {
      (void)time;
      math::HierarchyCellDataOpsReal<double> cellops(patch_hierarchy);

      return cellops.integral(d_temperature_scratch_id, d_weight_id) /
             cellops.sumControlVolumes(d_temperature_scratch_id, d_weight_id);
   }

   virtual void setCurrentTemperature(
       boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time) = 0;

   void resetSolversState(
       const boost::shared_ptr<hier::PatchHierarchy> hierarchy);

   void solveSystem()
   {
      d_temperature_sys_solver->solveSystem(d_temperature_scratch_id, d_rhs_id,
                                            d_weight_id);
   }

 protected:
   TemperatureFACSolver* d_temperature_sys_solver;
   int d_temperature_scratch_id;
   int d_rhs_id;

 private:
   int d_weight_id;

   double getCurrentTemperature(const double time);
};


#endif
