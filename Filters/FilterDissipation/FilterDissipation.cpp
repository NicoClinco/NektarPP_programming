///////////////////////////////////////////////////////////////////////////////
//
// File: FilterDissipation.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Output The balance of the kinetic energy.
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <SolverUtils/Filters/FilterDissipation.h>

using namespace std;

namespace Nektar
{
  namespace SolverUtils
  {
    std::string FilterDissipation::className =
      SolverUtils::GetFilterFactory().RegisterCreatorFunction("Dissipation",
							      FilterDissipation::create);

    FilterDissipation::FilterDissipation(const LibUtilities::SessionReaderSharedPtr &pSession,
					 const std::weak_ptr<EquationSystem> &pEquation,
					 const ParamMap &pParams)
      : Filter(pSession, pEquation), m_index(-1), m_homogeneous(false), m_planes()
    {
      // OutputFile
      auto it = pParams.find("OutputFile");
      if (it == pParams.end())
	{
	  m_outputFile = m_session->GetSessionName();
	}
      else
	{
	  ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
	  m_outputFile = it->second;
	}
      m_outputFile += ".avg";

      // OutputFrequency
      it = pParams.find("OutputFrequency");
      ASSERTL0(it != pParams.end(), "Missing parameter 'OutputFrequency'.");
      LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
      m_outputFrequency = round(equ.Evaluate());

      
    }// end constructor

    FilterDissipation::~FilterDissipation()
    {
    } //end virtual destructor

    void FilterDissipation::v_Initialise(
					 const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
					 const NekDouble &time)
    {
      MultiRegions::ExpListSharedPtr areaField;
    
      // Lock equation system pointer
      auto equ = m_equ.lock();
      ASSERTL0(equ, "Weak pointer expired");

   
      // Open OutputFile
      LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
      if (vComm->GetRank() == 0)
	{
	  m_outputStream.open(m_outputFile.c_str());
	  ASSERTL0(m_outputStream.good(),
		   "Unable to open: '" + m_outputFile + "'");
	  m_outputStream.setf(ios::scientific, ios::floatfield);
	  m_outputStream << "# Time,epsilon" << "\n";
 
	}

      // Output values at initial time.
      m_index = 0;
      v_Update(pFields, time);
    } //end initialization

    void FilterDissipation::v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
				     const NekDouble &time)
    {
      // Only output every m_outputFrequency
      if ((m_index++) % m_outputFrequency)
	{
	  return;
	}

      auto equ = m_equ.lock();
      ASSERTL0(equ,"Weak pointer expired");

      //Ensure that we are using a fluid-interface as system:
      auto fluidEqu = std::dynamic_pointer_cast<FluidInterface>(equ);
      ASSERTL0(fluidEqu,"The filter is not compatible with this solver");


      //Check the dimensions of the problem:
      unsigned int dim =0;
      usinged int nPoints = pFields[0]->GetNpoints();

      if(pFields[0]->GetExpType() == MultiRegions::e1D)
	{
	  dim=1;
	}
      if(pFields[0]->GetExpType() == MultiRegions::e2D)
	{
	  dim=2;
	}
      if(pFields[0]->GetExpType() == MultiRegions::e3D)
	{
	  dim=3;
	}
      
      //Get the thermodynamc variables, density,pressure,velocities:

      // Velocities:
      Array<OneD, Array<OneD, NekDouble>> u(dim);
      for (unsigned int i = 0; i < dim; ++i)
	{
	  u[i] = Array<OneD, NekDouble>(nPoints);
	}
      fluidEqu->GetVelocity(physfields, u);

      //density:
      Array<OneD,NekDouble> rho(nPoints,0.0);
      fluidEqu->GetDensity(physfields,rho);

      //pressure:
      Array<OneD,NekDouble> pressure(nPoints,0.0);
      fluidEqu->GetDensity(physfields,pressure);
      
      // TO-DO:

      // Just create a wrapper for the CompressibleFlowSystem
      // because there are methods that are not public...
       

 
  
      LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
      if(vComm->GetRank()==0)
	{
	  //Open a file at this time:
    
	}

    }//end update

    void FilterDissipation::v_Finalise(
				       const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
				       const NekDouble &time)
    {
      boost::ignore_unused(pFields, time);

      if (pFields[0]->GetComm()->GetRank() == 0)
	{
	  m_outputStream.close();
	}
    }

    bool FilterDissipation::v_IsTimeDependent()
    {
      return true;
    }

  } // namespace SolverUtils
} // namespace Nektar
