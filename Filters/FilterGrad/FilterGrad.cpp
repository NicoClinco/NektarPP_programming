///////////////////////////////////////////////////////////////////////////////
//
// File: FilterGrad.cpp
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
// Description: Output mean.
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <SolverUtils/Filters/FilterGrad.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
std::string FilterGrad::className =
    SolverUtils::GetFilterFactory().RegisterCreatorFunction("Grad",
                                                            FilterGrad::create);

FilterGrad::FilterGrad(const LibUtilities::SessionReaderSharedPtr &pSession,
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

    pSession->LoadParameter("LZ", m_homogeneousLength, 0.0);
}

FilterGrad::~FilterGrad()
{
}

void FilterGrad::v_Initialise(
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
      m_outputStream << "# Time,ux,uy,uz" << "\n";
      /*
      for (int i = 0; i < pFields.size(); ++i){
	m_outputStream << setw(22) << equ->GetVariable(i);
      }
      */
    }

    // Output values at initial time.
    m_index = 0;
    v_Update(pFields, time);
}

void FilterGrad::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
  // Only output every m_outputFrequency
  if ((m_index++) % m_outputFrequency)
    {
      return;
    }
  
  auto equ = m_equ.lock();
  ASSERTL0(equ, "Weak pointer expired");

  // Get the fluid-interface as object for obtain the velocities and pressure:
  auto fluidEqu = std::dynamic_pointer_cast<FluidInterface>(equ);
  ASSERTL0(fluidEqu, "Grad-filter is incompatible with this solver.");


  // Get the physical-fields:
  Array<OneD, Array<OneD, NekDouble>> physfields(pFields.size());
  for (unsigned int i = 0; i < pFields.size(); ++i)
    {
      physfields[i] = pFields[i]->GetPhys();
    }
  
  unsigned int dim = 0;
  unsigned int nPoints = pFields[0]->GetNpoints();
  // Get the expansion type:
  if (pFields[0]->GetExpType() == MultiRegions::e1D)
    {
      // 1D
      dim = 1;      
    }
  if (pFields[0]->GetExpType() == MultiRegions::e2D)
    {
      // 2D
      dim=2;
    }
  if (pFields[0]->GetExpType() == MultiRegions::e3D)
    {
      // 3D
      dim=3;
    }

  // Velocities:
  Array<OneD, Array<OneD, NekDouble>> u(dim);
  for (unsigned int i = 0; i < dim; ++i)
    {
      u[i] = Array<OneD, NekDouble>(nPoints);
    }
  fluidEqu->GetVelocity(physfields, u);

  using doubleArray = Array<OneD,NekDouble>;
  using vecArray = std::vector<doubleArray>;
  //Compute the gradients:
  // ux,uy,uz,vx,vy,vz,wx,wy,wz
  std::vector<Array<OneD,NekDouble>> gUcmpt(dim,doubleArray(nPoints,0.0));
  //std::vector<Array<OneD,NekDouble>> gradV;
  //std::vector<Array<OneD,NekDouble>> gradW;

  // vector which contains: [ [ux,uy,uz],[vx,vy,vz],[wx,wy,wz] ]
  //initialized:
  std::vector<vecArray> gradU(dim,gUcmpt);


  //Store the derivatives (gradients):
  for(unsigned int uCMPT =0;uCMPT<dim;++uCMPT)
    {
      //derivative of u:
      for(unsigned int gu=0;gu<dim;++gu)
	{
	  pFields[0]->PhysDeriv(gu,u[uCMPT],gradU[uCMPT][gu]);
	}
    }
  
  LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
  
  
  //Open a file at this time:
  ofstream current_file;
  std::string current_filename = "TGV_uvw_"+std::to_string(time)+".csv";
  m_outputStream << time << " " << current_filename << "\n";
  //Open-File:
  // ux[pi] uy[pi]
  current_file.open(current_filename,std::ios_base::app);

  for(unsigned int p=0;p<nPoints;++p)
    {
      for(unsigned int gu=0;gu<dim;++gu)
	{
	  for(unsigned int ucmpt=0;ucmpt<dim;++ucmpt)
	    {
	      current_file << gradU[ucmpt][gu][p] <<","; 
	    }
	}
      //END_LINE
      current_file << "\n"; 
    }
   
      
  current_file.close();
   
	
  /*
    if(std::all_of(found.begin(),found.end(),[](bool f){ return f;}))
    {
    /// Write to a certain file:
    if(vComm->GetRank() == 0)
    {
		
    ///Open a file and save at this time what happen:
    ofstream current_file;
    std::string current_filename = "TGV_uvw_"+std::to_string(time) + ".csv";
		
    m_outputStream << time << " " << current_filename << "\n";
		
    current_file.open(current_filename);

    // Check if the file is good: TO DO
		
    for(unsigned int i = 0;i<3;++i)
    {
    for(unsigned int j = 0; j<V_xyz[0].size();++j)
    {
    /// Save ux,uy,uz
    current_file << V_xyz[0][j] << "," << V_xyz[1][j] << "," << V_xyz[2][j] << "\n";
    }
    }
    current_file.close();
    }
	    
    }
  */
}

void FilterGrad::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(pFields, time);

    if (pFields[0]->GetComm()->GetRank() == 0)
    {
        m_outputStream.close();
    }
}

bool FilterGrad::v_IsTimeDependent()
{
    return true;
}

} // namespace SolverUtils
} // namespace Nektar
