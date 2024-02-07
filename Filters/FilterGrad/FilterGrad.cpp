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

    int i;
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    const Array<OneD,NekDouble> field0 = pFields[0]->GetPhys();

    ///Create a vector that contains the the array ux,uy,uz:
    std::vector<Array<OneD,NekDouble>> V_xyz(3,Array<OneD,NekDouble>(field0.size(),0.0));

    /// Boolean-vector to check if something is wrong.
    std::vector<bool> found{false,false,false};
    
    // Get the physical fields:
    for(unsigned int i = 0;i<pFields.size();++i)
      {
	/// Get the field at the quadrature point:
	const Array<OneD,NekDouble> field = pFields[i]->GetPhys();

	/// Get the field definition:
	std::vector<LibUtilities::FieldDefinitionsSharedPtr> fieldDefs = pFields[i]->GetFieldDefinitions();
	
	std::string nameField = fieldDefs[0]->m_fields[0];
	if(nameField == "u")
	  {
	    Array<OneD,NekDouble> ux(field.size(),0.0);
	    pFields[i]->PhysDeriv(0,field,ux);
	    ///Copy:
	    V_xyz[0] = ux;
	    found[0] = true;
	  }
	if(nameField == "v")
	  {
	    Array<OneD,NekDouble> uy(field.size(),0.0);
	    pFields[i]->PhysDeriv(1,field,uy);
	    V_xyz[1] = uy;
	    found[1] = true;
	  }
	if(nameField == "w")
	  {
	    Array<OneD,NekDouble> uz(field.size(),0.0);
	    pFields[i]->PhysDeriv(2,field,uz);
	    V_xyz[2] = uz;
	    found[2] = true;
	  }
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
	
      }
    
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
