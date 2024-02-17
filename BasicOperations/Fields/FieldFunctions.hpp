/**
 * This file contains various functions useful for reading and operating
 * on fields (post-processing utilities)
 *
 */

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/BasicUtils/FieldIOXml.h>
#include <FieldUtils/Field.hpp>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

using namespace Nektar;

namespace MathLab
{
  namespace Utilities
  {
/**
 * @brief Extract the physical field
 * specified by the user from a file
 *
 */
Array<OneD,double>
ExtractPhysFromFile
(
 const LibUtilities::SessionReaderSharedPtr& pSession,
 const MultiRegions::ExpListSharedPtr& pExp,
 const std::string& fieldname,
 const std::string& filename
 )
{
  try{
    //
    Array<OneD,double> coeffs(pExp->GetNcoeffs());
    pExp->ExtractCoeffsFromFile(filename,pSession->GetComm(),fieldname,coeffs);
    return coeffs;
  }
  catch(...)
    {
      std::cout << "It was not possible to extract "<< fieldname <<" from ";
      std::cout << filename << "\n";
      return Array<OneD,double>();
    } 
}


using StringList = std::vector<std::string>;


/**
 * @brief Extract the physical fields
 * specified by the user from a file
 *
 * 
 */
std::vector<Array<OneD,double>>
ExtractPhysFromFile
(
 const LibUtilities::SessionReaderSharedPtr& pSession,
 const MultiRegions::ExpListSharedPtr& pExp,
 const StringList& fieldname,
 const std::string& filename
 )
{
  std::vector<Array<OneD,double>> fields(fieldname.size(),Array<OneD,double>(pExp->GetNcoeffs(),0.0));
  unsigned int i=0;
  for(const auto& f : fieldname)
    {
      try
	{
	  Array<OneD,double> current_array(pExp->GetNcoeffs());
	  pExp->ExtractCoeffsFromFile(filename,pSession->GetComm(),fieldname[i],current_array);
	  fields[i] = current_array;
	}
      catch(...)
	{
	  std::cout << "It was not possible to extract the physical coefficients for ";
	  std::cout << f << " field\n";
	  std::cout << "Returning an empty vector instead "<<"\n";
	  return
	    std::vector<Array<OneD,double>>(1,Array<OneD,double>());
	}
      ++i;
    }//end cycle

  //Return a copy of the physical fields:
  return fields;
}//end function

using ArrayList =  std::vector<Array<OneD,double>>;

/**
 * @brief Compute physical gradients in a specific
 *        direction set by "dir"
 * 
 */
void ComputePhysGradients
(
 const MultiRegions::ExpListSharedPtr& pExp,
 const ArrayList& in_fields,
 ArrayList& out_fields,
 unsigned int dir
 )
{
  if(out_fields.size()!=in_fields.size())
    out_fields.resize(in_fields.size());
  
  for(unsigned int i=0;i<in_fields.size();++i)
    {
      
      //Array<OneD,NekDouble> in_array = in_fields[i];
      //Array<OneD,NekDouble> out_array(in_array.size());
      //std::cout << "size :"<< out_array.size() << " "<< in_array.size() <<  "\n"; 
      pExp->PhysDeriv(MultiRegions::Direction::eX,in_fields[i],out_fields[i]);
    } 
  return;
}


  }//end namespace utilities

}//end namespace mathlab


