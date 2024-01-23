#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/BasicUtils/FieldIOXml.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
  std::cout << "Reading fields from xml files"<<std::endl;

  // Creating the Session-reader pointer
  // used for taking a communication object:
  LibUtilities::SessionReaderSharedPtr pcurrentSession;
  pcurrentSession =  LibUtilities::SessionReader::CreateInstance(argc,argv);
  LibUtilities::CommSharedPtr pComm = pcurrentSession->GetComm();
  
  
  // Creating a shared pointer to a FieldIOXml
  LibUtilities::FieldIOSharedPtr pFieldIOxml;


  // Get the pointer for the common comunicator:
  pFieldIOxml = LibUtilities::FieldIOXml::create(pComm,false);


  //Import an Xml file:
  std::string Inputfile("TGV6p3_ALIGNED_1.chk");
  std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDefs;
  std::vector<std::vector<NekDouble>> FieldData;
  LibUtilities::FieldMetaDataMap FieldInfoMap;

  pFieldIOxml->Import(Inputfile,FieldDefs,FieldData,FieldInfoMap);

  // Perform a cycle for every fld files
  // that is present in the domain:
  for(auto pfd : FieldDefs)
    {
      auto fd = *pfd;
      
      std::vector<std::string> fields = fd.m_fields;
      for(auto field : fields)
	std::cout << field << " ";
      std::cout << "\n";
    }
  
  
  /*
  for(auto row :FieldData)
    {
      
    for(auto col : row)
      std::cout << col << " ";
    std::cout << "\n";
    }
  */
  
  
  
  
  
  return 0;
}
