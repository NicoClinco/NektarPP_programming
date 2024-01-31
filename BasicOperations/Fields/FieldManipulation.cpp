#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <vector>

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


  //Import a .chk file:
  std::string Inputfile("TGV6p3_ALIGNED_1.chk");
  std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDefs;
  std::vector<std::vector<NekDouble>> FieldData;
  LibUtilities::FieldMetaDataMap FieldInfoMap;

  pFieldIOxml->Import(Inputfile,FieldDefs,FieldData,FieldInfoMap);


  // Read the field Meta data:
  /*
   for (const auto& pair : FieldInfoMap) {
        std::cout << " " << pair.first << " " << pair.second << "\n";
    }
  */

  /**
   * FieldData will contain fields defined in FieldDefs.
   */
  
  /// Perform a cycle for every fld file
  /// that is contained in the directory:
  for(auto pfd : FieldDefs)
    {
      auto fd = *pfd;
      ///Every single file contains different fields:
      std::vector<std::string> field_names = fd.m_fields;

      /// Field names:
      ///for(auto fieldname : field_names)
      ///	std::cout << fieldname << " ";

      /// Number of modes per direction:
      ///std::vector<unsigned int> nummodes = fd.m_numModes;
      
      /*
      std::cout << "Number of points x-y-z:" << fd.m_numPoints[0] << " "
       	<< fd.m_numPoints[1] << " "<< fd.m_numPoints[2] <<"\n";
      */
      }

  /// TO DO: Read the mesh file and introduce an ordering for
  /// understanding the coordinates.
  
  /// Printing the fields as vectors: size = 20:
  std::cout << "Field-data-size :" <<  FieldData.size() <<"\n";

  /// Note: Every field is stored in a contiguos vector
  /// for each fld file.
  for(const auto& phi : FieldData)
    std::cout << phi.size() << " ";
  std::cout << "\n";

  
  // Creating a specific class for read multifld files:
  LibUtilities::FieldIOXml FieldIOxml(pComm,false);

  

  
  //In this case we read directly the .fld files and we get the
  //elements partition:
  LibUtilities::FieldMetaDataMap FieldInfoMapPar;

  //ElementsID
  std::vector< std::vector< unsigned int >> elementsID;
  //.fld files present in the directory
  std::vector<std::string> fileNames;
  
  std::string InfoFile("./TGV6p3_ALIGNED_1.chk/Info.xml");
  FieldIOxml.ImportMultiFldFileIDs(InfoFile, fileNames, elementsID,FieldInfoMapPar);

  //Print the .fld files:
  for(const auto& file : fileNames)
    std::cout << file << " ";

  //Print the elements id:
  //Here we have the elements id.
  unsigned int total_elems = 0;  
  for(const auto& els : elementsID)
    total_elems+=els.size();
  

  std::cout << total_elems << "\n";

  
  
  
  
  return 0;
}
