#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <vector>


#include <FieldUtils/Field.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>

using namespace Nektar;

/*
 Beside (or on the top of FieldIO) there is the Field struct.

 This class is the "top class" for the Field-utility module
 that allows us to post processing with the FieldConvert 
 module.
 */


int main(int argc, char *argv[])
{
  std::cout << "### EXAMPLES ON THE FIELD-UTILS CLASS ####";

  //Get the session-reader to get the mesh and the session-file
  LibUtilities::SessionReaderSharedPtr pcurrentSession;
  pcurrentSession =  LibUtilities::SessionReader::CreateInstance(argc,argv);
  LibUtilities::CommSharedPtr pComm = pcurrentSession->GetComm();
  
  
  
  //Basic constructor:
  auto FieldObj = FieldUtils::Field();

  //Since the FieldObj is a struct, its members are public by default,
  //thus we can write internal member simply by accessing:

  //Copy the shared pointer of the session-file (copy)
  FieldObj.m_session = pcurrentSession;

  //Copy the shared pointer of the mesh-graph:
  FieldObj.m_graph = SpatialDomains::MeshGraph::Read(pcurrentSession);

  //Create another ExpList:

  //Note: I suppose that "Map" will be given from the Input module
  //in the future
  std::map<std::string,bool> mapStandard = {{"useSessionExpansion",true}};
  auto new_exp = FieldObj.CreateExp(mapStandard);
  

  return 0;
}
