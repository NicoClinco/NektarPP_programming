#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <vector>

#include <LibUtilities/Foundations/ManagerAccess.h>
//#include <LibUtilities/Polylib/Polylib.h>
//#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
//#include <LibUtilities/BasicUtils/FieldIOXml.h>
#include <SpatialDomains/MeshGraph.h>
/// Read a mesh and put everything i a mesh graph

using namespace Nektar;

int main(int argc, char *argv[])
{
  /// Pointer to a session reader:  
  LibUtilities::SessionReaderSharedPtr pcurrentSession;
  pcurrentSession = LibUtilities::SessionReader::CreateInstance(argc,argv);

  /// Create the mesh graph
  SpatialDomains::MeshGraphSharedPtr graph = SpatialDomains::MeshGraph::Read(pcurrentSession);
  

  /// Read the geometry given the filename:
  graph->ReadExpansionInfo();
  
  /// Mesh-dimension:
  std::cout << "Mesh dimension :" << graph->GetMeshDimension() << "\n";
  std::cout << "Space dimension :" << graph->GetSpaceDimension() << "\n";
  
  
  

  
  
  return 0;
}
