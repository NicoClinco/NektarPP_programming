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

  /// Remember that the session reader will read <mesh> <session.xml>
  /// as arguments:
  LibUtilities::SessionReaderSharedPtr pcurrentSession;
  pcurrentSession = LibUtilities::SessionReader::CreateInstance(argc,argv);

  /// Create the mesh graph from Session-Reader:
  SpatialDomains::MeshGraphSharedPtr graph = SpatialDomains::MeshGraph::Read(pcurrentSession);
  

  /// Read the geometry given the filename:
  graph->ReadExpansionInfo();
  
  /// Mesh-dimension:
  std::cout << "Mesh dimension :" << graph->GetMeshDimension() << "\n";
  std::cout << "Space dimension :" << graph->GetSpaceDimension() << "\n";

  /// Get the total number of elements (Should be 216) 
  std::cout << "Total number of elements: " << graph->GetNumElements() <<"\n";
  
  /// Create the geometry which contains vertexes, elements id and so on.

  /// Here i return a "geometry" object that identify the zero element
  /// of the composite number = 0. GetCompositeItem(composite,elemofcomposite)
  const auto& geometry = graph->GetCompositeItem(0,0);

  /// For what i have understood the geometry is the base class object
  /// that has data-members related to the vertex-coordinates, vertex
  /// ids and so on.

  ///Plot some information about number of vertexes, ids and other related topics.

  std::cout<< "The element "<< geometry->GetGlobalID() << " has:\n"
	   << geometry->GetNumFaces() << " faces\n"
	   << geometry->GetNumVerts() << " vertexes\n";

  ///Loop trough the element faces:
  std::cout << "Global faces ids \n";
  for(unsigned int i=0;i< geometry->GetNumFaces();++i)
    std::cout << geometry->GetFid(i) << " ";
  std::cout << "\n";
  ///Get the LOCAL COORDINATES: in the element:
  Array<OneD, NekDouble> glo_coord(3,-3.14);
  Array<OneD, NekDouble> loc_coord(3,0.0);

  /// If the element contains the global
  /// coordinate, we are ok :D
  if(geometry->ContainsPoint(glo_coord))
    {
      geometry->GetLocCoords(glo_coord,loc_coord);
      for(const auto& lc : loc_coord)
	std::cout << lc << " ";
    }
  ///At the end of the story, the geometry object is an object which
  /// contains a bunch of vertexes, faces id.

  /*
  std::cout << "\n";
  Array<OneD,NekDouble> a(3,2.0);
  std::vector<Array<OneD,NekDouble>> b(2,Array<OneD,NekDouble>(3,1.0));
  /// Try to copy the array:
  b[0] = a;
  for(const auto& v : b)
    for(const auto& e : v)
      std::cout << e << " ";
  */
  
  
  
  return 0;
}
