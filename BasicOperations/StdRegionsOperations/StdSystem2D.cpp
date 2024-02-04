#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <StdRegions/StdQuadExp.h>
#include <LocalRegions/QuadExp.h>

using namespace Nektar;
using namespace std;

int main(int, char **)
{

  /**
   * Various examples of the usage of 2d expansion (in the StdQuad region)
   */

  /// Create a 2D StdQuad expansion in the standard
  /// region:

  // Create a basiskey:
  unsigned int NUMMODES = 5;
  unsigned int NQx=10; unsigned int NQy =10; // Number of points of quadrature along the x direction
  
  LibUtilities::PointsKey pkeyDir1(NQx,LibUtilities::PointsType::eGaussLobattoLegendre);
  LibUtilities::BasisKey bDir1(LibUtilities::BasisType::eOrtho_A,NUMMODES,pkeyDir1);

  // We will use the same basis in both directions.
  StdRegions::StdQuadExp qExp(bDir1,bDir1);

  /// At this point, we can integrate out function directly with the Basis object.

  auto toIntegrate = [](Array<OneD,NekDouble>& f,const StdRegions::StdQuadExp& exp){
    /// Get the coordinates of the quadrature points:
    Array<OneD,NekDouble> xi1 = exp.GetPoints(0);
    Array<OneD,NekDouble> xi2 = exp.GetPoints(1);
  
    
    for(unsigned i =0;i<xi1.size();++i)
      for(unsigned j =0;j<xi2.size();++j)
	f[j+xi1.size()*i] = std::pow(xi1[i],1)*std::pow(xi2[j],1);

    return;
  };

  Array<OneD,NekDouble> f(NQx*NQy,0.0);
  toIntegrate(f,qExp);
  assert(f.size( ) == NQx*NQy);
  /// Integrate the function directly:
  const  auto basisArraySharedPtr = qExp.GetBase();

  // Integrate the functions over the 2D domain:
  NekDouble result = qExp.Integral(f,basisArraySharedPtr[0]->GetW(),basisArraySharedPtr[1]->GetW());

  std::cout << " The integral of (x1^6)*(x2^6) in the element region is " <<
    "the following: "<< result << "\n";

  /**
   * In the following example we consider the
   * projected problem onto the basis functions
   * (Discrete galerkin projection)
   */

  /// Step 1: Get the mass matrix:
  StdRegions::StdMatrixKey massKey(StdRegions::MatrixType::eMass,
				   LibUtilities::ShapeType::eQuadrilateral,qExp);

  StdRegions::StdMatrixKey InvmassKey(StdRegions::MatrixType::eInvMass,
				      LibUtilities::ShapeType::eQuadrilateral,qExp);
  DNekMatSharedPtr mass = qExp.GenMatrix(massKey);
  /// Invert:
  mass->Invert();

  const auto& InvMass = *mass;
  
  /// Step 2: We project the functions onto the basis:
  
  Array<OneD,NekDouble> fHat;
  qExp.IProductWRTBase_SumFac(f,fHat);

  // Get a nek-vector:
  NekVector<double> VecFhat(fHat);
  
  // Step3: We compute uHat:
  NekVector<double> VecuHat= Multiply(InvMass,VecFhat);
  std::cout << VecuHat << "\n";
  
  /// Perform the backward transformation
  /// and verify that is coherent with the nodal values
  /// in the quadrature points:
  Array<OneD,NekDouble> f_nodal(f.size(),0.0);
  qExp.BwdTrans(VecuHat.GetPtr(),f_nodal);

  std::cout << "Result of the Discrete Galerkin projection: \n";
  std::cout << "Qx,Qy,NmodesX,NmodesY: "<< NQx <<" " << NQy << " " << NUMMODES << "\n";
  std::cout << "Error L^2 norm: "<< qExp.L2(f_nodal,f) << "\n";

 
  /*  
  for(unsigned int i=0;i<f_nodal.size();++i)
    std::cout <<f[i] << " " << f_nodal[i]<<"\n";
  */

  /// Moving from the StandardRegion to the LocalRegion:

  // Define a geometry2D:

  using namespace SpatialDomains;
  // Generate a points geom:
  PointGeom p1(2,0,-0.5,-0.5,0.0);
  PointGeom p2(2,1,+0.5,0.0,0);
  PointGeom p3(2,2,-0.2,0.2,0);
  PointGeom p4(2,3,-0.4,0.1,0);


  PointGeomSharedPtr e1[2];
  PointGeomSharedPtr e2[2];
  PointGeomSharedPtr e3[2];
  PointGeomSharedPtr e4[2];
  e1[0] = std::make_shared<PointGeom>(p1);
  e1[1] = std::make_shared<PointGeom>(p2);
  e2[0] = std::make_shared<PointGeom>(p2);
  e2[1] = std::make_shared<PointGeom>(p3);
  e3[0] = std::make_shared<PointGeom>(p3);
  e3[1] = std::make_shared<PointGeom>(p4);
  e4[0] = std::make_shared<PointGeom>(p4);
  e4[1] = std::make_shared<PointGeom>(p1);

  // Generate 4 different edges:
  SpatialDomains::SegGeom s1(0,2,e1);
  SpatialDomains::SegGeom s2(1,2,e2);
  SpatialDomains::SegGeom s3(2,2,e3);
  SpatialDomains::SegGeom s4(3,2,e4);
  
  SpatialDomains::SegGeomSharedPtr edges[4];

  edges[0] = std::make_shared<SpatialDomains::SegGeom>(s1);
  edges[1] = std::make_shared<SpatialDomains::SegGeom>(s2);
  edges[2] = std::make_shared<SpatialDomains::SegGeom>(s3);
  edges[3] = std::make_shared<SpatialDomains::SegGeom>(s4);
  
  SpatialDomains::QuadGeom my_geometry(0,edges);

  SpatialDomains::QuadGeomSharedPtr pGeom = std::make_shared<SpatialDomains::QuadGeom>(my_geometry);
  
  LocalRegions::QuadExp locQuadExp(bDir1,bDir1,pGeom);

  // Here, we have the expansion in the local region and
  // we can use functions to integrate directly in
  // this element.
  
  
  
  

  

  

  
  
  
  

}
