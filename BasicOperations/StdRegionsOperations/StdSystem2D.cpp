#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <StdRegions/StdQuadExp.h>

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
  unsigned int NUMMODES = 3;
  unsigned int NQx=5; unsigned int NQy =5; // Number of points of quadrature along the x direction
  
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
	f[j+xi1.size()*i] = std::pow(xi1[i],6)*std::pow(xi2[j],6);

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
  DNekMatSharedPtr mass = qExp.GetStdMatrix(massKey);
  /// Invert:
  mass->Invert();
  
  /// Step 2: We project the functions onto the basis:
  
  Array<OneD,NekDouble> fHat;
  qExp.IProductWRTBase(f,fHat);

  
  // Get a nek-vector:
  NekVector<double> VecFhat(fHat);
  
  // Step3: We compute uHat:
  NekVector<double> VecuHat= Multiply(*mass,VecFhat);

  /// To do: Perform the backward transformation
  /// and verify that is coherent with the nodal values.
  
  std::cout << VecuHat ;

  

  
  
  
  

}
