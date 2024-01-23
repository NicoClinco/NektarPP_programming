#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdMatrixKey.h>

#include <Eigen/Dense>

using namespace Nektar;
using namespace std;

int main(int, char **)
{
    cout << "======================================================" << endl;
    cout << "|           Quadrature in a Standard Region          |" << endl;
    cout << "======================================================" << endl;

    cout << "Integrate the function f(xi) = x1^5*x2^4 in the " << endl;
    cout << "standard segment xi=[-1,1]x[-1,1] using quadrature points" << endl;

    // Specify the number of quadrature points
    int nQuadPoints = 7;

    // Specify the type of quadrature points. This is done using the proper
    // Nektar++ syntax.
    //LibUtilities::PointsType quadPointsType = LibUtilities::eGaussGaussLegendre;
    LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;
    Array<OneD, NekDouble> quadZeros(nQuadPoints);
    Array<OneD, NekDouble> quadWeights(nQuadPoints);

    // Step 1: Declare a PointsKey which uniquely defines the
    // quadrature points
    const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);

    // Step 2: Using this key, the quadrature zeros and the differentiation
    // matrix can now be retrieved through the PointsManager in namespace
    // LibUtilities
    quadZeros   = (LibUtilities::PointsManager()[quadPointsKey])->GetZ();
    quadWeights = (LibUtilities::PointsManager()[quadPointsKey])->GetW();
    

    auto fx1 = [](const double& xi){ return pow(xi,2);};
    auto fx2 = [](const double& xi2){ return pow(xi2,2);};

    double I1,I2 = 0.0;
    I1 = I2 = 0;
    //Evaluate the integrals
    for(unsigned int i=0;i<quadZeros.size();++i)
      {
	std::cout << quadZeros[i] << " ";
	I1+=fx1(quadZeros[i])*quadWeights[i];
      }
    std::cout <<"\n";
    for(unsigned int i=0;i<quadZeros.size();++i)
      {
	I2+=fx2(quadZeros[i])*quadWeights[i];
      }
    std::cout << "Value of the integrals :" << I1*I2 << "\n";
    
    /**
     * @brief Anothere way to to the things, 
     * thanks to the Basis class
     */
    const unsigned int Nmodes = 3;
    const LibUtilities::PointsKey basisPointsKey(nQuadPoints, quadPointsType);
    const LibUtilities::BasisKey basisKey(LibUtilities::BasisType::eOrtho_A,
					  Nmodes,basisPointsKey);
    
    auto pBasis = LibUtilities::Basis::Create(basisKey);

    // Check if the basis is evaluated at the specified set of points:
    Array<OneD, NekDouble> basisValues = pBasis->GetBdata();

    /**
     *
     * Note: the basis are generated in a matrix with the following format:
     * [  | ------ Npoints ------  ]
     * [Nmodes                     ]
     * [  |                        ]
     */

    auto basisVal = [basisValues,nQuadPoints]
      (const unsigned int& r, const unsigned int& c)
    {
      return basisValues[r*nQuadPoints+c];
    };

    // The same is true for r,s:
    auto n = [Nmodes](const unsigned int&  p,const unsigned int&  q)
    {
      return p*Nmodes+q;
    };
    
    // Create the mass matrix and store it into eigen:
    const unsigned int P = Nmodes-1;
    
    Eigen::MatrixXd mass=Eigen::MatrixXd::Zero(Nmodes*Nmodes,Nmodes*Nmodes);
    
    for(unsigned int p=0;p<Nmodes;++p)
      {
	for(unsigned int q=0;q<Nmodes;++q)
	  {
	    for(unsigned int r=0;r<Nmodes;++r)
	      {
		for(unsigned int s=0;s<Nmodes;++s)
		  {
		    // Evaluate the first formula (p,q)
		    double Sum1 = 0.0;
		    for(unsigned int i=0;i<nQuadPoints;++i)
		      {
			Sum1+= quadWeights[i]*basisVal(p,i)*basisVal(r,i);
		      }

		    //Evaluate the second formula (r,s)
		    double Sum2 = 0.0;
		    for(unsigned int i=0;i<nQuadPoints;++i)
		      {
			Sum2+= quadWeights[i]*basisVal(q,i)*basisVal(s,i);
		      }
		    mass(n(p,q),n(r,s)) = Sum1*Sum2;
		   
		  }
	      }
	  }// End complete loop with q
      }// End complete loop with p

    // In this case we should obtain a diagonal
    // matrix :D
   
    //std::cout << mass << "\n";

    /**
     *@brief In the following example, 
     *we generate an StdQuadExpansion
     */
    StdRegions::StdQuadExp Quad_exp(basisKey,basisKey);
    //StdRegions::DNekMatSharedPtr matrix(StdRegions::StdMatrixKey::e

    //For obtaining the mass matrix, we must set the stdMatrixKey
    StdRegions::StdMatrixKey mkey(StdRegions::MatrixType::eMass,
			     LibUtilities::ShapeType::eQuadrilateral,
			     Quad_exp);
    
    DNekMatSharedPtr pmatrix = Quad_exp.GetStdMatrix(mkey);
    const auto& matrix = *pmatrix;

    //In this case we can obtain the Physical derivatives:
    //at the quadrature points:

    Array<OneD,NekDouble> Derivatives(2);
    Quad_exp.PhysDeriv(1,quadZeros,Derivatives);

    for(const auto& x : Derivatives)
      std::cout << x << " ";
   
    /*
    for(unsigned int r=0;r<matrix.GetRows();++r)
      {
	for(unsigned int c=0;c<matrix.GetColumns();++c)
	  {
	    std::cout << matrix(r,c) << " ";
	  }
	std::cout << "\n";
      }
    */
    
    //std::cout << "Number of rows: " << matrix.GetRows() << "\n";
    //std::cout << "Number of columns: "<< matrix.GetColumns() << "\n";
    
    
    
    
    
    
    
}
