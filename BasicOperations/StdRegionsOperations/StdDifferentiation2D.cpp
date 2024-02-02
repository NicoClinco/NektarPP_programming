#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>

/*
 Tutorial for differentiation in Nektar
 StdDifferentiation : Refered to the standard region

 TO DO :
 Try to differentiation in the local region
*/


using namespace Nektar;
using namespace std;

int main(int, char **)
{
    cout << "===========================================================" << endl;
    cout << "|    DIFFERENTIATION IN 2D ELEMENT in Standard Region     |" << endl;
    cout << "===========================================================" << endl;
    cout << endl;
    cout << "Differentiate the function f(x1,x2) = (x1)^7*(x2)^9" << endl;
    cout << "in the standard quadrilateral element:" << endl;

    // Specify the number of quadrature points in both directions
    int nQuadPointsDir1 = 4;
    int nQuadPointsDir2 = 4;

    // Specify the type of quadrature points in both directions
    LibUtilities::PointsType quadPointsTypeDir1 =
        LibUtilities::eGaussLobattoLegendre;
    LibUtilities::PointsType quadPointsTypeDir2 =
        LibUtilities::eGaussLobattoLegendre;

    // Declare variables (of type Array) to hold the quadrature zeros
    // and the values of the derivative at the quadrature points in both
    // directions
    Array<OneD, NekDouble> quadZerosDir1(nQuadPointsDir1);
    Array<OneD, NekDouble> quadZerosDir2(nQuadPointsDir2);
    Array<TwoD, NekDouble> quadDerivsDir1(nQuadPointsDir1, nQuadPointsDir2);
    Array<TwoD, NekDouble> quadDerivsDir2(nQuadPointsDir1, nQuadPointsDir2);

    // Declare pointers (to type NekMatrix<NekDouble>) to hold the
    // differentiation matrices
    DNekMatSharedPtr derivMatrixDir1;
    DNekMatSharedPtr derivMatrixDir2;

    // Calculate the GLL-quadrature zeros and the differentiation
    // matrices in both directions. This is done in 2 steps.

    // Step 1: Declare the PointsKeys which uniquely defines the
    // quadrature points
    const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1,
                                                    quadPointsTypeDir1);
    const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2,
                                                    quadPointsTypeDir2);

    // Step 2: Using this key, the quadrature zeros and differentiation
    // matrices can now be retrieved through the PointsManager
    quadZerosDir1   = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetZ();
    derivMatrixDir1 = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetD();

    quadZerosDir2   = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetZ();
    derivMatrixDir2 = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetD();


    /// Create the stiffness-matrix <gradPhii,gradPhij>

    //Number of modes in the elements:
    unsigned int NUMMODES = 2;

    // Create a basiskey by the Gauss-Legendre basis:
    const LibUtilities::BasisKey basisKey(LibUtilities::BasisType::eOrtho_A,
					  NUMMODES,quadPointsKeyDir1);

    auto pBasisDir1 = LibUtilities::Basis::Create(basisKey);
    // Create basis for the other directions:
    auto pBasisDir2 = LibUtilities::Basis::Create(basisKey);
    
    // Create NekMatrixes  which contain in every row the derivatives
    // of the basis in the quadrature points:
    NekMatrix<double> dPhix(nQuadPointsDir1,NUMMODES,pBasisDir1->GetDbdata());
    NekMatrix<double> Phix(nQuadPointsDir1,NUMMODES,pBasisDir1->GetBdata());
    NekVector<double> Wx(pBasisDir1->GetW()); //Weights along the x direction
    
    NekMatrix<double> dPhiy(nQuadPointsDir2,NUMMODES,pBasisDir2->GetDbdata());
    NekMatrix<double> Phiy(nQuadPointsDir2,NUMMODES,pBasisDir2->GetBdata());
    NekVector<double> Wy(pBasisDir2->GetW()); // Weights along the y-direction

    ASSERTL0(Wx.GetDimension()==nQuadPointsDir1 ,"WARNING, THE DIMENSIONS OF THE WEIGHTS ARE != FROM THE NODES dir1");
    ASSERTL0(Wy.GetDimension()==nQuadPointsDir2 ,"WARNING, THE DIMENSIONS OF THE WEIGHTS ARE != FROM THE NODES dir2");

    /// Note: for linear-algebra operations, we need a NekVector
    /// that is conceived for linear-algebra operations.

    /**
     * @ brief Perform the scalar product
     *
     */
    auto IntegrateInNodes = [](const NekMatrix<double>& phi1,const NekMatrix<double>& phi2,
			       unsigned int i,unsigned int j,const NekVector<double>& W)
    {
      ASSERTL0(phi1.GetRows()==phi2.GetRows(),"The matrixes must have the same number of quadrature points!!");
      double Value = 0.0;
      for(unsigned int nq =0;nq<phi1.GetRows();++nq)
	Value+=phi1(nq,i)*phi2(nq,j)*W[nq];
      return Value;
    };
    
    // This time we use a NekMatrix<double>:
    NekMatrix<double> K(NUMMODES*NUMMODES,NUMMODES*NUMMODES);
    
    for(unsigned int i=0;i<NUMMODES;++i)
      {
	for(unsigned int j=0;j<NUMMODES;++j)
	  {
	    unsigned int a = i+j*NUMMODES;
	    for(unsigned int r=0;r<NUMMODES;++r)
	      {
		for(unsigned int k=0;k<NUMMODES;++k)
		  {
		    unsigned int b = r+k*NUMMODES;
		    auto grgrx = IntegrateInNodes(dPhix,dPhix,i,r,Wx);
		    auto grgry = IntegrateInNodes(dPhiy,dPhiy,j,k,Wy);
		    auto phixphix = IntegrateInNodes(Phix,Phix,i,r,Wx);
		    auto phiyphiy = IntegrateInNodes(Phiy,Phiy,j,k,Wy);
		    K(a,b) = grgrx*phiyphiy+grgry*phixphix;
		  }
	      }
	  }
      }

    // Iterate through the matrix using iterators
    std::cout << K << "\n";

    /// We can create a diagonal matrix:
    /// See blockMatrix.hpp for further references:
    NekMatrix<double> K_diag(K.GetRows(),K.GetColumns(),MatrixStorage::eDIAGONAL);
    for(unsigned int i =0;i<K_diag.GetRows();++i)
      K_diag(i,i) = K(i,i);

    std::cout << K_diag << "\n";
    
    

    
    
    
    // Now you have the quadrature zeros and the differentiation matrix,
    // apply the Gaussian quadrature technique to differentiate the function
    // f(x_1,i,x_2,j) = x_1,i^7 * x2,j^9 on the standard
    // quadrilateral.  To do so, write a (double) loop which performs
    // the summation.
    //
    // Store the solution in the matrices 'quadDerivsDir1' and 'quadDerivsDir2'

    /*
    // Analytic function:
    auto f = [] (NekDouble x, NekDouble y)
    {
      return pow(x,7)*pow(y,9);
    };

    // Dir1
    for(size_t j=0;j<nQuadPointsDir2;++j){
      
      for(size_t i =0;i<nQuadPointsDir1;++i)
	{
	  for(size_t k =0;k<nQuadPointsDir1;++k)
	    {
	      quadDerivsDir1[i][j]+=(*derivMatrixDir1)(i,k) *
		f(quadZerosDir1[k],quadZerosDir2[j]);
	    }
	}
    }
    
     for(size_t i=0;i<nQuadPointsDir1;++i){
      
      for(size_t j=0;j<nQuadPointsDir2;++j)
	{
	  for(size_t k=0;k<nQuadPointsDir2;++k)
	    {
	      quadDerivsDir2[i][j]+= (*derivMatrixDir2)(j,k) *
		f(quadZerosDir1[i],quadZerosDir2[k]);
	    }
	}
    }

    // Compute the total error
    NekDouble error = 0.0;
    for (size_t i = 0; i < nQuadPointsDir1; ++i)
    {
        for (size_t j = 0; j < nQuadPointsDir2; ++j)
        {
            error +=
                fabs(quadDerivsDir1[i][j] -
                     7 * pow(quadZerosDir1[i], 6) * pow(quadZerosDir2[j], 9));
            error +=
                fabs(quadDerivsDir2[i][j] -
                     9 * pow(quadZerosDir1[i], 7) * pow(quadZerosDir2[j], 8));
        }
    }

    // Display the output
    cout << "\t q1 = " << nQuadPointsDir1 << ", q2 = " << nQuadPointsDir2;
    cout << ": Error = " << error << endl;
    cout << endl;
    */
}
