#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>

using namespace Nektar;


int main(int argc, char *argv[])
{
  std::cout << "Application to read the SessionReader file as an option\n"<<std::endl;

  // Create a SessionReader file (thanks to the factories)
  LibUtilities::SessionReaderSharedPtr pcurrentSession;

  try
    {
      std::vector<std::string> vecFiles{"session.xml"};
      pcurrentSession =  LibUtilities::SessionReader::CreateInstance(argc,argv,vecFiles);
      //Initialize the session reader:
      pcurrentSession->InitSession();

      LibUtilities::SessionReader& currentSession = *pcurrentSession;
      // Read parameter from file: it will print in the command line:
      if(currentSession.DefinesParameter("num0"))
	{
	  NekDouble num0 = currentSession.GetParameter("num0");
	  std::cout << "Parameter num0 specified :\n";
	  std::cout << num0 <<" ";
	  //Alternative:
	  NekDouble num0_;
	  currentSession.LoadParameter("num0",num0_);
	  std::cout << num0_ <<"\n";
	}

      // Get the solver information:
      if(currentSession.DefinesSolverInfo("EquationOfState"))
	{
	  std::cout <<"EOS : " <<  currentSession.GetSolverInfo("EquationOfState") << "\n";
	}

      // View if a function is defined in the file:
      if(currentSession.DefinesFunction("my_beauty_function"))
	{
	  LibUtilities::EquationSharedPtr pFun0;
	  pFun0 = currentSession.GetFunction("my_beauty_function","rho");
	  // Print the expression defined for "rho"
	  std::cout <<"Expression for rho :" <<  pFun0->GetExpression() << std::endl;


	  // With the equation specified we can evaluate basically everything we want:
	  double array[4] ={0.0,0.1,0.2,0.3};
	  Array<OneD,const NekDouble> x(4,array);
	  Array<OneD,const NekDouble> y(4,array);
	  Array<OneD,const NekDouble> z(4,array);
	  Array<OneD,NekDouble> res(x.size(),0.0);
	  pFun0->Evaluate(x,y,z,res);
	  for(const auto& el : res)
	    std::cout << el << std::endl;
	  
	}
      
 
    }
  catch(const std::runtime_error &)
    {
      return 1;
    }

  // Loop trough the parameters

  // Cout some interesting results
  
}
