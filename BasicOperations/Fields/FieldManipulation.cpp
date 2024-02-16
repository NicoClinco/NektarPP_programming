#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <vector>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/BasicUtils/FieldIOXml.h>
#include <FieldUtils/Field.hpp>
#include <MultiRegions/ExpList.h>

using namespace Nektar;


std::ostream& operator<<(std::ofstream& ostream,
			 const MultiRegions::ExpList& expansion)
{
  ostream << "This file contains information about an expansion list\n";
  ostream << "Total number of quadrature-points :"<< expansion.GetTotPoints() <<"\n";
  ostream << "Total number of degrees of freedom :"<< expansion.GetNcoeffs() << "\n";
  return ostream;
  //const auto& phys_array = expansion.GetPhys()
};


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
  
  /// Perform a cycle for every fld file contained in the .chk directory
  /// that is contained in the directory:
  /*
  for(auto pfd : FieldDefs)
    {
      auto fd = *pfd;
      ///Every single "FieldDef" contains different physical-fields:
      std::vector<std::string> field_names = fd.m_fields;

      // Field names:
      unsigned int counter =0;
      for(auto fieldname : field_names)
	{
	  std::cout << fieldname << " ";
	  ++counter;
	}
      //Here we have 14 fields.
      //std::cout << counter <<"\n";

    }
  */
  /// TO DO: Read the mesh file and introduce an ordering for
  /// understanding the coordinates.
  
  /// Printing the fields as vectors: size = 20:
  std::cout << "Field-data-size :" <<  FieldData.size() <<"\n";

  /// Note: Every field is stored in a contiguos vector
  /// for each fld file.
  /*
  for(const auto& phi : FieldData)
    std::cout << phi.size() << " ";
  std::cout << "\n";
  */
  
  FieldUtils::FieldSharedPtr pGenField
    = std::make_shared<FieldUtils::Field>(FieldUtils::Field());
  
  //Copy of the mesh-graph:
  pGenField->m_graph = SpatialDomains::MeshGraph::Read(pcurrentSession);

  
   
  //Resize the ExpList vector:
  //pGenField->m_exp.resize(1);
  pGenField->m_exp.emplace_back(new MultiRegions::ExpList(pcurrentSession,pGenField->m_graph));
  //pGenField->m_exp.push_back(pGenField->SetUpFirstExpList(0));
  //auto shared_ptrExp = pGenField->SetUpFirstExpList(1,true);
  //Create the expansion-list:
  //pGenField->m_exp[0] = pGenField->SetUpFirstExpList(0);
  
  std::ofstream infoExpansion("expansion_info.txt");

 
  
  infoExpansion << (*pGenField->m_exp[0]);
  infoExpansion.close();
  
  Array<OneD,double> global_coeffs(pGenField->m_exp[0]->GetNcoeffs());

  //Get the coefficients of the expansion:
  pGenField->m_exp[0]->ExtractCoeffsFromFile("TGV6p3_ALIGNED_1.chk",pComm,"rhou",global_coeffs);
  
  std::cout << "Number of coefficients of the solution :"<< global_coeffs.size() <<"\n";
  
  std::ofstream outputRho("rhou.csv");

  std::vector<Array<OneD,double>> coordinates(3,Array<OneD,double>(pGenField->m_exp[0]->GetTotPoints(),0.0));
  //pGenField->m_exp[0]->GetCoords(coordinates[0],coordinates[1],coordinates[2]);
  Array<OneD,double> coordX(pGenField->m_exp[0]->GetTotPoints());
  Array<OneD,double> coordY(pGenField->m_exp[0]->GetTotPoints());
  Array<OneD,double> coordZ(pGenField->m_exp[0]->GetTotPoints());
  pGenField->m_exp[0]->GetCoords(coordX,coordY,coordZ);
  
  size_t counter = 0;
  for(counter=0;counter<global_coeffs.size();++counter)
    {
      outputRho <<coordX[counter] <<"," << global_coeffs[counter] << ",\n";
    }
  outputRho.close();
   /*
  std::cout <<" ########################################### \n";
  std::cout << "#### TESTING HOW TO READ MULTIPLE FLD FILES ###\n";

  
  //Creating a specific class for read multifld files:
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
  std::cout << "Filenames : ";
  for(const auto& file : fileNames)
    std::cout << file << " ";
  std::cout << "\n";
  
  unsigned int total_elems = 0;  
  for(const auto& els : elementsID)
    total_elems+=els.size();
  

  // Check if we have printed all the elements:
  std::cout << total_elems << "\n";

  */
  
  
  
  return 0;
}
