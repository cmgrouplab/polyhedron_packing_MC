//obtain statistics of the packing configurations (g2, allignment, face_normal, average contact)
//print out a file of particle centers in sphere packing format for other order metrics using Aleks codes


//author: Yang JIAO, yjiao@princeton.edu
//started: Dec. 9, 2010

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using namespace std;

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <time.h>
#include "const.h"

#include "Geometry.h"
#include "Polyhedron.h"
#include "Packing.h"
#include "Stat.h"
#include "Contact.h"


//Cells.h is included in Packing.h


//print the centroids as the LSD format, and using Aleks exe to obtain Q6...
void Write_LSD_Format(Packing &PK)
{
  ofstream fout;

  PK.GetGlobalPosition();

  fout.open("LSD_Packing.dat");

  fout<<"3 HS mono"<<endl;
  fout<<PK.N<<"\t 1"<<endl;
  fout<<PK.N<<endl;
  fout<<setprecision(12)<<PK.d0<<endl;

  for(int i=0; i<SD; i++)
    {
      for(int j=0; j<SD; j++)
	fout<<setprecision(12)<<PK.Lambda[j][i]<<"\t";
      fout<<endl;
    }

  fout<<"T  T  T"<<endl;

  for(int i=0; i<PK.N; i++)
    {
      fout<<setprecision(12)<<PK.CenterE[i][0]<<"\t";
      fout<<setprecision(12)<<PK.CenterE[i][1]<<"\t";
      fout<<setprecision(12)<<PK.CenterE[i][2]<<endl;	  
    }

  fout.close();
}



int main(int argc, char **argv)
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //cout<<argv[0]<<endl<<argv[1]<<endl;
  //argv[0] should be the name of the exe...
  Packing PK(argv[1]);
  double rho = PK.GetDensity();
  cout<<"The packing density is rho = "<<rho<<endl;


  cout<<"Overlap_flag = "<<PK.GlobalOverlapCheck(0)<<endl;

   //*************************************************************************
  //Get the contact numbers...

  Contact Cont(argv[1], PK);

  Cont.GlobalGetContact();

  Cont.Print_Stat();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //get statistics (g2, allignment, face_normal correlation...)
  int option = 0;
  cout<<"Comput Correlation functions? Yes - 1; No - 0 "<<endl;
  cin>>option;
  if(option == 1)
    {

      Stat Sample_Stat(1, PK);
      
      Sample_Stat.Get_PairCorr(PK);
      
      cout<<"~~~~~~~~~~~~~~~~~~~~"<<endl;

      Sample_Stat.Get_LocalNormOrder(PK);

      cout<<"~~~~~~~~~~~~~~~~~~~~"<<endl;

      Sample_Stat.Get_LocalAllign(PK);

      
    }
 

  cout<<"Print out LSD format for Q6? Yes - 1; No - 0 "<<endl;
  cin>>option;
  if(option == 1)
    {
      Write_LSD_Format(PK);
    }

}
