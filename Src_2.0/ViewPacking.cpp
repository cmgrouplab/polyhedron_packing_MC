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



void Print_Packing(int LNx, int LNy, int LNz, Packing& PK) //print a single large input file for antiprism
{
  //cout<<"here1"<<endl;

  FILE* fp = fopen("Packing.View.off","w");

  //fprintf(fp, "\n\n\n");
  fprintf(fp, "OFF\n");
  fprintf(fp, "%d %d 0\n", PK.PackPoly[0].n_vertex*PK.N*LNx*LNy*LNz, PK.PackPoly[0].n_face*PK.N*LNx*LNy*LNz);


  //cout<<"here1"<<endl;
  
  for(int lx=0; lx<LNx; lx++)
    for(int ly=0; ly<LNy; ly++)
      for(int lz=0; lz<LNz; lz++)
	{
	  for(int n=0; n<PK.N; n++)
	    {
	      double TempCenterE[SD];
	      double TempCenterL[SD];
	      
	      for(int u=0; u<SD; u++)
		TempCenterE[u] = 0;
	      
	      TempCenterL[0] = PK.CenterL[n][0] + lx;
	      TempCenterL[1] = PK.CenterL[n][1] + ly;
	      TempCenterL[2] = PK.CenterL[n][2] + lz;

	      for(int u=0; u<SD; u++)
		for(int v=0; v<SD; v++)
		  TempCenterE[u] += PK.Lambda[u][v]*TempCenterL[v];

	      for(int v = 0; v<PK.PackPoly[n].n_vertex; v++)
		fprintf(fp, "%.12f  %.12f  %.12f\n", PK.PackPoly[n].Vertex[v][0]+TempCenterE[0], PK.PackPoly[n].Vertex[v][1]+TempCenterE[1], PK.PackPoly[n].Vertex[v][2]+TempCenterE[2]);
    
	    }
	}
  
  
  //cout<<"here2"<<endl;

  int temp_counter = 0;

  int num_vertex = PK.PackPoly[0].n_vertex;
  
  for(int lx=0; lx<LNx; lx++)
    for(int ly=0; ly<LNy; ly++)
      for(int lz=0; lz<LNz; lz++)
	{
	  int Cell_Index = temp_counter%(LNx*LNy*LNz);
	  printf("Cell_Index = %d\n", Cell_Index);

	  for(int n=0; n<PK.N; n++)
	    {

	      for(int f=0; f<PK.PackPoly[0].n_face; f++)
		{
		  fprintf(fp, "%d     ", PK.PackPoly[0].vert_num_fspecie[PK.PackPoly[0].face_type[f]]);
		  for(int h=0; h<PK.PackPoly[0].vert_num_fspecie[PK.PackPoly[0].face_type[f]]; h++)
		    fprintf(fp, "%d ", PK.PackPoly[0].Face[f][h]+num_vertex*(n+Cell_Index*PK.N));
		  fprintf(fp, "\n");
		}

	      ///$$$$$$$$$$$$$$$$$$$$####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
	    
	    }

	  temp_counter++;
	}
     
  

  fclose(fp);

}







int main(int argc, char **argv)
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout<<argv[0]<<endl<<argv[1]<<endl;
  //argv[0] should be the name of the exe...

  Packing PK(argv[1]); //this function always reads in Iconfig.txt
  double rho = PK.GetDensity();
  cout<<"The packing density is rho = "<<rho<<endl;


  cout<<"Overlap_flag = "<<PK.GlobalOverlapCheck(0)<<endl;

   //*************************************************************************

  int LNx, LNy, LNz;
  cout<<"please specify the system size:"<<endl;
  cout<<"LNx = "; cin>>LNx;
  cout<<"LNy = "; cin>>LNy;
  cout<<"LNz = "; cin>>LNz;
  
 

  Print_Packing(LNx, LNy, LNz, PK); //print a single large input file for antiprism

}
