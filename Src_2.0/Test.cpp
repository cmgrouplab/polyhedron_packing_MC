//test all new programmed classes...


using namespace std;

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "const.h"
#include "Geometry.h"
#include "Polyhedron.h"
//#include "Cells.h"
#include "Packing.h"
#include "Moves.h"
#include "ReadParam.h"
#include "Contact.h"




int main()
{
  //Test class Geometry.h
  /*
  //******************************************************************
  double vect1[SD], vect2[SD], vect3[SD];
  vect1[0] = vect1[1] = 1; vect1[2] = 1;
  vect2[0] = vect2[1] = 0; vect2[2] = 1;
  
  Geometry Comput;
  cout<<"Inner Product = "<<Comput.GetInnerProduct(vect1, vect2)<<endl;
  Comput.GetCrossProduct(vect1, vect2, vect3);
  cout<<vect3[0]<<vect3[1]<<vect3[2]<<endl;
  //******************************************************************
  */
  //*******************************************************************
  
  //Test class Polyhedron.h
  /*
  Polyhedron A("P1");
  Polyhedron B("P1");

  A.Vertex[0][0] = 0.0; A.Vertex[0][1] = 0.0; A.Vertex[0][2] = 0.0;
  A.Vertex[1][0] = 1.0; A.Vertex[1][1] = 0.0; A.Vertex[1][2] = 0.0;
  A.Vertex[2][0] = 0.0; A.Vertex[2][1] = 1.0; A.Vertex[2][2] = 0.0;
  A.Vertex[3][0] = 0.0; A.Vertex[3][1] = 0.0; A.Vertex[3][2] = 1.0;


  A.ShiftVert();
  A.GetNormal();
  A.GetCentDist();

  B = A;

  cout<<"B.Vertex[1] = ("<<B.Vertex[0][0]<<","<<B.Vertex[0][1]<<","<<B.Vertex[0][2]<<")"<<endl;
  cout<<"B.Normal[1] = ("<<B.Normal[1][0]<<","<<B.Normal[1][1]<<","<<B.Normal[1][2]<<")"<<endl;
  cout<<"B.d_cf = ("<<B.d_cf[0]<<","<<B.d_cf[1]<<","<<B.d_cf[2]<<","<<B.d_cf[3]<<")"<<endl;
  cout<<"B.d_cv = ("<<B.d_cv[0]<<","<<B.d_cv[1]<<","<<B.d_cv[2]<<","<<B.d_cv[3]<<")"<<endl;
  */
  //*************************************************************************
  /*
  //Test class Cells.h
  cout<<"Get here"<<endl;

  double Lambda[SD][SD] = {{1,0,0},{0,1,0},{0,0,1}};
  double CenterL[3][SD] = {{0.2, 0.4, 0.2},{0.4, 0.7, 0.9},{0.6, 0.2, 0.6}};

  Cells BoxCell;
  BoxCell.GetCellList(0, 0.3, Lambda, CenterL, 3, 0.5);
  cout<<"Ncell = "<<BoxCell.Ncell<<endl;
  BoxCell.PrintCellList();

  double OldC[SD] = {0.2, 0.4, 0.2}; double NewC[SD] = {0.7, 0.2, 0.4};
  BoxCell.UpdateCellList(0, OldC, NewC);
  BoxCell.PrintCellList();
  */
  //*************************************************************************
  //Test class Packing.h
  
  Packing PK("P1");
  PK.GetGlobalPosition();
  double rho = PK.GetDensity();
  cout<<"The packing density is rho = "<<rho<<endl;

  Cells BoxCell;
  BoxCell.GetCellList(0, PK.d0, PK.Lambda, PK.CenterL, PK.N, rho);
  cout<<"Ncell = "<<BoxCell.Ncell<<endl;
  BoxCell.PrintCellList();
 

  cout<<"Overlap_flag = "<<PK.GlobalOverlapCheck(0)<<endl;
  cout<<"Overlap_flag = "<<PK.GlobalOverlapCheck(BoxCell)<<endl;

  double overlap_flag = PK.GlobalOverlapCheck(0);
  while(overlap_flag == 1)
    {
      PK.Rescale(1.001);
      overlap_flag = PK.GlobalOverlapCheck(0);
      
    }
  cout<<"The packing density is rho = "<<PK.GetDensity()<<endl;
  
  /*
  Moves MC("P1");

  for(int i=0; i<100; i++)
    {
      int temp_index = rand()%PK.N;
      
      cout<<MC.ParticleMove(0.5, 0.1, 0.1, temp_index, PK, BoxCell)<<endl;

      cout<<"Overlap_flag = "<<PK.GlobalOverlapCheck(BoxCell)<<endl;
      cout<<"Overlap_flag = "<<PK.GlobalOverlapCheck(0)<<endl;
    }
  cout<<endl;
  MC.BoundaryMove(20, 0.15, 0.85, -0.01, 0.01, PK, BoxCell);
  cout<<"Overlap_flag = "<<PK.GlobalOverlapCheck(0)<<endl;
  cout<<"The packing density is rho = "<<PK.GetDensity()<<endl;

  ReadParam Param;
  cout<<"Nstage = "<<Param.Get_Nstage()<<endl;

  PK.PrintDensity();
  PK.PrintDensity(PK.GetDensity());
  */  

  //*************************************************************************
  //Test class Contact.h

  Contact Cont("P1", PK);

  Cont.GlobalGetContact();

  Cont.Print_Stat();

}
