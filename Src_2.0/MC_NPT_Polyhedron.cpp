//A MC quasi-NPT simulation of hard polyhedron system
//deformable boundary conditions and rombhedral unit cells are allowed
//using separating axies theorom as rigorous non-overlapping conditions
//ANY convex polyhedral shape can be handled by the code: simply program your shape in Polyhedron.h

//author: Yang JIAO, yjiao@princeton.edu
//started: July 1st, 2010

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using namespace std;

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "const.h"

#include "Geometry.h"
#include "Polyhedron.h"
#include "Packing.h"
#include "Moves.h"
#include "Pressure.h"
#include "Stat.h"
#include "ReadParam.h"


//Cells.h is included in Packing.h

int main(int argc, char **argv)
{
  //first define all the objects for the classes
  srand(time(NULL));

  ReadParam Param; //read the parameters first

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //cout<<argv[0]<<endl<<argv[1]<<endl;
  Packing PK(argv[1]);
  double rho = PK.GetDensity();
  cout<<"The initial packing density is rho_I = "<<rho<<endl;

 
  Cells BoxCell;
  BoxCell.GetCellList(0, PK.d0, PK.Lambda, PK.CenterL, PK.N, rho);
  cout<<"The total number of cells Ncell = "<<BoxCell.Ncell<<endl;
  BoxCell.PrintCellList();

  //cout<<"here 1"<<endl;

  //check overlapping of the intial configuration
  PK.Check_Iconfig(Param.Get_check_option(), Param.Get_relax_ratio());
  rho = PK.GetDensity();
  cout<<"The rescaled packing density (if applicable) is rho_R = "<<rho<<endl;

  //cout<<"here 2"<<endl;

  if(rho<Param.Get_Starting_Density())
    {
      printf("The initial packing is less dense than rho_start = %f\n", Param.Get_Starting_Density());
      printf("Start the simulation at rho_int = %f\n", rho);
    }
  else{
    while(rho>=Param.Get_Starting_Density())
      {
	PK.Rescale(Param.Get_relax_ratio());
	rho = PK.GetDensity();
      }
    printf("The specified starting density is rho_start = %f\n", Param.Get_Starting_Density());
    printf("Start the simulation at rho_int = %f\n", rho);
  }


  //for the statiscis ...
  Pressure Sample_Pressure;
  Stat Sample_Stat(1, PK);

  //*****************************************************************************
  //starting the MC simulations...


  //initialize SEP_List and NNL list
  // PK.GetNNL(Param.nnl_cut_dist); NNL is not rigorous
  //PK.GetSEP_List(argv[1]);


  Moves MC(argv[1]);
  
  printf("MC moves start...\n");

  //int method_flag = 0;
  
  for(int i=0; i<Param.Get_Nstage(); i++) //the stages...
    {
      printf("The %dth stage starts...\n", i+1);

    
      //initilize the property statistics
      //for the pressure
      if(Param.flag_pressure() == 1)
	{
	  Sample_Pressure.GetPressure();
	}   
    
      /* 
      if(rho > Param.rho_nnl_th && method_flag == 0)
	{
	 
	  cout<<"Switch to the SEP_List check...."<<endl;

	  method_flag = 1;
	  //PK.UpdateNNL(Param.nnl_cut_dist);

	  
	}
      */

      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      int sum_res = 0;
      int sample_start = (int)floor((double)Param.Get_Ncycle()*PK.N/2.0);

      for(int j=0; j<Param.Get_Ncycle()*PK.N; j++) //the cycles...
	{
	  //int index = rand()%PK.N; //this is randomly picking up th particles	  
	  int index = j%PK.N; //this is a linear sequential picking up...
	  
	  int res;
	  
	  res = MC.ParticleMove(Param.Get_p_trans(), Param.Get_TransMod(), Param.Get_RotMod(), index, PK, BoxCell);
	  //else if(method_flag == 1) res = MC.ParticleMoveII(Param.Get_p_trans(), Param.Get_TransMod(), Param.Get_RotMod(), index, PK, BoxCell);
	  //printf("%d", res);
	  if(res>0) sum_res++;
	 
	  //collect statistics after system has been equilibrated for a while, e.g., t>sample_start
	  if(j>sample_start && j%(Param.Get_Npc()*PK.N)==0)
	    {
	      
	      if(Param.flag_pressure() == 1)
		{
		  Sample_Pressure.GetPressure();
		}   
	      //could also sample other properties, like, g2 and order parameters... 
	    }

	} //Finishing a stage, i.e., a Ncycle number of sequential moves of each particle
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //compute the acc_rate and rescale the step size when necessary, adjust the dynamical parameters
      
      


      double R_Acc = (double)(sum_res+1)/(double)(Param.Get_Ncycle()*PK.N);
      
      printf("Acceptance Rate: Acc = %.10f\n", R_Acc);
      
      //permannetly adjust the step sizes if necessary...
      int Hnew_Acc = (int)log(1.0/R_Acc);
      printf("Hnew_Acc = %d\n", Hnew_Acc);
      //for adaptive step size, in both directions...
      if(Hnew_Acc > 1) //the step size is too large...
	{
	  Param.Rescale_TransMod();
	  Param.Rescale_RotMod();
	  Param.Rescale_Strain();
      	  //it is better to keep the StrainMod consistent with the trial move step size...
	  cout<<"The TransMod, RotMod and StrainMod are rescaled!"<<endl;
	}
   
      printf("\nThe %dth stage ends...\n", i+1);
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
      //Now we deform the boundary and report density...
      MC.BoundaryMove(Param.Get_Ndeform(),Param.Get_p_uphill(),Param.Get_Strain_Rescale(),Param.Get_Global_StrainMod(),Param.Get_Shear_StrainMod(),PK,BoxCell);
      //if we are unlucky, i.e., limit>Ndeform, we do not compress at the current stage...


      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //Now compute and report the density and other statistics for the current stage...

      rho = PK.GetDensity();
      PK.PrintDensity(rho);
      printf("The density of %dth stage: rho = %f\n", (i+1), rho);


      //the density change is larger than 5%, re-set the cell-list
      //cout<<"(rho, rho_cell) = "<<rho<<", "<<rho_cell<<endl;
      if(fabs(rho-BoxCell.rho_cell)>0.05)
	{
	  cout<<"(rho, rho_cell) = "<<rho<<", "<<BoxCell.rho_cell<<endl;
	  BoxCell.GetCellList(1, PK.d0, PK.Lambda, PK.CenterL, PK.N, rho);
	  cout<<"The CellList has been rebulit!"<<endl;

	  //check possible overlaps and rescale configuration if necessary
	  PK.Check_Iconfig(0, Param.Get_relax_ratio());
	}


      //get the statistics after each stage
      //for the pressure 
      if(Param.flag_pressure() == 1)
	{
	  Sample_Pressure.GetPressure();
	}   
      //for the g2
      if(Param.flag_g2() == 1)
	{
	  Sample_Stat.Get_PairCorr(PK);
	}
      //for the the order metric
      if(Param.flag_ordermetric() == 1)
	{
	  if(Param.flag_g2() == 1)
	    {
	      Sample_Stat.Get_LocalNormOrder(PK);
	      Sample_Stat.Get_LocalAllign(PK);
	    }
	  else if(Param.flag_g2() == 0)
	    {
	      Sample_Stat.Get_PairCorr(PK); //g2 has to be sampled
	      Sample_Stat.Get_LocalNormOrder(PK);
	      Sample_Stat.Get_LocalAllign(PK);
	    }
	}
      

      printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n");

      if(i%Param.Get_Npc()==0) PK.PrintPacking(); //print out the packing every 50 stages... 

      if(rho>=Param.Get_Terminate_Density())
	{
	  PK.PrintPacking();
	  printf("The density has reached rho = %f\n!", Param.Get_Terminate_Density());
	  printf("The program is terminated as requested!\n");
	  exit(1);
	}

      //PrintCellList();


    }//Finishing the current stage...

  //*****************************************************************************
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //post-processing...

  
  //compute the statistics for the final configuration...
  Sample_Stat.Get_PairCorr(PK); //g2 has to be sampled
  Sample_Stat.Get_LocalNormOrder(PK);
  Sample_Stat.Get_LocalAllign(PK);

  if(PK.GlobalOverlapCheck(1)==1)
    {
      printf("The Final Configuration Contains Overlapping Pairs!\n");
      exit(1);
    }


  PK.PrintPacking();



  return 0;
}


//EOF
