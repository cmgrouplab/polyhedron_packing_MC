//A MC simulation of hard tetrahedron system
//deformable boundary conditions and rombhedral unit cells are allowed
//using separating axies theorom as rigorous non-overlapping conditions
//subroutines dealing with the geometries of tetrahedron are developed



//author: Yang JIAO, yjiao@princeton.edu
//started: Jun. 15, 2010
//based on tetrahedron code MC_NVT_Tetrah.16.v2.cpp
//      and included all the updated features
//      (i) sample g2 (ii) sample pressure <not accurate at this stage> 
//      (iii) improved overlap_check(), which only looks at immediate neighbors
//Further more, new features added
//      (iv) all parameters can be read in or default
//      (v) allow RSA inital configurations
//      (vi) the partition of cells are adaptive, i.e., the number of cell along
//           each dimension is computed as N_x = max{[L_x/(1.5*d_max)], 3}, with
//           the density rho_part at current partition saved. If the density is 
//           increased by 5%, the partition is re-done


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using namespace std;

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "GEsolver.h"
#include "const.h" 

#define MAXX 50000
#define MAXY 100000
#define Delta 0.00001

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int N = 300; //number of particles, to be read in latter unless default will be used

int Nstage = 10; //on each the box is fixed, at the end of a stage deformation is performed
int Ncycle = 100; //number of cycles per stage
int Ndeform = 20; //number of deformation performed at the end of each stage


//specify the magnitude of the displacements...
double pi = 3.141592654;
double TransMod = 0.1; //associated with the relative coordinates...no unity needs to be specified...
double RotMod = 0.1*pi;
double Global_StrainMod = -0.001; // the default strain mod: >0 is expansion and <0 is compression
double Shear_StrainMod = 0.00; //separate shear strain and principle strain

double Strain_Rescale = 0.75; //if compression fails, rescale the strain mod by this amount each time
double Trans_Rescale = 0.5; //for adpative step sizes, make sure acceptance rate is around 0.5
double Rot_Rescale = 0.5;
//when rescaling the step size, also rescaling the strain...
//so the specified step size and strain should be comparable, otherwise the real strain would small than expected

double p_uphill = 0.35; //the generation rate of uphill moves...
double p_trans = 0.5; //probabiliy of the current move is a translation 

double relax_ratio = 1.001; //if the initial configuration contains overlapping, rescale the lattice

double Terminate_Density = 0.6;
double Starting_Density = 0.6;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double d0 = 1.0;  
//for random tetrahedron, this is a charactersitic length
//   useful when determining the cell numbers
//   normally choose to be the largest edge length, computed when read in the config.
                

struct Tetrah{
  double Vertex[4][SD];
  //the coordinates of the 4 vertex of a tetrah. relative to the center O = 0
  //the order is A, B, C, D

  double Normal[4][SD];
  //unit outward normal of the four faces of a tetrah.
  //the order is ABC, ABD, ACD, BCD

  double d_cf[4];
  //for a general shape, this gives the 4 center to face distance for overlap check
  //   it's computed when read in the config.
  //   by doing this, we allow each particle possing a different shape

  double d_cv[4];
  //this is the centroid to vertex distances, useful for overlapping check
};

double Lambda[SD][SD]; //the lattice vectors...
double Strain[SD][SD]; //the strain tensor...
double Old_Lambda[SD][SD]; //the old lattice vectors...

double** CenterE; // the centers in Eculdean coordinates...
double** CenterL; // the centers in Relative coordinates...

Tetrah* PackTetrah; //provide the orentations...

//temp varables for the MC moves...
Tetrah Old_Tetrah;
double Old_CenterL[SD]; //store the old position info. of the moved particle


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//for the cells

//the number of cells along each direction, determined automatically later
//    these are also updated adaptively
int LXcell = 3; 
int LYcell = 3;
int LZcell = 3;
int Ncell; //the total number of cells

double rho_cell; //the density associated with the last cell update

struct node{
  int index_Tetrah; //the index of the Tetrah. in a certain cell...
  node* next;
};

node** CellList; //the cell list containing the particles in the cell...
//the corner locations are the index...
//implemented as a one-dimensional array

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//for statistics

int Npc = 10; //the number of trial moves between each property collection

//for the pressure
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int flag_comput_pressure = 0; //by default, do not collect pressure

double Stage_Pressure = 0; 
int pressure_ct = 0;

double GapMod = 0.00005;
//this should be a very small constant
//delta_V = 3*GapMod
double V_sys = 0; //the volume of the system...
double rho_s;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//for g_2

int flag_comput_g2 = 0;

double N_bin = 8.0; //the number of bins per insphere radius
double bin; //the bin width...
double MaxL; //the characteristic length of the simulation box, i.e., the length of the shortest lattice vector
int NLcounter; //the actual sample length
#define MaxNL 200 // for the sampling length, the upper limit

//the pair distribution function...
double g[MaxNL];
double g_ave[MaxNL]; //get the sum and average at the end of each stage
int g2_ct = 0;//the times g2 is sampled per stage...
double rho_n; //number density

//*********************************************************************
//declaration of all the functions...









//*********************************************************************
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//subroutines for pre-processing....

void Int_Data()
{
  cout <<" Please specify all the parameters for the packing program: "<<endl;

  cout <<" The total number of particles -  N := "; cin >> N; cout<<endl;
  cout <<" The number of stages -  Nstage := "; cin >> Nstage; cout<<endl;
  cout <<" The number of cycles per stage - Ncycle := "; cin >> Ncycle; cout<<endl;
  cout <<" The number of deformation per stage - Ndeform := "; cin>>Ndeform; cout<<endl;

  cout <<" The translation magnitude - TransMod := "; cin>>TransMod; cout<<endl;
  cout <<" The rotation magnitude (in unit pi) - RotMod := "; cin>>RotMod; RotMod = RotMod*pi; cout<<endl;
  cout <<" Normal strain magnitude - Global_StrainMod := "; cin>>Global_StrainMod; cout<<endl;
  cout <<" Shear strain magnitude - Shear_StrainMod := "; cin>>Shear_StrainMod; cout<<endl;

  cout <<" Strain rescale magnitude - Strain_Rescale := "; cin>>Strain_Rescale; cout<<endl;
  cout <<" Translation rescale magnitude - Trans_Rescale := "; cin>>Trans_Rescale; cout<<endl;
  cout <<" Rotation rescale magnitude - Rot_Rescale := "; cin>>Rot_Rescale; cout<<endl;

  cout <<" Probability of uphill moves - p_uphill := "; cin>>p_uphill; cout<<endl;
  cout <<" Probability of translation for a trail move - p_trans := "; cin>>p_trans; cout<<endl;

  cout <<" Rescale paramter for initial configuration - relax_ratio :="; cin>>relax_ratio; cout<<endl;
  cout <<" Termination density - Terminate_Density := "; cin>>Terminate_Density; cout<<endl;
  cout <<" Starting density - Starting_Density := "; cin>>Starting_Density; cout<<endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
  cout <<" The number of trial moves between each property collection - Npc := "; cin>>Npc; cout<<endl;

  cout <<" Collecting Pressure? 1 - Yes; 0 - No." <<endl;
  cin >> flag_comput_pressure; cout<<endl;
  if(flag_comput_pressure == 1)
    {
      cout <<" Pressure will be collected and control paramters are tunable in source code" <<endl;
    }

  cout <<" Computing g2? 1 - Yes; 0 - No." <<endl;
  cin >> flag_comput_g2; cout<<endl;
  if(flag_comput_g2 == 1)
    {
      cout <<" g2 will be computed and control paramters are tunable in source code" <<endl;
    }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //now initalize the variables

  CenterE = new double*[N];
  CenterL = new double*[N];

  PackTetrah = new Tetrah[N];
  
  for(int i=0; i<N; i++)
    {
      CenterE[i] = new double[SD];

      CenterL[i] = new double[SD];

    }

}




void Rescale() //if packing containing overlapping pairs, rescale the packing
{
  for(int i=0; i<SD; i++)
    for(int j=0; j<SD; j++)
      Lambda[i][j] = relax_ratio*Lambda[i][j];
}

void GetStrainIso(double P_StrainMod, double S_StrainMod) //modified to allow uphill moves...
{
  
  //double Rand = (double)(rand()%MAXX)/(double)MAXX;
  //Strain[0][0] = Strain[1][1] = Strain[2][2] = fabs(P_StrainMod); 
  //this guarantee the deformation is isotropic, while cannot fully explore the strain space....
  //only for isotropic expansion...

  double Rand;

  Rand = (double)(rand()%MAXX)/(double)MAXX;
  Strain[0][0] = (Rand-p_uphill)*P_StrainMod;

  Rand = (double)(rand()%MAXX)/(double)MAXX;
  Strain[1][1] = (Rand-p_uphill)*P_StrainMod;

  Rand = (double)(rand()%MAXX)/(double)MAXX;
  Strain[2][2] = (Rand-p_uphill)*P_StrainMod;


  //Sal's suggestion, double the size of the uphill moves at middle stages..
  Rand = (double)(rand()%MAXX)/(double)MAXX;
  if(Rand<0.0005)//use supper large uphill moves
    {
      Strain[0][0] = -10.0*fabs(Strain[0][0]);
      Strain[1][1] = -10.0*fabs(Strain[1][1]);
      Strain[2][2] = -10.0*fabs(Strain[2][2]);      
    }


  Rand = (double)(rand()%MAXX)/(double)MAXX;

  Strain[0][1] = Strain[1][0] = (Rand-0.5)*S_StrainMod;

  Rand = (double)(rand()%MAXX)/(double)MAXX;

  Strain[0][2] = Strain[2][0] = (Rand-0.5)*S_StrainMod;

  Rand = (double)(rand()%MAXX)/(double)MAXX;

  Strain[2][1] = Strain[1][2] = (Rand-0.5)*S_StrainMod;
}


void ResetStrain(int H) //decrease the magnitude by certain amount..
{
  double Temp_PStrainMod = pow(Strain_Rescale, H)*Global_StrainMod;
  double Temp_SStrainMod = pow(Strain_Rescale, H)*Shear_StrainMod;

  GetStrainIso(Temp_PStrainMod, Temp_SStrainMod);

}




void GetNormal(int);

//no need to change the configuration format, the full geo charact. is contained in ths shape defined via vertex
void Read_Config()
{
  FILE* fp;

  if((fp=fopen("Mconfig.txt","r"))==NULL)
    {
      if((fp=fopen("Iconfig.txt","r"))==NULL)
	{
	  printf("Cannot open file Iconfig.txt!\n");
	  exit(1);
	}
      else //read the initial configuration...
	{
	  int tempm;
	  double dI;
	  double val, xt, yt, zt;
	  fscanf(fp, "%d", &tempm); // read in the number of particles...
	  if(tempm!=N) 
	    {
	      printf("Bad particle numbers! Recheck!\n");
	      //exit(1);
	    }
	  fscanf(fp, "%lf", &dI); // read in the length...
	  d0 = dI; //should be the maximal length for all tetrah.
	  
	  for(int i=0; i<SD; i++) // read in the lattice vectors...
	    for(int j=0; j<SD; j++)
	      {
		fscanf(fp, "%lf", &val);
		Lambda[j][i] = val;
	      }

	  for(int i=0; i<N; i++) // read the posistions of the center and the vertex
	    {
	      //read in the center of mass
	      fscanf(fp, "%lf", &val); CenterL[i][0] = val; 
	      if(CenterL[i][0]>=1.0) CenterL[i][0] = CenterL[i][0]-1.0;
	      else if(CenterL[i][0]<0) CenterL[i][0] = CenterL[i][0]+1.0;
	      fscanf(fp, "%lf", &val); CenterL[i][1] = val;
	      if(CenterL[i][1]>=1.0) CenterL[i][1] = CenterL[i][1]-1.0;
	      else if(CenterL[i][1]<0) CenterL[i][1] = CenterL[i][1]+1.0;
	      fscanf(fp, "%lf", &val); CenterL[i][2] = val;
	      if(CenterL[i][2]>=1.0) CenterL[i][2] = CenterL[i][2]-1.0;
	      else if(CenterL[i][2]<0) CenterL[i][2] = CenterL[i][2]+1.0;

	   
	      //read in vertex with O = 0;
	      //for vertex A...
	      for(int j=0; j<SD; j++)
		{
		  fscanf(fp, "%lf", &xt);
		  PackTetrah[i].Vertex[0][j] = xt;
		}
              //for vertex B...
	      for(int j=0; j<SD; j++)
		{
		  fscanf(fp, "%lf", &xt);
		  PackTetrah[i].Vertex[1][j] = xt;
		}
              //for vertex C...
	      for(int j=0; j<SD; j++)
		{
		  fscanf(fp, "%lf", &xt);
		  PackTetrah[i].Vertex[2][j] = xt;
		}
              //for vertex D...
	      for(int j=0; j<SD; j++)
		{
		  fscanf(fp, "%lf", &xt);
		  PackTetrah[i].Vertex[3][j] = xt;
		}

	      //we shift the particle so that the center is at the origin of the local system
	      double PO[SD];
	      for(int j=0; j<SD; j++)
		{
	   	  PO[j] = (PackTetrah[i].Vertex[0][j]+PackTetrah[i].Vertex[1][j]+PackTetrah[i].Vertex[2][j]+PackTetrah[i].Vertex[3][j])/4.0;
		 }

	      for(int j=0; j<4; j++)
		{
		  double temp_dist2 = 0.0;

		  for(int v=0; v<SD; v++)
		    {
		      PackTetrah[i].Vertex[j][v] = PackTetrah[i].Vertex[j][v] - PO[v];
		      
		      temp_dist2 += PackTetrah[i].Vertex[j][v]*PackTetrah[i].Vertex[j][v];
		    }

		  PackTetrah[i].d_cv[j] = sqrt(temp_dist2);

		  //cout << "Tetrah["<<i<<"].d_cv["<<j<<"] = "<<PackTetrah[i].d_cv[j]<<endl;
		}

	      //compute the normal 
	      //also compute the ceter to face distance
	      GetNormal(i);

	    }

	  fclose(fp);
	}

    }
  else // read in the intermediate configruations...
    { 
      printf("No Mconfig.txt file is needed here. Change or check the files!\n");  
      
    }
}


void GetGlobalPosition()
{
  //compute the global posistions of the center of mass
  for(int n=0; n<N; n++)
    {
      //for(int i=0; i<SD; i++)
      //CenterE[n][i] = 0;

      for(int i=0; i<SD; i++)
	{
	  //for(int j=0; j<SD; j++)
	    CenterE[n][i] = Lambda[i][0]*CenterL[n][0]+Lambda[i][1]*CenterL[n][1]+Lambda[i][2]*CenterL[n][2];
	}
    
    }
}

int NumMax(int a, int b)
{
  if(a>b) return a;
  else return b;
}

void GetCellList(int flag)
{
  
  if(flag != 0) //not the first time to establish CellList, so need to clean the list first
    {
      for(int i=0; i<Ncell; i++)
	{ 
	  while(CellList[i]!=NULL)
	    {
	      node* temp_pt = CellList[i];
	      CellList[i] = CellList[i]->next;

	      delete temp_pt;
	    }
	}
    }

  //now need to compute the number of cell along each direction
  double LatX, LatY, LatZ; //the length of lattice vectors along each direction
  
  LatX = sqrt(Lambda[0][0]*Lambda[0][0]+Lambda[1][0]*Lambda[1][0]+Lambda[2][0]*Lambda[2][0]);
  LatY = sqrt(Lambda[0][1]*Lambda[0][1]+Lambda[1][1]*Lambda[1][1]+Lambda[2][1]*Lambda[2][1]);
  LatZ = sqrt(Lambda[0][2]*Lambda[0][2]+Lambda[1][2]*Lambda[1][2]+Lambda[2][2]*Lambda[2][2]);

  LXcell = (int)NumMax(3, (int)floor(LatX/(1.2*d0)));
  LYcell = (int)NumMax(3, (int)floor(LatY/(1.2*d0)));
  LZcell = (int)NumMax(3, (int)floor(LatZ/(1.2*d0)));
  //the size of cell has to be large enough such that next neighbor cell can never contain possibly overlapped particles

  Ncell = LXcell*LYcell*LZcell;
  
  cout<<"The total number of cells Ncell = "<<Ncell<<endl;

  CellList = new node*[Ncell]; //re-setup the cell list
  for(int i=0; i<Ncell; i++)
    CellList[i] = NULL;

  for(int n=0; n<N; n++)
    {
      int temp_indexI = floor(CenterL[n][0]*LXcell);
      int temp_indexJ = floor(CenterL[n][1]*LYcell);
      int temp_indexK = floor(CenterL[n][2]*LZcell);

      if(temp_indexI>=LXcell || temp_indexJ>=LYcell || temp_indexK>=LZcell)
	{
	  printf("CenterL includes invalid/un-scaled coordinates! Recheck!\n");
	  exit(1);
	}

      int cell_index = LXcell*LYcell*temp_indexK + LXcell*temp_indexJ + temp_indexI;

      //cout<<"Tetrah"<<n<<" cell_index = "<<cell_index<<endl;

      node* pt = new node;
      pt->index_Tetrah = n;
      pt->next = CellList[cell_index];
      CellList[cell_index] = pt;

    }
}


//subroutines dealing with geometries of the particles
//********************************************************
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GetCrossProduct(double Vect1[SD], double Vect2[SD], double Product[SD])
{
  Product[0] = Vect1[1]*Vect2[2] - Vect1[2]*Vect2[1];
  Product[1] = Vect1[2]*Vect2[0] - Vect1[0]*Vect2[2];
  Product[2] = Vect1[0]*Vect2[1] - Vect1[1]*Vect2[0];

}

double GetInnerProduct(double Vect1[SD], double Vect2[SD])
{
  double sum = 0;

  for(int i=0; i<SD; i++)
    sum += Vect1[i]*Vect2[i];

  return sum;
}

//get the common 
void GetPLine(double Direct[SD], double Vect[SD], double PLine[SD])
{
  double sum = sqrt(Direct[0]*Direct[0]+Direct[1]*Direct[1]+Direct[2]*Direct[2]);
  for(int k=0; k<SD; k++)
    Direct[k] = Direct[k]/sum;

  double temp_length = GetInnerProduct(Direct, Vect);

  for(int k=0; k<SD; k++)
    Direct[k] = temp_length*Direct[k];

  for(int k=0; k<SD; k++)
    PLine[k] = Vect[k] - Direct[k];
}

//the parameter is the index of a tetrahedron...
//the normal is unitary
void GetNormal(int m)
{
  double v1[SD], v2[SD];
  double temp_norm[SD];
  double temp_vect[SD];

  //for ABC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(int i=0; i<SD; i++)
    {
      v1[i] = PackTetrah[m].Vertex[0][i] - PackTetrah[m].Vertex[1][i]; //v1 = A-B
      v2[i] = PackTetrah[m].Vertex[0][i] - PackTetrah[m].Vertex[2][i]; //v2 = A-C
    }

  GetCrossProduct(v1, v2, temp_norm);

  for(int i=0; i<SD; i++)
    {
      //PackTetrah[m].Normal[0][i] = temp_norm[i];
      temp_vect[i] = PackTetrah[m].Vertex[3][i]; // the rest vertex...
    }

  if(GetInnerProduct(temp_vect, temp_norm)>0) //this is the inner direction....
    {
      for(int i=0; i<SD; i++)
	PackTetrah[m].Normal[0][i] = -temp_norm[i];
    }
  else //make sure it is the outward normal...
    {
      for(int i=0; i<SD; i++)
	PackTetrah[m].Normal[0][i] = temp_norm[i];
    }
 
  //for ABD~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(int i=0; i<SD; i++)
    {
      v1[i] = PackTetrah[m].Vertex[0][i] - PackTetrah[m].Vertex[1][i]; //v1 = A-B
      v2[i] = PackTetrah[m].Vertex[0][i] - PackTetrah[m].Vertex[3][i]; //v2 = A-D
    }

  GetCrossProduct(v1, v2, temp_norm);

  for(int i=0; i<SD; i++)
    {
      //PackTetrah[m].Normal[0][i] = temp_norm[i];
      temp_vect[i] = PackTetrah[m].Vertex[2][i]; // the rest vertex...
    }

  if(GetInnerProduct(temp_vect, temp_norm)>0) //this is the inner direction....
    {
      for(int i=0; i<SD; i++)
	PackTetrah[m].Normal[1][i] = -temp_norm[i];
    }
  else //make sure it is the outward normal...
    {
      for(int i=0; i<SD; i++)
	PackTetrah[m].Normal[1][i] = temp_norm[i];
    }


 //for ACD~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(int i=0; i<SD; i++)
    {
      v1[i] = PackTetrah[m].Vertex[0][i] - PackTetrah[m].Vertex[2][i]; //v1 = A-C
      v2[i] = PackTetrah[m].Vertex[0][i] - PackTetrah[m].Vertex[3][i]; //v2 = A-D
    }

  GetCrossProduct(v1, v2, temp_norm);

  for(int i=0; i<SD; i++)
    {
      //PackTetrah[m].Normal[0][i] = temp_norm[i];
      temp_vect[i] = PackTetrah[m].Vertex[1][i]; // the rest vertex...
    }

  if(GetInnerProduct(temp_vect, temp_norm)>0) //this is the inner direction....
    {
      for(int i=0; i<SD; i++)
	PackTetrah[m].Normal[2][i] = -temp_norm[i];
    }
  else //make sure it is the outward normal...
    {
      for(int i=0; i<SD; i++)
	PackTetrah[m].Normal[2][i] = temp_norm[i];
    }


 //for BCD~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(int i=0; i<SD; i++)
    {
      v1[i] = PackTetrah[m].Vertex[1][i] - PackTetrah[m].Vertex[2][i]; //v1 = B-C
      v2[i] = PackTetrah[m].Vertex[1][i] - PackTetrah[m].Vertex[3][i]; //v2 = B-D
    }

  GetCrossProduct(v1, v2, temp_norm);

  for(int i=0; i<SD; i++)
    {
      //PackTetrah[m].Normal[0][i] = temp_norm[i];
      temp_vect[i] = PackTetrah[m].Vertex[0][i]; // the rest vertex...
    }

  if(GetInnerProduct(temp_vect, temp_norm)>0) //this is the inner direction....
    {
      for(int i=0; i<SD; i++)
	PackTetrah[m].Normal[3][i] = -temp_norm[i];
    }
  else //make sure it is the outward normal...
    {
      for(int i=0; i<SD; i++)
	PackTetrah[m].Normal[3][i] = temp_norm[i];
    }

  //Now change them to unit normal...
  for(int i=0; i<4; i++)
    {
      double temp_length = 0;

      for(int j=0; j<SD; j++)
	temp_length += PackTetrah[m].Normal[i][j]*PackTetrah[m].Normal[i][j];

      temp_length = sqrt(temp_length);

      for(int j=0; j<SD; j++)
	PackTetrah[m].Normal[i][j] = PackTetrah[m].Normal[i][j]/temp_length;
    }


  //now get the center2face distance
  //for ABC
  PackTetrah[m].d_cf[0] = fabs(GetInnerProduct(PackTetrah[m].Vertex[0], PackTetrah[m].Normal[0]));

  //for ABD
  PackTetrah[m].d_cf[1] = fabs(GetInnerProduct(PackTetrah[m].Vertex[0], PackTetrah[m].Normal[1]));
 
  //for ACD
  PackTetrah[m].d_cf[2] = fabs(GetInnerProduct(PackTetrah[m].Vertex[0], PackTetrah[m].Normal[2]));
 
  //for BCD
  PackTetrah[m].d_cf[3] = fabs(GetInnerProduct(PackTetrah[m].Vertex[1], PackTetrah[m].Normal[3]));
 

}




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//*****************************************************************
//subroutines related to the MC moves...


double MinDis(int m, int n, double temp_disE[SD], int indexI, int indexJ, int indexK)
{
  //find the minimal distance between the centers of two tetrah. in Ecludean space...
  //by checking all the images of Tetrah[m], while keeping Tetrah[n] in the central box
  //record the index of the box which Tetrah[m] is in...
  //the order (m,n) is important, especially when getting the NNL.... 
  double tempLd[SD];
  tempLd[0] = (CenterL[m][0]-CenterL[n][0])+indexI;
  tempLd[1] = (CenterL[m][1]-CenterL[n][1])+indexJ;
  tempLd[2] = (CenterL[m][2]-CenterL[n][2])+indexK;
  
  for(int q=0; q<SD; q++)
    temp_disE[q] = 0;
  
  for(int q=0; q<SD; q++)
    for(int h=0; h<SD; h++)
      temp_disE[q] += Lambda[q][h]*tempLd[h];
  
  double dist = temp_disE[0]*temp_disE[0]+temp_disE[1]*temp_disE[1]+temp_disE[2]*temp_disE[2]; 

  //printf("Dist = %f\n", dist);

  return sqrt(dist); //this is center-to-center distance...
}



double MinDis(int m, int n, double temp_disE[SD])
{
  //find the minimal distance between the centers of two tetrah. in Ecludean space...
  //by checking all the images of Tetrah[m], while keeping Tetrah[n] in the central box
  //record the index of the box which Tetrah[m] is in...
  //the order (m,n) is important, especially when getting the NNL.... 
  double dx = (CenterL[m][0]-CenterL[n][0]);
  double dy = (CenterL[m][1]-CenterL[n][1]);
  double dz = (CenterL[m][2]-CenterL[n][2]);

  double dist = 10000000.0; //just a large number...


  //loop over all possible images of Point m, keep n fixed in the center simulation box....
  for(int i=-1; i<=1; i++)
    for(int j=-1; j<=1; j++)
      for(int k=-1; k<=1; k++)
	{
	  double tempLd[SD]; 
	  tempLd[0] = dx + i;
	  tempLd[1] = dy + j;
	  tempLd[2] = dz + k;

	  double tempGd[SD];
	  for(int q=0; q<SD; q++)
	    tempGd[q] = 0;

	  for(int q=0; q<SD; q++)
	    for(int h=0; h<SD; h++)
	      tempGd[q] += Lambda[q][h]*tempLd[h];


	  //including the translational displacement...
	  double tempdist = tempGd[0]*tempGd[0]+tempGd[1]*tempGd[1]+tempGd[2]*tempGd[2]; 

	  //printf("TempDist[%d][%d][%d] = %f\n", i, j, k, tempdist);

	  if(tempdist<dist) // store the smallest distance...
	    {
	      dist = tempdist;

     	      //IndexI = i; IndexJ = j; IndexK = k;

	      for(int q=0; q<SD; q++)
		temp_disE[q] = tempGd[q];
	      
	    }

	}

  // printf("IMAmndx = %f\t IMAmndy = %f\t IMAmndz = %f\n", IMAmn_Ldx, IMAmn_Ldy, IMAmn_Ldz);

  //printf("Dist = %f\n", dist);

  return sqrt(dist); //this is center-to-center distance...
}








double GetMin(double Arr[], int Size)
{
  double temp_min = 1000000000.0;

  for(int i=0; i<Size; i++)
    {
      if(Arr[i]<temp_min) temp_min = Arr[i];
    }

  return temp_min;
}

double GetMax(double Arr[], int Size)
{
  double temp_max = -1000000000.0;

  for(int i=0; i<Size; i++)
    {
      if(Arr[i]>temp_max) temp_max = Arr[i];
    }

  return temp_max;
}




int CheckOverlap(int m, int n, int indexI, int indexJ, int indexK) //this is a higly non-trivial task...
{
  double temp_VdisE[SD]; //the relative and Ecludean separation distance...
  
  //first get the center-to-center distance of the particles...
  double temp_dis = MinDis(m, n, temp_VdisE, indexI, indexJ, indexK);
  //recall n is fixed in the central box, the images of m are checked...

  //the mininmal and maximal possible separations
  double d_cvmax = (GetMax(PackTetrah[m].d_cv, 4)+GetMax(PackTetrah[n].d_cv, 4));
  double d_cfmin = (GetMin(PackTetrah[m].d_cf, 4)+GetMin(PackTetrah[n].d_cf, 4));

  if(temp_dis>d_cvmax || temp_dis < Delta) return 0; //greater than largest possible separation, definitely not overlapping or checking the particle itself, i.e., m = n
  else if(temp_dis<d_cfmin &&  temp_dis>=Delta) return 1; //smaller than the smallest possible separation, definitely overlapping


  //if in the intermediate distance, need to check with the separation axis theorem
  //need to do the complete check if overlapping: the most inefficient cases...

 
  //Now we work with TetrahM and TetrahN, equvalently the image of m and n; NO further modifications concerning the global coordiations of the vertex and normal are needed
  
  //*******************************************************************
  int overlap_flag = 1;

  //first, check the faces of TetrahN, i.e., loop over all faces of N and vertex of M...
  for(int i=0; i<4; i++) //the face normal
    {
      double temp_dist[4];
     
      for(int j=0; j<4; j++) //the vertex
	{
	  double temp_vertex[SD];

	  for(int k=0; k<SD; k++)
	    temp_vertex[k] = PackTetrah[m].Vertex[j][k] + temp_VdisE[k];

	  temp_dist[j] = GetInnerProduct(PackTetrah[n].Normal[i], temp_vertex);
	}

      //each face may have a different center2face distance, make sure they match
      if(GetMin(temp_dist, 4)>PackTetrah[n].d_cf[i]) //the particles do not overlap...
	{
	  overlap_flag = 0;
	  return overlap_flag;
	}
    }//end of loop for the faces...


  //second, check the faces of TetrahM
  
  for(int i=0; i<4; i++)
   {
     double temp_dist[4];
   
      for(int j=0; j<4; j++) //the vertex
	{
	  double temp_vertex[SD];

	  for(int k=0; k<SD; k++)
	    temp_vertex[k] = PackTetrah[n].Vertex[j][k] - temp_VdisE[k];

	  temp_dist[j] = GetInnerProduct(PackTetrah[m].Normal[i], temp_vertex);
	}

      if(GetMin(temp_dist, 4)>PackTetrah[m].d_cf[i]) //the particles do not overlap...
	{
	  overlap_flag = 0;
	  return overlap_flag;
	}

   }//end of loop for the faces...


  //*************************************************************
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //next we need to check all pair of edges to find the separating axis...

  double edgeM[6][SD], edgeN[6][SD];

  for(int j=0; j<SD; j++)
    {
      edgeM[0][j] = PackTetrah[m].Vertex[0][j] - PackTetrah[m].Vertex[1][j];
      edgeM[1][j] = PackTetrah[m].Vertex[0][j] - PackTetrah[m].Vertex[2][j];
      edgeM[2][j] = PackTetrah[m].Vertex[0][j] - PackTetrah[m].Vertex[3][j];
      edgeM[3][j] = PackTetrah[m].Vertex[1][j] - PackTetrah[m].Vertex[2][j];
      edgeM[4][j] = PackTetrah[m].Vertex[1][j] - PackTetrah[m].Vertex[3][j];
      edgeM[5][j] = PackTetrah[m].Vertex[2][j] - PackTetrah[m].Vertex[3][j];

      edgeN[0][j] = PackTetrah[n].Vertex[0][j] - PackTetrah[n].Vertex[1][j];
      edgeN[1][j] = PackTetrah[n].Vertex[0][j] - PackTetrah[n].Vertex[2][j];
      edgeN[2][j] = PackTetrah[n].Vertex[0][j] - PackTetrah[n].Vertex[3][j];
      edgeN[3][j] = PackTetrah[n].Vertex[1][j] - PackTetrah[n].Vertex[2][j];
      edgeN[4][j] = PackTetrah[n].Vertex[1][j] - PackTetrah[n].Vertex[3][j];
      edgeN[5][j] = PackTetrah[n].Vertex[2][j] - PackTetrah[n].Vertex[3][j];
    }


  for(int i=0; i<6; i++)
    for(int j=0; j<6; j++)
      {
	double temp_axis[SD];
	
	GetCrossProduct(edgeM[i], edgeN[j], temp_axis);

	//see whether the pair of edges are parallel...
	double temp_sum = 0;

	for(int k=0; k<SD; k++)
	  temp_sum += temp_axis[k]*temp_axis[k];

	//dealing with the parallel case
	if(temp_sum==0)
	  {
	    double temp_vect1[SD], temp_vect2[SD];
	    //for Tetrah M
	    if(i<=2)
	      {
		for(int v=0; v<SD; v++)
		  temp_vect1[v] = PackTetrah[m].Vertex[0][v]+temp_VdisE[v]; 
	      }
	    else if(i==3)
	      {
        	for(int v=0; v<SD; v++)
		  temp_vect1[v] = PackTetrah[m].Vertex[1][v]+temp_VdisE[v]; 
	      }
            else if(i>3)
	      {
        	for(int v=0; v<SD; v++)
		  temp_vect1[v] = PackTetrah[m].Vertex[3][v]+temp_VdisE[v]; 
	      }
	    //for Tetrah N
            if(j<=2)
	      {
		for(int v=0; v<SD; v++)
		  temp_vect2[v] = PackTetrah[n].Vertex[0][v]; 
	      }
	    else if(j==3)
	      {
        	for(int v=0; v<SD; v++)
		  temp_vect2[v] = PackTetrah[n].Vertex[1][v]; 
	      }
            else if(j>3)
	      {
        	for(int v=0; v<SD; v++)
		  temp_vect2[v] = PackTetrah[n].Vertex[3][v]; 
	      }


	    for(int v=0; v<SD; v++)
	      temp_vect1[v] = temp_vect1[v] - temp_vect2[v];

	    GetPLine(edgeM[i], temp_vect1, temp_axis);
	  }


	//the separating axis has been obtained, now check the vertices...

	double temp_disM[4], temp_disN[4];

	for(int k=0; k<4; k++)
	  {
	    temp_disM[k] = 0;
	    temp_disN[k] = 0;
	  }

	for(int u=0; u<4; u++)
	  {
	    for(int k=0; k<SD; k++)
	      {
		temp_disM[u] += temp_axis[k]*(PackTetrah[m].Vertex[u][k] + temp_VdisE[k]);
		temp_disN[u] += temp_axis[k]*PackTetrah[n].Vertex[u][k];
	      }
	  }

	if(GetMax(temp_disM,4)<GetMin(temp_disN,4) || GetMax(temp_disN,4)<GetMin(temp_disM,4))
	  {
	    overlap_flag = 0;

	    return overlap_flag;
	  }
      }//end of the loop for the edge pairs...


  //if we get here, then the M and N overlap...

  return overlap_flag;

}


//work with the relative coordinates,
void GetTranslation(int m)
{
  //srand(time(NULL));

  double p[SD], dL[SD];

  for(int k=0; k<SD; k++)
    {
      p[k] = (double)(rand()%MAXX)/(double)MAXX;
      dL[k] = (p[k]-0.5)*TransMod;
    }

  for(int k=0; k<SD; k++)
    {
      Old_CenterL[k] = CenterL[m][k];

      CenterL[m][k] = CenterL[m][k] + dL[k];
      //update the relative position and make sure the particle is in the unit box...
      if(CenterL[m][k]>=1.0) CenterL[m][k] = CenterL[m][k]-1.0;
      else if(CenterL[m][k]<0) CenterL[m][k] = CenterL[m][k]+1.0;
    }

  //Note the cell which the particle belongs to may have changed, we will check it in MCmove()

}


void GetRotation(int m)
{
  //srand(time(NULL));

  //store the old orenation info. in Old_Tetrah...
  for(int i=0; i<4; i++)
    {
      for(int k=0; k<SD; k++)
	{
	  Old_Tetrah.Vertex[i][k] = PackTetrah[m].Vertex[i][k];
	  Old_Tetrah.Normal[i][k] = PackTetrah[m].Normal[i][k];
	}
    }


  //now we compute the rotations for each vetex of the tetrah...
   double e1[SD], e2[SD], e3[SD];

   double pr = (double)(rand()%MAXX)/(double)MAXX;

   double Phi = (pr-0.5)*RotMod; //the rotation angle...

   //the rotation axis...
   for(int k=0; k<SD; k++)
     e1[k] = (double)(rand()%MAXX)/(double)MAXX - 0.5;

   double temp_sum = sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);

   for(int k=0; k<SD; k++)
     e1[k] = e1[k]/temp_sum;
   
   //printf("e1 = (%f, %f, %f)\n", e1[0], e1[1], e1[2]);


   //for the vertices of the Tetrah[m]
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   for(int i=0; i<4; i++)
     {
       double rx = PackTetrah[m].Vertex[i][0];
       double ry = PackTetrah[m].Vertex[i][1];
       double rz = PackTetrah[m].Vertex[i][2];


       //printf("Vertex[%d] = (%f, %f, %f)\n", i, PackTetrah[m].Vertex[i][0], PackTetrah[m].Vertex[i][1], PackTetrah[m].Vertex[i][2]);

       double R1 = rx*e1[0] + ry*e1[1] + rz*e1[2];

       //get the perpendicular components and e2
       e2[0] = rx-R1*e1[0]; e2[1] = ry-R1*e1[1]; e2[2] = rz-R1*e1[2];

       double tempd = sqrt(e2[0]*e2[0]+e2[1]*e2[1]+e2[2]*e2[2]);

       e2[0] = e2[0]/tempd; e2[1] = e2[1]/tempd; e2[2] = e2[2]/tempd;

       //printf("e2[%d] = (%f, %f, %f)\n", i, e2[0], e2[1], e2[2]);

       //get e3 as follows or GetCrossProduct(e1, e2, e3)
       e3[0] = e1[1]*e2[2]-e1[2]*e2[1];
       e3[1] = e1[2]*e2[0]-e1[0]*e2[2];
       e3[2] = e1[0]*e2[1]-e1[1]*e2[0];

       //printf("e3[%d] = (%f, %f, %f)\n", i, e3[0], e3[1], e3[2]);

       //get the components..
       double R2 = tempd*cos(Phi);
       double R3 = tempd*sin(Phi);

       //the new posistion of the Vertex
       PackTetrah[m].Vertex[i][0] = R1*e1[0] + R2*e2[0] + R3*e3[0];
       PackTetrah[m].Vertex[i][1] = R1*e1[1] + R2*e2[1] + R3*e3[1];
       PackTetrah[m].Vertex[i][2] = R1*e1[2] + R2*e2[2] + R3*e3[2];

       //printf("Vertex[%d] = (%f, %f, %f)\n", i, PackTetrah[m].Vertex[i][0], PackTetrah[m].Vertex[i][1], PackTetrah[m].Vertex[i][2]);
       //printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

     } //end of loop over all the 4 vertices of tetrah[m]

   /*
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   //This is a check....
   double OT[SD];
   for(int k=0; k<SD; k++)
     OT[k] =0;
   for(int i=0; i<4; i++)
     for(int k=0; k<SD; k++)
       OT[k] += PackTetrah[m].Vertex[i][k];

   printf("(O.x, O.y, O.z) of [%d] = %f\t%f\t%f\n", m, OT[0], OT[1], OT[2]);
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */

   //IMPORTANT: Updat the Normal!!!!!
   GetNormal(m);

   //printf("*****************************************************\n");

}

void RetainTranslation(int m)
{
  for(int k=0; k<SD; k++)
    {
      CenterL[m][k] =  Old_CenterL[k];  
    }
}

void RetainRotation(int m)
{
  for(int i=0; i<4; i++)
    {
      for(int k=0; k<SD; k++)
	{
	   PackTetrah[m].Vertex[i][k] = Old_Tetrah.Vertex[i][k];
	   PackTetrah[m].Normal[i][k] = Old_Tetrah.Normal[i][k];
	}
    }
}


void PrintCellList()
{
   for(int i=0; i<Ncell; i++)
     {
       node* pt = CellList[i];

       int temp_indexK = floor((double)i/(double)(LXcell*LYcell));
       int temp_indexJ = floor((double)(i%(LXcell*LYcell))/(double)LXcell);
       int temp_indexI = (i%(LXcell*LYcell))%LXcell;

       printf("CellList(%d,%d,%d) = ", temp_indexI, temp_indexJ, temp_indexK);
       
       while(pt!=NULL)
	 {
	   printf(" %d ", pt->index_Tetrah);
	   pt = pt->next;
	 }
       
       printf("\n");
       
     }
   
}

void UpdateCellList(int m)
{
  //delete the particle from the old list...
  int temp_indI = floor(Old_CenterL[0]*LXcell);
  int temp_indJ = floor(Old_CenterL[1]*LYcell);
  int temp_indK = floor(Old_CenterL[2]*LZcell);

  int cell_index = LXcell*LYcell*temp_indK + LXcell*temp_indJ + temp_indI;

  node* pt1 = CellList[cell_index];
  node* pt2 = CellList[cell_index];

  if(pt1 == NULL)
    { 
       printf("Bugs exist in CellList! Re-check!\n");
       //print out the list..
       printf("The particle is %d\n", m);
       printf("It should be in (%d, %d, %d)\n", temp_indI, temp_indJ, temp_indK);
       PrintCellList();
    
       exit(1);
    }
  else
    {
      while(pt1!=NULL)
	{
	  if((pt1->index_Tetrah)==m)
	    {
	      if(pt2==CellList[cell_index] && pt1==CellList[cell_index])
		CellList[cell_index] = pt1->next;
	      else
                pt2->next = pt1->next;
	      free(pt1);
	      break;
	    }

	  pt2 = pt1;
	  pt1 = pt1->next;

	}
      /*
      if(pt1==NULL)
	{
	  printf("Bugs exist in celllist! Recheck!\n");
	  exit(1);
	}
      */
    }


  //insert the particle in the new list...
  temp_indI = floor(CenterL[m][0]*LXcell);
  temp_indJ = floor(CenterL[m][1]*LYcell);
  temp_indK = floor(CenterL[m][2]*LZcell);

  if(temp_indI>=LXcell || temp_indJ>=LYcell || temp_indK>=LZcell)
	{
	  printf("CenterL includes invalid/un-scaled coordinates! Recheck!\n");
	  exit(1);
	}

  int new_cell_index = LXcell*LYcell*temp_indK + LXcell*temp_indJ + temp_indI;

  node* pt = (node *)malloc(sizeof(node));
  pt->index_Tetrah = m;
  pt->next = CellList[new_cell_index];

  CellList[new_cell_index] = pt;

}



//each move could be a translation or a rotation with equal probability
//return 1 if the trial traslation is successful
//retrun 2 if the trial rotation is sucessful
//return 0 if fails....
int MCmove(int m)
{

  //srand(time(NULL));

  double pm = (double)(rand()%MAXX)/(double)MAXX;

  if(pm<=p_trans) //translational trail move
    {
      //translate the particle
      GetTranslation(m);

      //printf("Translation is performed!\n");

      //need to check overlapping...
      int temp_indI = floor(CenterL[m][0]*LXcell);
      int temp_indJ = floor(CenterL[m][1]*LYcell);
      int temp_indK = floor(CenterL[m][2]*LZcell);

     
     
      for(int i=-1; i<=1; i++)
	for(int j=-1; j<=1; j++)
	  for(int k=-1; k<=1; k++)
	    {
	      //for the box index of temp_P, they are the same for each cell
	      int indexI = 0;
	      int indexJ = 0; 
	      int indexK = 0;

	      int tempI = temp_indI + i;
	      if(tempI>=LXcell) 
		{
		  tempI = tempI-LXcell;
		  indexI = 1;
		}
	      else if(tempI<0) 
		{
		  tempI = tempI+LXcell;
		  indexI = -1;
		}
	      int tempJ = temp_indJ + j;
	      if(tempJ>=LYcell) 
		{
		  tempJ = tempJ-LYcell;
		  indexJ = 1;
		}
	      else if(tempJ<0) 
		{
		  tempJ = tempJ+LYcell; 
		  indexJ = -1;
		}
	      int tempK = temp_indK + k;
	      if(tempK>=LZcell) 
		{ 
		  tempK = tempK-LZcell; 
		  indexK = 1;
		}
	      else if(tempK<0) 
		{ 
		  tempK = tempK+LZcell; 
		  indexK = -1;
		}

	      int cell_index = LXcell*LYcell*tempK + LXcell*tempJ + tempI;

	      node* pt = CellList[cell_index];

	      while(pt!=NULL)
		{
		  int temp_P = pt->index_Tetrah;

		  if(CheckOverlap(temp_P, m, indexI, indexJ, indexK)==1) //the particles overlap
		    {
		      RetainTranslation(m);

		      return 0; //indicating an unsuccessful trial...
		    }

		  pt = pt->next;
		}
	    }//end of loop over all neighboring cells...

      //if we reach here, then means no overlap...
      //only need to check and update the celllist...
         
      int temp2_indI = floor(Old_CenterL[0]*LXcell);
      int temp2_indJ = floor(Old_CenterL[1]*LYcell);
      int temp2_indK = floor(Old_CenterL[2]*LZcell);

      if(temp2_indI==temp_indI && temp2_indJ==temp_indJ && temp2_indK==temp_indK)
	{
	  //no need to update the CellList
	  return 1;
	}
      else
	{
	  UpdateCellList(m);
	  return 1;
	}
        //indicating a successful trail move....

    } //ENDIF
  else //rotational trail move...
    {
      //make a rotation...
      GetRotation(m);

      //printf("Rotation is performed!\n");

      //need to check overlapping...
      int temp_indI = floor(CenterL[m][0]*LXcell);
      int temp_indJ = floor(CenterL[m][1]*LYcell);
      int temp_indK = floor(CenterL[m][2]*LZcell);

 
      for(int i=-1; i<=1; i++)
	for(int j=-1; j<=1; j++)
	  for(int k=-1; k<=1; k++)
	    {
	      int indexI = 0;
	      int indexJ = 0;
	      int indexK = 0;
	     
	      int tempI = temp_indI + i;
	      if(tempI>=LXcell) 
		{
		  tempI = tempI-LXcell;
		  indexI = 1;
		}
	      else if(tempI<0) 
		{
		  tempI = tempI+LXcell;
		  indexI = -1;
		}
	      int tempJ = temp_indJ + j;
	      if(tempJ>=LYcell) 
		{
		  tempJ = tempJ-LYcell;
		  indexJ = 1;
		}
	      else if(tempJ<0) 
		{
		  tempJ = tempJ+LYcell; 
		  indexJ = -1;
		}
	      int tempK = temp_indK + k;
	      if(tempK>=LZcell) 
		{ 
		  tempK = tempK-LZcell; 
		  indexK = 1;
		}
	      else if(tempK<0) 
		{ 
		  tempK = tempK+LZcell; 
		  indexK = -1;
		}
	      
	      int cell_index = LXcell*LYcell*tempK + LXcell*tempJ + tempI;

	      //node* pt = CellList[tempI][tempJ][tempK];
	      node* pt = CellList[cell_index];

	      

	      /*
	      if(err_flag>2000)
		printf("(I, J, K) = (%d, %d, %d)\n", tempI, tempJ, tempK);
	      */

	      while(pt!=NULL)
		{
		  int temp_P = pt->index_Tetrah;

		  /*
		  if(err_flag>2000 && tempI==4 && tempJ==4 && tempK==3)
		    printf("Temp_P = %d\n", temp_P);
		  */

		  if(CheckOverlap(temp_P, m, indexI, indexJ, indexK)==1) //the particles overlap
		    {
		      RetainRotation(m);

		      return 0; //indicating an unsuccessful trial...
		    }

		  pt = pt->next;
		}
	    }//end of loop over all neighboring cells...


      return 2; //indicating a successful rotation trail 


    }//ENDELSE

}

void PrintOverlapPair(int, int);

int GlobalOverlapCheck(int a)
{
  int overlap_flag = 0;

  for(int m=0; m<N; m++)
    for(int n=m; n<N; n++)
      {
	
	for(int i=-1; i<=1; i++)
	  for(int j=-1; j<=1; j++)
	    for(int k=-1; k<=1; k++)
	      {
		
		if(CheckOverlap(n, m, i, j, k) == 1) 
		  {
		    printf("Tetrah %d <%d, %d, %d> and Tetrah %d overlap!\n", n, i, j, k, m);
		    //PrintOverlapPair(n, m);
		    overlap_flag = 1;
		  }
	      }
      }


  return overlap_flag;
}

int GlobalOverlapCheck()
{
  int overlap_flag = 0; //assuming no overlapping...

  for(int m=0; m<N; m++)
    {
      
      int temp_indI = floor(CenterL[m][0]*LXcell);
      int temp_indJ = floor(CenterL[m][1]*LYcell);
      int temp_indK = floor(CenterL[m][2]*LZcell);

     
      //cout<<" Neighbor of ["<<m<<"] = ( ";

      for(int i=-1; i<=1; i++)
	for(int j=-1; j<=1; j++)
	  for(int k=-1; k<=1; k++)
	    {
	      int indexI = 0;
	      int indexJ = 0;
	      int indexK = 0;
	      
	      int tempI = temp_indI + i;
	      if(tempI>=LXcell) 
		{
		  tempI = tempI-LXcell;
		  indexI = 1;
		}
	      else if(tempI<0) 
		{
		  tempI = tempI+LXcell;
		  indexI = -1;
		}
	      int tempJ = temp_indJ + j;
	      if(tempJ>=LYcell) 
		{
		  tempJ = tempJ-LYcell;
		  indexJ = 1;
		}
	      else if(tempJ<0) 
		{
		  tempJ = tempJ+LYcell; 
		  indexJ = -1;
		}
	      int tempK = temp_indK + k;
	      if(tempK>=LZcell) 
		{ 
		  tempK = tempK-LZcell; 
		  indexK = 1;
		}
	      else if(tempK<0) 
		{ 
		  tempK = tempK+LZcell; 
		  indexK = -1;
		}

	      

	      int cell_index = LXcell*LYcell*tempK + LXcell*tempJ + tempI;
	      
	      
	      //node* pt = CellList[tempI][tempJ][tempK];
	      node* pt = CellList[cell_index];

	    
              
	      while(pt!=NULL)
		{
		  int temp_P = pt->index_Tetrah;
		  
		  //cout<<temp_P<<"<"<<indexI<<","<<indexJ<<","<<indexK<<">"<<"\t";

		  if(CheckOverlap(temp_P, m, indexI, indexJ, indexK)==1) //the particles overlap
		    {
		      overlap_flag = 1;

		      printf("Tetrah %d and Tetrah %d overlap!\n", temp_P, m);
                      //PrintOverlapPair(temp_P, m);		    

		      return overlap_flag; //indicating overlapping...
		    }

		  pt = pt->next;
		}
	    }//end of loop over all neighboring cells...

      
      //cout<<")"<<endl;

    }//end of loop over all particles...

 
  
  //if we reach here, then no overlap...

  return overlap_flag;

}

void BoundaryDeform()
{
  //store the old lambda first...
  for(int i=0; i<SD; i++)
    for(int j=0; j<SD; j++)
      Old_Lambda[i][j] = Lambda[i][j];


  //chack the boundary...
  double D_Lambda[SD][SD];

  for(int i=0; i<SD; i++)
    for(int j=0; j<SD; j++)
      D_Lambda[i][j]  =0;

  for(int i=0; i<SD; i++)
    for(int j=0; j<SD; j++)
      {
	for(int k=0; k<SD; k++)
	  D_Lambda[i][j] += Strain[i][k]*Lambda[k][j];
      }

  for(int i=0; i<SD; i++)
    for(int j=0; j<SD; j++)
      Lambda[i][j] += D_Lambda[i][j];
  
}


void BoundaryRetain()
{
  //store the old lambda first...
  for(int i=0; i<SD; i++)
    for(int j=0; j<SD; j++)
      Lambda[i][j] = Old_Lambda[i][j];
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//*****************************************************************
//subroutines for post-processing

void GetVolume()
{
   V_sys =  Lambda[0][2]*(Lambda[1][0]*Lambda[2][1]-Lambda[1][1]*Lambda[2][0]) + Lambda[1][2]*(Lambda[2][0]*Lambda[0][1] - Lambda[0][0]*Lambda[2][1]) + Lambda[2][2]*(Lambda[0][0]*Lambda[1][1] - Lambda[1][0]*Lambda[0][1]);
}

double GetTetrahVol(int m)
{
  double vect1[SD], vect2[SD], vect3[SD], vect4[SD];

  for(int i=0; i<SD; i++)
    {
      vect1[i] = PackTetrah[m].Vertex[0][i] - PackTetrah[m].Vertex[3][i];
      vect2[i] = PackTetrah[m].Vertex[1][i] - PackTetrah[m].Vertex[3][i];
      vect3[i] = PackTetrah[m].Vertex[2][i] - PackTetrah[m].Vertex[3][i];
    }

  GetCrossProduct(vect2, vect3, vect4);

  double vol = fabs(GetInnerProduct(vect1, vect4))/6.0;

  return vol;
}

double GetDensity()
{
  double VLambda;
  VLambda = Lambda[0][2]*(Lambda[1][0]*Lambda[2][1]-Lambda[1][1]*Lambda[2][0]) + Lambda[1][2]*(Lambda[2][0]*Lambda[0][1] - Lambda[0][0]*Lambda[2][1]) + Lambda[2][2]*(Lambda[0][0]*Lambda[1][1] - Lambda[1][0]*Lambda[0][1]);
 
  double VTetrah = 0.0;

  for(int i=0; i<N; i++)
    VTetrah += GetTetrahVol(i);

  //VTetrah = (d0*d0*d0)/sqrt(72.0); //for a regular tetrah

  return VTetrah/VLambda;

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//collect pressure...

void Boundary_Shrink(double Temp_GapMod)
{
  //change the boundary...
  double D_Lambda[SD][SD];

  //get the temp strain..
  double Temp_Strain[SD][SD];
  
   //store the old lambda first...
  for(int i=0; i<SD; i++)
    for(int j=0; j<SD; j++)
      {
	Old_Lambda[i][j] = Lambda[i][j];
	D_Lambda[i][j] = 0;
	Temp_Strain[i][j] = 0;
      }
 
  for(int i=0; i<SD; i++)
    {
      Temp_Strain[i][i] = -Temp_GapMod;
      //printf("Strain[%d][%d] = %f\n", i, i, Temp_Strain[i][i]);
    }

  //get the new vector...
  for(int i=0; i<SD; i++)
    for(int j=0; j<SD; j++)
      {
	for(int k=0; k<SD; k++)
	  D_Lambda[i][j] += Temp_Strain[i][k]*Lambda[k][j];
      }

  for(int i=0; i<SD; i++)
    for(int j=0; j<SD; j++)
      Lambda[i][j] += D_Lambda[i][j]; 
}

void GetPressure()
{
  //**************************************************
  //perform the virtual volume change...
  Boundary_Shrink(GapMod);

  int overlap_flag = GlobalOverlapCheck();
  Stage_Pressure += (1-overlap_flag); //computing the acceptance rate...
  
  BoundaryRetain(); //these are just test deform, boundaried should always be retained...
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//compute g2

double GetNDensity() // the number density, for computing g2...
{
  double VLambda;
  VLambda = Lambda[0][2]*(Lambda[1][0]*Lambda[2][1]-Lambda[1][1]*Lambda[2][0]) + Lambda[1][2]*(Lambda[2][0]*Lambda[0][1] - Lambda[0][0]*Lambda[2][1]) + Lambda[2][2]*(Lambda[0][0]*Lambda[1][1] - Lambda[1][0]*Lambda[0][1]);

  //double VTetrah;

  //VTetrah = (d0*d0*d0)/sqrt(72.0);

  return (double)N/VLambda;
}


void Get_SampleLength() //only need to be done once for each stage...
{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //get the required the parameters... one tenth of the insphere radius...
  //double N_bin = 10.0;

  bin = d0/(N_bin*sqrt(24.0)); //sqrt(24.0) is for the scaling purpose...

  //compute the minimum length of the lattice vectors
  double temp_d;
  double temp_dmin = 100000.0;

  for(int i=0; i<SD; i++)
    {
      temp_d = 0;

      for(int j=0; j<SD; j++)
	temp_d += Lambda[i][j]*Lambda[i][j];

      if(temp_d < temp_dmin) temp_dmin = temp_d;
    }

  temp_dmin = sqrt(temp_dmin);

  //obtain the sample length....
  NLcounter = floor(2.0*temp_dmin/(2.0*bin));

  if(NLcounter>MaxNL) NLcounter = MaxNL;

}

void Get_g2(double rho_n)
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  for(int i=0; i<MaxNL; i++)
     {
        g[i] = 0;
     }

  //printf("Computing G2 and g2 now....\n");

  //loop over all particles
  for(int i=0; i<N; i++)
    {
      for(int j=0; j<N; j++)
	{

	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          //Note we didn't truncate the distances that are greater than half the box length...
	  //for the radiual g2...
	  double temp_Vdis[SD];
	  double temp_dis = MinDis(j, i, temp_Vdis)/bin;

	  if(temp_dis>0 && temp_dis<NLcounter) //this means the particle itself is excluded...
	    {
	      int int_dis = floor(temp_dis);
	      double delta_dis = temp_dis - int_dis;
	      if(delta_dis>0.5) int_dis++;

	      g[int_dis]++;
	    }

	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	}//end of the inner loop
    }//end of the outer loop

 
  //now we normalize g2 and print the results....
  //rho_n = GetNDensity();
  //double pi = 3.1415926;

  for(int r=0; r<NLcounter; r++)
    {
      g[r] = g[r]/(N*rho_n*4*pi*(r+0.5)*(r+0.5)*bin*bin*bin);
    }

}



//**************************************************************
//

void PrintOverlapPair(int m, int n)
{
  FILE* fp = fopen("OverlapPair.txt","w");

  fprintf(fp, "%d\n", 2);
  fprintf(fp, "%.12f\n", d0);

  for(int i=0; i<SD; i++)
    {
    for(int j=0; j<SD; j++)
      fprintf(fp, "%.12f\t", Lambda[j][i]);

    fprintf(fp, "\n");
    }
  fprintf(fp, "\n");


 
      fprintf(fp, "%.12f\t%.12f\t%.12f\n", CenterL[m][0], CenterL[m][1], CenterL[m][2]);

      for(int j=0; j<4; j++)
          fprintf(fp, "%.12f\t%.12f\t%.12f\n", PackTetrah[m].Vertex[j][0], PackTetrah[m].Vertex[j][1], PackTetrah[m].Vertex[j][2]);
     
      fprintf(fp, "\n");

 fprintf(fp, "%.12f\t%.12f\t%.12f\n", CenterL[n][0], CenterL[n][1], CenterL[n][2]);

      for(int j=0; j<4; j++)
          fprintf(fp, "%.12f\t%.12f\t%.12f\n", PackTetrah[n].Vertex[j][0], PackTetrah[n].Vertex[j][1], PackTetrah[n].Vertex[j][2]);
     
      fprintf(fp, "\n");
   

  fclose(fp);
}

void PrintTetrah()
{
  FILE* fp = fopen("Fconfig.txt","w");

  fprintf(fp, "%d\n", N);
  fprintf(fp, "%.12f\n", d0);

  for(int i=0; i<SD; i++)
    {
    for(int j=0; j<SD; j++)
      fprintf(fp, "%.12f\t", Lambda[j][i]);

    fprintf(fp, "\n");
    }
  fprintf(fp, "\n");


  for(int i=0; i<N; i++)
    {
      fprintf(fp, "%.12f\t%.12f\t%.12f\n", CenterL[i][0], CenterL[i][1], CenterL[i][2]);

      for(int j=0; j<4; j++)
          fprintf(fp, "%.12f\t%.12f\t%.12f\n", PackTetrah[i].Vertex[j][0], PackTetrah[i].Vertex[j][1], PackTetrah[i].Vertex[j][2]);
     
      fprintf(fp, "\n");
    }

  fclose(fp);
}




//*****************************************************************
//checking the subroutines...
void FunctionCheck()
{
   Read_Config();

   //GetGlobalPosition();

  GetCellList(0);

  /*
  //print the normal ...
  for(int i=0; i<N; i++)
    {
      printf("Normal_Tetrah[%d] = { ", i);

      for(int j=0; j<4; j++)
	{
	  printf("( ");
	  for(int k=0; k<SD; k++)
	    printf("%f ", PackTetrah[i].Normal[j][k]);
	  printf("), ");
	}

      printf("}\n");
    }
  */

  //print the cell list...
  PrintCellList();
  
  /*
  //print check results
  for(int m=0; m<N; m++)
    for(int n=0; n<N; n++)
      {
	int res = CheckOverlap(m,n);

	printf("CheckOverlap(%d, %d) = %d\n", m, n, res);
      }


  //check GetPLine(Direct, Vect, PLine)
  double Temp_PLine[SD];
  double Temp_Direct[SD];
  Temp_Direct[0] = Temp_Direct[1] = Temp_Direct[2] = 1.0;
  double Temp_Vect[SD];
  Temp_Vect[0] = Temp_Vect[1] = 1.0; Temp_Vect[2] = -1.0;

  GetPLine(Temp_Direct, Temp_Vect, Temp_PLine);

  printf("PLine = (%f, %f, %f)\n", Temp_PLine[0], Temp_PLine[1], Temp_PLine[2]);
  

  //GetRotation(0);

  //GetRotation(1);


  printf("The MC moves start:\n");
  for(int i=0; i<100; i++)
    {
      int m = rand()%N;
      int succ = MCmove(m);
      printf("%d", succ);
   
    }
  printf("\n");
  */

}


//***********************************************************************
main()
{
  //printf("check!\n");

  srand(time(NULL));

  //FunctionCheck();

  FILE* fp;

  fp = fopen("rho.txt","w");
  fclose(fp);
  
  if(flag_comput_pressure == 1)
    {
      fp = fopen("pressure.txt","w");
      fclose(fp);
    }
  
  if(flag_comput_g2 == 1)
    {
      fp = fopen("g2.txt","w");
      fclose(fp);
    }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //pre-processing...
  
  Int_Data();

  Read_Config();

  //PrintTetrah();

  //GetStrain();

  GetCellList(0);
  rho_cell = GetDensity();

  PrintCellList();

  //check the initial configuration...
  int temp_option = -1;
  cout<<" Rescale the packing or Exit? Rescale - 0; Exit - 1"<<endl;
  cin >> temp_option; //this also read in from the input.txt file
  
  int temp_flag = GlobalOverlapCheck(1);
  if(temp_flag == 1)
    {
  
      if(temp_option == 0)
	{
	  while(GlobalOverlapCheck(1)==1)
	    {
	      printf("The Initial Configuration Contains Overlapping Pairs!\n");
	      Rescale();
	      //exit(1);
	    }
	}
      else if(temp_option == 1)
	{ 
	  exit(1);
	}
    }

 
  
  //rescale the packing if it is denser than the specified starting density
  double rho_int = GetDensity();
  if(rho_int<Starting_Density)
    {
      printf("The initial packing is less dense than rho_start = %f\n", Starting_Density);
      printf("Start the simulation at rho_int = %f\n", rho_int);
    }
  else{
    while(rho_int>=Starting_Density)
      {
	Rescale();
	rho_int = GetDensity();
      }
    printf("The specified starting density is rho_start = %f\n", Starting_Density);
    printf("Start the simulation at rho_int = %f\n", rho_int);
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //perform the moves..

  printf("MC moves start...\n");

  for(int i=0; i<Nstage; i++) //the stages...
    {
      printf("The %dth stage starts...\n", i+1);

      //initilize the property statistics

      //for the pressure
      if(flag_comput_pressure == 1)
	{
	  
	  rho_s = GetDensity();
	  GetVolume();
	  Stage_Pressure = 0; //initialize ...
	  pressure_ct = 0;
	}
      
      
      //for the g2
      if(flag_comput_g2 == 1)
	{
	  rho_n = GetNDensity();
	  Get_SampleLength();
	  g2_ct = 0;
	  
	  for(int m=0; m<MaxNL; m++)
	    g_ave[m] = 0;
	}
      

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      int sum_res = 0;

      for(int j=0; j<Ncycle*N; j++) //the cycles...
	{
	  //int index = rand()%N; //this is randomly picking up th particles

	  int index = j%N; //this is a linear sequential picking up...

	  int res = MCmove(index);
	  
	  if(res>0)
	    sum_res++;

	  //printf("%d", res);
	  
	  //collect statistics after system has been equilibrated for a while
	  if(j>floor((double)Ncycle*N/2.0) && j%(Npc*N)==0)
	    {
	      //collect pressure...
	      if(flag_comput_pressure == 1)
		{
		  GetPressure();
		  pressure_ct++;
		}
	      
	      //compute g_2
	      if(flag_comput_g2 == 1)
		{
		  Get_g2(rho_n);
		  for(int m=0; m<MaxNL; m++)
		    g_ave[m] += g[m];     
		  g2_ct++;
		}
	      
	    }
	}

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //compute the acc_rate and rescale the step size when necessary
      double R_Acc = (double)sum_res/(double)(Ncycle*N);

      printf("Acceptance Rate: Acc = %.10f\n", R_Acc);

      //permannetly adjust the step sizes if necessary...
      int Hnew_Acc = log(1.0/R_Acc);
      printf("Hnew_Acc = %d\n", Hnew_Acc);
      //for adaptive step size, in both directions...

      if(Hnew_Acc > 1) //the step size is too large...
	{
	  TransMod = TransMod*Trans_Rescale;
	  RotMod = RotMod*Rot_Rescale;

	  Global_StrainMod = Global_StrainMod*Strain_Rescale; //in the meantime, scale the strain magnitude
	  Shear_StrainMod = Shear_StrainMod*Strain_Rescale;
	  //thus, it is better to keep the StrainMod consistent with the trial move step size...
	  cout<<"The TransMod, RotMod and StrainMod are rescaled!"<<endl;
	}
   

      printf("\nThe %dth stage ends...\n", i+1);

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //Now we perform a GlobalOverlapCheck()...
      //int overlap_flag = GlobalOverlapCheck();

      /*
      if(GlobalOverlapCheck() == 1)
	{
	  printf("The current configuration contains overlapping particles! Recheck!\n");
	  exit(1);
	}
      */


      //Now we deform the boundary and report density...
      GetStrainIso(Global_StrainMod, Shear_StrainMod); //each time we use the orginally specified strain rate

      BoundaryDeform();

      int overlap_flag = GlobalOverlapCheck();
      int limit = 0;

      while(overlap_flag==1 && limit<Ndeform)
	{
	  BoundaryRetain(); //retain the old boundary...

	  ResetStrain(limit);

	  //GetStrainIso();

	  BoundaryDeform();

	  overlap_flag = GlobalOverlapCheck();

	  limit++;
	}

      if(limit==Ndeform) 
        {
          BoundaryRetain();
	  printf("Failed to compress at stage %d\n", i+1);
	}

      /*
      //we put this double check here....
      overlap_flag = GlobalOverlapCheck();

      if(overlap_flag == 1)
	{
	  printf("The current configuration contains overlapping particles! Recheck!\n");
	  exit(1);
	}
      */
      //if we are unlucky, i.e., limit>Ndeform, we do not compress at the current stage...

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //Now compute and report the density and other statistics...

      double rho = GetDensity();

      printf("The density of %dth stage: rho = %f\n", (i+1), rho);


      //the density change is larger than 5%, re-set the cell-list
      //cout<<"(rho, rho_cell) = "<<rho<<", "<<rho_cell<<endl;
      if(fabs(rho-rho_cell)>0.05)
	{
	  cout<<"(rho, rho_cell) = "<<rho<<", "<<rho_cell<<endl;
	  GetCellList(1);
	  rho_cell = rho;

	  cout<<"The CellList has been rebulit!"<<endl;
	}


      //print the density vs stage to file
      fp = fopen("rho.txt","a");
      fprintf(fp,"%d\t%f\n", i+1, rho);
      fclose(fp);

      if(flag_comput_pressure == 1)
	{
	  //Now compute and print the stage pressure...
	  Stage_Pressure = -V_sys*log(Stage_Pressure/(double)pressure_ct)/(N*3*GapMod);
	  printf("The pressure of %dth stage: pressure = %f\n", (i+1), Stage_Pressure);
	  printf("pressure_ct = %d\n", pressure_ct);

	  //print the density vs pressure to file
	  fp = fopen("pressure.txt","a");
	  fprintf(fp,"%f\t%f\n", rho_s, Stage_Pressure);
	  fclose(fp);
	}

      if(flag_comput_g2 == 1)
	{
	  //print the g2 of this stage to file...
	  fp = fopen("g2.txt","a");
	  fprintf(fp, "#~~~~~~~~~~~~ Stage# %d ~~~~~~~~~~~~\n", (i+1));
	  for(int r=0; r<NLcounter; r++)
	    fprintf(fp, "%f\t%f\n", (r)/(2*N_bin), g_ave[r]/g2_ct);
	  fprintf(fp, "#**************************************\n");
	  fclose(fp);
	}


      printf("**********************************\n");

      if(i%50==0) PrintTetrah(); //print out the packing every 50 stages... 

      if(rho>=Terminate_Density)
	{
	  PrintTetrah();
	  printf("The density has reached rho = %f\n!", Terminate_Density);
	  printf("The program is terminated as requested!\n");
	  exit(1);
	}

      //PrintCellList();

    }



  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //post-processing...

  if(GlobalOverlapCheck(1)==1)
    {
      printf("The Final Configuration Contains Overlapping Pairs!\n");
      exit(1);
    }


  PrintTetrah();

 
}
