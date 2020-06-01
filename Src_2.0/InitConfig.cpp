//produce relative dense initial configurations, from MRJ packings of spheres
//      this is to maximize translational disorder and see how orientation can be correlated
//such configurations can be rescaled to a low-density and system will quickly loose memory




//author: Yang JIAO, yjiao@princeton.edu
//started: Sep. 03, 2010

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using namespace std;

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "const.h"


#include "Geometry.h"
#include "Polyhedron.h"
#include "GEsolver.h"

ifstream fin;
ofstream fout;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//double pi = 3.141592654;


int N; //the number of particles
double d0; //the diameter of the sphere
double d_cv_max = 0.0; //the maximal centroid to vertice distance
double d_edge; //the edge length of the shape

double relax_ratio = 1.01; //rescale the packing if the current density is larger than required
double rho_tar; //the target density
double rescale_ratio; //to rescale the sphere packing, i.e. D = 2*circums-radius of the particle

double VParticle = 0.0;

double Lambda[SD][SD];

double** CenterE;
double** CenterL;

Polyhedron* PackPoly; 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void Read_SpherePacking()
{
  fin.open("MRJ.txt");

  if(!fin)
    {
      cout<<"Cannot open file <MRJ.txt> to read! Check again!"<<endl;
      exit(1);
    }
  else
    {
      cout<<"Initializing polyhedron packing from MRJ sphere packing..."<<endl;
    }

  char dummy[100];
  int temp_int;
  double temp_val;

  fin>>temp_int>>dummy>>dummy;
  fin>>N>>temp_int;
  fin>>temp_int;
  fin>>d0;

  CenterE = new double*[N];
  CenterL = new double*[N];
  for(int i=0; i<N; i++)
    {
      CenterE[i] = new double[SD];
      CenterL[i] = new double[SD];
    }

  for(int i=0; i<SD; i++) // read in the lattice vectors...
    for(int j=0; j<SD; j++)
      {
	fin>>Lambda[j][i];
      }

  fin>>dummy>>dummy>>dummy;

  for(int i=0; i<N; i++)
    fin>>CenterE[i][0]>>CenterE[i][1]>>CenterE[i][2];

  fin.close();

}


//compute the local coordinates
void Get_CenterL()
{
  //obtain and normalize the local coordinates..
  for(int n=0; n<N; n++)
    {
      GESOL solver;
      solver.GEsolver(SD, SD, Lambda, CenterE[n]);
      
      for(int i=0; i<SD; i++)
	{
	  CenterL[n][i] = solver.back[i];
	  
	  if(CenterL[n][i]>=1.0) CenterL[n][i] = CenterL[n][i] - 1.0;
	  else if(CenterL[n][i]<0) CenterL[n][i] = CenterL[n][i] + 1.0;
	}
    }
}


void Get_VParticle()
{
  VParticle = 0.0;

  for(int i=0; i<N; i++)
    {      
      VParticle += PackPoly[i].GetVol();  
      //cout<<"VParticle = "<<VParticle<<endl;
    }
}

double Get_Density()
{
  double VLambda;
  VLambda = Lambda[0][2]*(Lambda[1][0]*Lambda[2][1]-Lambda[1][1]*Lambda[2][0]) + Lambda[1][2]*(Lambda[2][0]*Lambda[0][1] - Lambda[0][0]*Lambda[2][1]) + Lambda[2][2]*(Lambda[0][0]*Lambda[1][1] - Lambda[1][0]*Lambda[0][1]);

  return VParticle/VLambda;

}


void Get_Rotation(double RotMod, Polyhedron &Particle)
{
  //now we compute the rotations for each vetex of the particle.
   double e1[SD], e2[SD], e3[SD];
   Geometry Comput;

   //the rotation angle...
   double pr = (double)(rand()%MAXX)/(double)MAXX;
   double Phi = (pr-0.5)*RotMod*pi; 
   double CosPhi = cos(Phi);
   double SinPhi = sin(Phi);

   //the rotation axis...
   for(int k=0; k<SD; k++)
     e1[k] = (double)(rand()%MAXX)/(double)MAXX - 0.5;
   Comput.Normalize(e1); 
   //printf("e1 = (%f, %f, %f)\n", e1[0], e1[1], e1[2]);

   //for the vertices of the 
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   for(int i=0; i<Particle.n_vertex; i++)
     {
       double Vr[SD]; //the coordinates of a vertex
       for(int k=0; k<SD; k++)
	 Vr[k] = Particle.Vertex[i][k];
       
       double R1 = Comput.GetInnerProduct(Vr, e1); 

       //get the perpendicular components and e2
       for(int k=0; k<SD; k++)
	 e2[k] = Vr[k]-R1*e1[k]; 
       double tempd = Comput.Normalize(e2, 0);
       //printf("e2[%d] = (%f, %f, %f)\n", i, e2[0], e2[1], e2[2]);

       //get e3 as follows or 
       Comput.GetCrossProduct(e1, e2, e3);
       //printf("e3[%d] = (%f, %f, %f)\n", i, e3[0], e3[1], e3[2]);

       //get the other components, rotate the vector on the plane...
       double R2 = tempd*CosPhi;
       double R3 = tempd*SinPhi;

       //the new posistion of the Vertex
       Particle.Vertex[i][0] = R1*e1[0] + R2*e2[0] + R3*e3[0];
       Particle.Vertex[i][1] = R1*e1[1] + R2*e2[1] + R3*e3[1];
       Particle.Vertex[i][2] = R1*e1[2] + R2*e2[2] + R3*e3[2];

       //printf("Vertex[%d] = (%f, %f, %f)\n", i, TParticle.Vertex[i][0], TParticle.Vertex[i][1], TParticle.Vertex[i][2]);
       //printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

     } //end of loop over all the vertices 

   //needed for computing the volume of the particle
    Particle.GetNormal();
    

   //printf("*****************************************************\n");

}

void Read_Polyhedron(char* temp_name)
{
  //read in the shape, check vertex number, if don't match, report error
  //read in the vertex, initialize the packing
  //resacle the sphere packing
  //finish constructing the polypacking...

  fin.open("Poly_Shape.txt");
  if(!fin)
    {
      cout<<"Can not open file <Poly_Shape.txt>! Recheck!"<<endl;
      exit(1);
    }

  PackPoly = new Polyhedron[N]; //new [] cannot call constructors with parameters
  for(int i=0; i<N; i++)
    {
      PackPoly[i].PolyConstr(temp_name);
    }
  
  char poly_name[100];
  int temp_nvertex;

  fin>>poly_name; cout<<"Read in "<<poly_name<<endl;
  fin>>temp_nvertex;

  if(temp_nvertex != PackPoly[0].n_vertex)
    {
      cout<<"The vertex for the input shape is not correct! Recheck!"<<endl;
      exit(1);
    }
  

  //initalize the first particle with read-in; and copy it to the others
  for(int v=0; v<PackPoly[0].n_vertex; v++)
    for(int j=0; j<SD; j++)
      {	
	fin>>PackPoly[0].Vertex[v][j];
      }
  PackPoly[0].ShiftVert(); //shift vertices to make the centroid at orgin
  PackPoly[0].GetNormal(); //get the face normals
  PackPoly[0].GetCentDist(); //get d_cf and d_cv

  fin.close();

  for(int i=1; i<N; i++)
    {
      for(int v=0; v<PackPoly[0].n_vertex; v++)
	for(int j=0; j<SD; j++)
	  {	
	    PackPoly[i].Vertex[v][j] = PackPoly[0].Vertex[v][j];
	  }

      PackPoly[i].ShiftVert(); //shift vertices to make the centroid at orgin
      PackPoly[i].GetNormal(); //get the face normals
      PackPoly[i].GetCentDist(); //get d_cf and d_cv
    }


  //initialize d_cv_max with d_cv; and find the rescale_ratio
  d_cv_max = PackPoly[0].d_cv[0];

  rescale_ratio = 2*d_cv_max/d0;

  for(int i=0; i<SD; i++)
    for(int j=0; j<SD; j++)
      Lambda[i][j] = rescale_ratio*Lambda[i][j];  

  d_edge = 0.0;
  for(int i=0; i<SD; i++)
    d_edge += (PackPoly[0].Vertex[PackPoly[0].Edge[0][0]][i]-PackPoly[0].Vertex[PackPoly[0].Edge[0][1]][i])*(PackPoly[0].Vertex[PackPoly[0].Edge[0][0]][i]-PackPoly[0].Vertex[PackPoly[0].Edge[0][1]][i]);
  d_edge = sqrt(d_edge);

  //get random orientations for the shapes
  for(int i=0; i<N; i++)
    {
      Get_Rotation(pi, PackPoly[i]);
    }
  
}






void Print_PolyPacking()
{
  //get the current density, rescale to make a dilute packing if necessary
  //print out packing


  //check the target density if larger than that, further rescale the packing...
  double temp_rho = Get_Density();
  while(temp_rho > rho_tar)
    {
      for(int i=0; i<SD; i++)
	for(int j=0; j<SD; j++)
	  Lambda[i][j] = relax_ratio*Lambda[i][j];

      temp_rho = Get_Density();
    }

  cout<<"Current density rho = "<<temp_rho<<endl;
  cout<<"Target density rho = "<<rho_tar<<endl;


  //now print the packing...

  FILE* fp = fopen("Iconfig.txt","w");
  
  fprintf(fp, "%d\n", N);
  fprintf(fp, "%.12f\n", d_edge);
  
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
      
      for(int j=0; j<PackPoly[i].n_vertex; j++)
	fprintf(fp, "%.12f\t%.12f\t%.12f\n", PackPoly[i].Vertex[j][0], PackPoly[i].Vertex[j][1], PackPoly[i].Vertex[j][2]);
      
      fprintf(fp, "\n");
    }
  
  fclose(fp);
}




int main(int argc, char **argv)
{
  srand(time(NULL));

  
  cout<<"Please specify a target density for the packing: rho_tar = ";
  cin>>rho_tar;


  Read_SpherePacking();
  
  Get_CenterL();

  Read_Polyhedron(argv[1]);

  Get_VParticle();

  Print_PolyPacking();

  return 0;
}
