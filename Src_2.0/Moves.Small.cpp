
using namespace std;

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "const.h"
#include "Geometry.h"
#include "Polyhedron.h"
#include "Packing.h"
//#include "Cells.h"
#include "Moves.Small.h"


Moves::Moves(char* temp_name)
{
  srand(time(NULL));

  TParticle.PolyConstr(temp_name); //now provide the shape info.
}


//work with the relative coordinates,
void Moves::GetTranslation(double TransMod, double* &CenterLm)
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
      TCenterL[k] = CenterLm[k];
      CenterLm[k] = CenterLm[k] + dL[k];
      
      //update the relative position and make sure the particle is in the unit box...
      if(CenterLm[k]>=1.0) CenterLm[k] = CenterLm[k]-1.0;
      else if(CenterLm[k]<0) CenterLm[k] = CenterLm[k]+1.0;
    }

  //Note the cell which the particle belongs to may have changed, we will check it in MCmove()

}

void Moves::RetainTranslation(double* &CenterLm)
{
  for(int k=0; k<SD; k++)
    {
      CenterLm[k] =  TCenterL[k]; //if successful, copy the temp. one to the real one...  
    }
}


//now we try a different way of getting rotations, using the rotation matrix
//IMPORTANT!!!!! Not tested before, probably not correct...
void Moves::GetRotation(double RotMod, Polyhedron &Particle)
{
   //the rotation angles...
   //computing sine and cosine are not efficient...
   double pr = (double)(rand()%MAXX)/(double)MAXX;
   double Phi = (pr-0.5)*RotMod*pi; 
   double CosPhi = cos(Phi);
   double SinPhi = sin(Phi);


   //the rotation axis...
   double e1[SD];
   for(int k=0; k<SD; k++)
     e1[k] = (double)(rand()%MAXX)/(double)MAXX - 0.5;
   Comput.Normalize(e1); 

    //first copy the current particle to Tparticle
   // TParticle.CopyPoly(Particle);


   for(int i=0; i<Particle.n_vertex; i++)
     {
       for(int k=0; k<SD; k++)
	 TParticle.Vertex[i][k] = Particle.Vertex[i][k];

       //the new posistion of the Vertex
       Particle.Vertex[i][0] = (CosPhi+e1[0]*e1[0]*(1-CosPhi))*TParticle.Vertex[i][0] + (e1[0]*e1[1]*(1-CosPhi)-e1[2]*SinPhi)*TParticle.Vertex[i][1] + (e1[0]*e1[2]*(1-CosPhi)+e1[1]*SinPhi)*TParticle.Vertex[i][2];
       Particle.Vertex[i][1] = (e1[0]*e1[1]*(1-CosPhi)+e1[2]*SinPhi)*TParticle.Vertex[i][0] + (CosPhi+e1[1]*e1[1]*(1-CosPhi))*TParticle.Vertex[i][1] + (e1[1]*e1[2]*(1-CosPhi)-e1[0]*SinPhi)*TParticle.Vertex[i][2];
       Particle.Vertex[i][2] = (e1[0]*e1[2]*(1-CosPhi)-e1[1]*SinPhi)*TParticle.Vertex[i][0] + (e1[1]*e1[2]*(1-CosPhi)+e1[0]*SinPhi)*TParticle.Vertex[i][1] + (CosPhi+e1[2]*e1[2]*(1-CosPhi))*TParticle.Vertex[i][2];
       
       //printf("Vertex[%d] = (%f, %f, %f)\n", i, TParticle.Vertex[i][0], TParticle.Vertex[i][1], TParticle.Vertex[i][2]);
       //printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

     } //end of loop over all the vertices 

   /*
   //finally, get the center to vertices and center to face distances
   //these are necessary for overlap check...
   for(int i=0; i<Particle.n_face; i++)
     TParticle.d_cf[i] = Particle.d_cf[i];
  
   for(int i=0; i<Particle.n_vertex; i++)
     TParticle.d_cv[i] = Particle.d_cv[i];
   */

   //IMPORTANT: Updat the Normal!!!!!
   //TParticle.GetNormal();

   //printf("*****************************************************\n");

}





//no need to copy everything again, only the updated quantities
void Moves::RetainRotation(Polyhedron &Particle)
{
  for(int i=0; i<Particle.n_vertex; i++)
    {
      for(int k=0; k<SD; k++)
	{
	   Particle.Vertex[i][k] = TParticle.Vertex[i][k];
	}
    }
}


//each move could be a translation or a rotation with specific probability
//return 1 if the trial traslation is successful
//retrun 2 if the trial rotation is sucessful
//return 0 if fails....
int Moves::ParticleMove(double p_trans, double TransMod, double RotMod, int m, Packing &PK)
{
  double pm = (double)(rand()%MAXX)/(double)MAXX;

  if(pm<=p_trans) //translational trail move
    {
      //translate the particle
      GetTranslation(TransMod, PK.CenterL[m]);

      //printf("Translation is performed!\n");

      if(PK.GlobalOverlapCheck(1.0) == 1)
	{
	  RetainTranslation(PK.CenterL[m]);  
	  return 0;
	}
      else
	{
	  //UpdateTranslation(PK.CenterL[m]);  
	  return 1;
	}
        //indicating a successful trail move...
    } //ENDIF
  else //rotational trail move...
    {
      //make a rotation...
      GetRotation(RotMod, PK.PackPoly[m]);

      //printf("Rotation is performed!\n");
      if(PK.GlobalOverlapCheck(1.0) == 1)
	{
	  RetainRotation(PK.PackPoly[m]);
	  return 0;
	}
      else
	{
	  return 2; //indicating a successful rotation trail 
	}     
      
    }//ENDELSE
  
}


void Moves::GetStrainIso(double p_uphill, double P_StrainMod, double S_StrainMod) //modified to allow uphill moves...
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


void Moves::ResetStrain(double p_uphill, double Strain_Rescale, int H,  double Global_StrainMod, double Shear_StrainMod) //decrease the magnitude by certain amount..
{
  double Temp_PStrainMod = pow(Strain_Rescale, H)*Global_StrainMod;
  double Temp_SStrainMod = pow(Strain_Rescale, H)*Shear_StrainMod;

  GetStrainIso(p_uphill, Temp_PStrainMod, Temp_SStrainMod);

}

void Moves::BoundaryDeform(double** &Lambda)
{
  //store the old lambda first...
  for(int i=0; i<SD; i++)
    for(int j=0; j<SD; j++)
      TLambda[i][j] = Lambda[i][j];


  //chack the boundary...
  double D_Lambda[SD][SD] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};

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


void Moves::BoundaryRetain(double** &Lambda)
{
  //store the old lambda first...
  for(int i=0; i<SD; i++)
    for(int j=0; j<SD; j++)
      Lambda[i][j] = TLambda[i][j];
}


//if a success, return 1; otherwise return 0
int Moves::BoundaryMove(int Ndeform, double p_uphill, double Strain_Rescale, double Global_StrainMod, double Shear_StrainMod, Packing &PK)
{
  //Now we deform the boundary and report density...
  GetStrainIso(p_uphill, Global_StrainMod, Shear_StrainMod); //each time we use the orginally specified strain rate
  
  BoundaryDeform(PK.Lambda);
  
  int overlap_flag = PK.GlobalOverlapCheck(1);
  int limit = 0;
  
  while(overlap_flag==1 && limit<Ndeform)
    {
      BoundaryRetain(PK.Lambda); //retain the old boundary...
      
      ResetStrain(p_uphill, Strain_Rescale, limit, Global_StrainMod, Shear_StrainMod);
      
      //GetStrainIso();
      
      BoundaryDeform(PK.Lambda);
      
      overlap_flag = PK.GlobalOverlapCheck(1);
      
      limit++;
    }
  
  if(limit==Ndeform) 
    {
      BoundaryRetain(PK.Lambda);
      printf("Failed to compress at this stage \n");
      return 0;
    }
  else
    return 1;
}


