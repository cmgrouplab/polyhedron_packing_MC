
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
#include "Moves.h"


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
      TCenterL[k] = CenterLm[k] + dL[k];
      
      //update the relative position and make sure the particle is in the unit box...
      if(TCenterL[k]>=1.0) TCenterL[k] = TCenterL[k]-1.0;
      else if(TCenterL[k]<0) TCenterL[k] = TCenterL[k]+1.0;
    }

  //Note the cell which the particle belongs to may have changed, we will check it in MCmove()

}

void Moves::UpdateTranslation(double* &CenterLm)
{
  for(int k=0; k<SD; k++)
    {
      CenterLm[k] =  TCenterL[k]; //if successful, copy the temp. one to the real one...  
    }
}


/*
//manipulate each vertex, general enough for any shape
void Moves::GetRotation(double RotMod, Polyhedron &Particle)
{
  //now we compute the rotations for each vetex of the particle.
   double e1[SD], e2[SD], e3[SD];

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

   //first copy the current particle to Tparticle
   TParticle.CopyPoly(Particle);

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
       TParticle.Vertex[i][0] = R1*e1[0] + R2*e2[0] + R3*e3[0];
       TParticle.Vertex[i][1] = R1*e1[1] + R2*e2[1] + R3*e3[1];
       TParticle.Vertex[i][2] = R1*e1[2] + R2*e2[2] + R3*e3[2];

       //printf("Vertex[%d] = (%f, %f, %f)\n", i, TParticle.Vertex[i][0], TParticle.Vertex[i][1], TParticle.Vertex[i][2]);
       //printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

     } //end of loop over all the vertices 

   //IMPORTANT: Updat the Normal!!!!!
   //TParticle.GetNormal();

   //printf("*****************************************************\n");

}
*/

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
       
       //the new posistion of the Vertex
       TParticle.Vertex[i][0] = (CosPhi+e1[0]*e1[0]*(1-CosPhi))*Particle.Vertex[i][0] + (e1[0]*e1[1]*(1-CosPhi)-e1[2]*SinPhi)*Particle.Vertex[i][1] + (e1[0]*e1[2]*(1-CosPhi)+e1[1]*SinPhi)*Particle.Vertex[i][2];
       TParticle.Vertex[i][1] = (e1[0]*e1[1]*(1-CosPhi)+e1[2]*SinPhi)*Particle.Vertex[i][0] + (CosPhi+e1[1]*e1[1]*(1-CosPhi))*Particle.Vertex[i][1] + (e1[1]*e1[2]*(1-CosPhi)-e1[0]*SinPhi)*Particle.Vertex[i][2];
       TParticle.Vertex[i][2] = (e1[0]*e1[2]*(1-CosPhi)-e1[1]*SinPhi)*Particle.Vertex[i][0] + (e1[1]*e1[2]*(1-CosPhi)+e1[0]*SinPhi)*Particle.Vertex[i][1] + (CosPhi+e1[2]*e1[2]*(1-CosPhi))*Particle.Vertex[i][2];
       
       //printf("Vertex[%d] = (%f, %f, %f)\n", i, TParticle.Vertex[i][0], TParticle.Vertex[i][1], TParticle.Vertex[i][2]);
       //printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

     } //end of loop over all the vertices 


   /*
   pr = (double)(rand()%MAXX)/(double)MAXX;
   double Theta = pr*RotMod*pi; 
   double CosTheta = cos(Theta);
   double SinTheta = sin(Theta);

   pr = (double)(rand()%MAXX)/(double)MAXX;
   double Psi = (pr-0.5)*RotMod*pi; 
   double CosPsi = cos(Psi);
   double SinPsi = sin(Psi);

   //first copy the current particle to Tparticle
   //TParticle.CopyPoly(Particle);

   //for the vertices of the 
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   for(int i=0; i<Particle.n_vertex; i++)
     {
     
       //the new posistion of the Vertex
       TParticle.Vertex[i][0] = CosPhi*CosPsi*Particle.Vertex[i][0] + (-CosPhi*SinPsi+SinPhi*SinTheta*CosPsi)*Particle.Vertex[i][1] + (SinPhi*SinPsi+CosPhi*SinTheta*CosPsi)*Particle.Vertex[i][2];
       TParticle.Vertex[i][1] = CosPhi*SinPsi*Particle.Vertex[i][0] + (CosPhi*CosPsi+SinPhi*SinTheta*SinPsi)*Particle.Vertex[i][1] + (-SinPhi*CosPsi+CosPhi*SinTheta*SinPsi)*Particle.Vertex[i][2];
       TParticle.Vertex[i][2] = -SinTheta*Particle.Vertex[i][0] + SinPhi*CosTheta*Particle.Vertex[i][1] + CosPhi*CosTheta*Particle.Vertex[i][2];
       
       //printf("Vertex[%d] = (%f, %f, %f)\n", i, TParticle.Vertex[i][0], TParticle.Vertex[i][1], TParticle.Vertex[i][2]);
       //printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

     } //end of loop over all the vertices 
   */

   //finally, get the center to vertices and center to face distances
   //these are necessary for overlap check...
   for(int i=0; i<Particle.n_face; i++)
     TParticle.d_cf[i] = Particle.d_cf[i];
  
   for(int i=0; i<Particle.n_vertex; i++)
     TParticle.d_cv[i] = Particle.d_cv[i];
   
   //IMPORTANT: Updat the Normal!!!!!
   //TParticle.GetNormal();

   //printf("*****************************************************\n");

}





//no need to copy everything again, only the updated quantities
void Moves::UpdateRotation(Polyhedron &Particle)
{
  for(int i=0; i<Particle.n_vertex; i++)
    {
      for(int k=0; k<SD; k++)
	{
	   Particle.Vertex[i][k] = TParticle.Vertex[i][k];
	}
    }
  
  /*
  for(int i=0; i<Particle.n_face; i++)
    {
      for(int k=0; k<SD; k++)
	{
	   Particle.Normal[i][k] = TParticle.Normal[i][k];
	}
    }
  */
}


//each move could be a translation or a rotation with specific probability
//return 1 if the trial traslation is successful
//retrun 2 if the trial rotation is sucessful
//return 0 if fails....
int Moves::ParticleMove(double p_trans, double TransMod, double RotMod, int m, Packing &PK, Cells &BoxCell)
{
  double pm = (double)(rand()%MAXX)/(double)MAXX;

  if(pm<=p_trans) //translational trail move
    {
      //translate the particle
      GetTranslation(TransMod, PK.CenterL[m]);

      //printf("Translation is performed!\n");

      //need to check overlapping...
      int temp_indI = (int)floor(TCenterL[0]*BoxCell.LXcell);
      int temp_indJ = (int)floor(TCenterL[1]*BoxCell.LYcell);
      int temp_indK = (int)floor(TCenterL[2]*BoxCell.LZcell);

      for(int i=-1; i<=1; i++)
	for(int j=-1; j<=1; j++)
	  for(int k=-1; k<=1; k++)
	    {
	      //for the box index of temp_P, they are the same for each cell
	      int indexI = 0;
	      int indexJ = 0; 
	      int indexK = 0;

	      int tempI = temp_indI + i;
	      if(tempI>=BoxCell.LXcell) 
		{
		  tempI = tempI-BoxCell.LXcell;
		  indexI = 1;
		}
	      else if(tempI<0) 
		{
		  tempI = tempI+BoxCell.LXcell;
		  indexI = -1;
		}

	      int tempJ = temp_indJ + j;
	      if(tempJ>=BoxCell.LYcell) 
		{
		  tempJ = tempJ-BoxCell.LYcell;
		  indexJ = 1;
		}
	      else if(tempJ<0) 
		{
		  tempJ = tempJ+BoxCell.LYcell; 
		  indexJ = -1;
		}

	      int tempK = temp_indK + k;
	      if(tempK>=BoxCell.LZcell) 
		{ 
		  tempK = tempK-BoxCell.LZcell; 
		  indexK = 1;
		}
	      else if(tempK<0) 
		{ 
		  tempK = tempK+BoxCell.LZcell; 
		  indexK = -1;
		}

	      int cell_index = BoxCell.LXcell*BoxCell.LYcell*tempK + BoxCell.LXcell*tempJ + tempI;

	      node* pt = BoxCell.CellList[cell_index]; //node is defined in const.h

	      while(pt!=NULL)
		{
		  int temp_P = pt->index_poly;

		  //the particles overlap
		  //the position of m has been updated, have to exclude the particle when checking
		  if(temp_P!=m && PK.CheckOverlap(PK.PackPoly[temp_P], PK.PackPoly[m], PK.CenterL[temp_P], TCenterL, indexI, indexJ, indexK)==1) 
		    {		     
		      return 0; //indicating an unsuccessful trial...
		    }

		  pt = pt->next;
		}
	    }//end of loop over all neighboring cells...

      //if we reach here, then means no overlap...
      //only need to check and update the position and the celllist...
         
      //for the old position
      int temp2_indI = (int)floor(PK.CenterL[m][0]*BoxCell.LXcell);
      int temp2_indJ = (int)floor(PK.CenterL[m][1]*BoxCell.LYcell);
      int temp2_indK = (int)floor(PK.CenterL[m][2]*BoxCell.LZcell);

      if(temp2_indI==temp_indI && temp2_indJ==temp_indJ && temp2_indK==temp_indK)
	{
	  //no need to update the CellList
	  UpdateTranslation(PK.CenterL[m]);
	  return 1;
	}
      else
	{
	  BoxCell.UpdateCellList(m, PK.CenterL[m], TCenterL);
	  UpdateTranslation(PK.CenterL[m]);
	  return 1;
	}
        //indicating a successful trail move....

    } //ENDIF
  else //rotational trail move...
    {
      //make a rotation...
      GetRotation(RotMod, PK.PackPoly[m]);

      //printf("Rotation is performed!\n");

      //need to check overlapping...
      int temp_indI = (int)floor(PK.CenterL[m][0]*BoxCell.LXcell);
      int temp_indJ = (int)floor(PK.CenterL[m][1]*BoxCell.LYcell);
      int temp_indK = (int)floor(PK.CenterL[m][2]*BoxCell.LZcell);

      for(int i=-1; i<=1; i++)
	for(int j=-1; j<=1; j++)
	  for(int k=-1; k<=1; k++)
	    {
	      int indexI = 0;
	      int indexJ = 0;
	      int indexK = 0;
	      
	      int tempI = temp_indI + i;
	      if(tempI>=BoxCell.LXcell) 
		{
		  tempI = tempI-BoxCell.LXcell;
		  indexI = 1;
		}
	      else if(tempI<0) 
		{
		  tempI = tempI+BoxCell.LXcell;
		  indexI = -1;
		}

	      int tempJ = temp_indJ + j;
	      if(tempJ>=BoxCell.LYcell) 
		{
		  tempJ = tempJ-BoxCell.LYcell;
		  indexJ = 1;
		}
	      else if(tempJ<0) 
		{
		  tempJ = tempJ+BoxCell.LYcell; 
		  indexJ = -1;
		}

	      int tempK = temp_indK + k;
	      if(tempK>=BoxCell.LZcell) 
		{ 
		  tempK = tempK-BoxCell.LZcell; 
		  indexK = 1;
		}
	      else if(tempK<0) 
		{ 
		  tempK = tempK+BoxCell.LZcell; 
		  indexK = -1;
		}
	      
	      int cell_index = BoxCell.LXcell*BoxCell.LYcell*tempK + BoxCell.LXcell*tempJ + tempI;

	      node* pt = BoxCell.CellList[cell_index];
	     
	      while(pt!=NULL)
		{
		  int temp_P = pt->index_poly;

		  //the particles overlap
		  //the particle configuration is updated
		  if(temp_P!=m && PK.CheckOverlap(PK.PackPoly[temp_P], TParticle, PK.CenterL[temp_P], PK.CenterL[m], indexI, indexJ, indexK)==1)
		    {
		      return 0; //indicating an unsuccessful trial...
		    }

		  pt = pt->next;
		}
	    }//end of loop over all neighboring cells...
      

      UpdateRotation(PK.PackPoly[m]);
      return 2; //indicating a successful rotation trail 


    }//ENDELSE

}



//using SEP_List and rigorous cell method
int Moves::ParticleMoveII(double p_trans, double TransMod, double RotMod, int m, Packing &PK, Cells &BoxCell)
{
  double pm = (double)(rand()%MAXX)/(double)MAXX;

  if(pm<=p_trans) //translational trail move
    {
      //translate the particle
      GetTranslation(TransMod, PK.CenterL[m]);

      //printf("Translation is performed!\n");

      //need to check overlapping...
      int temp_indI = (int)floor(TCenterL[0]*BoxCell.LXcell);
      int temp_indJ = (int)floor(TCenterL[1]*BoxCell.LYcell);
      int temp_indK = (int)floor(TCenterL[2]*BoxCell.LZcell);

      for(int i=-1; i<=1; i++)
	for(int j=-1; j<=1; j++)
	  for(int k=-1; k<=1; k++)
	    {
	      //for the box index of temp_P, they are the same for each cell
	      int indexI = 0;
	      int indexJ = 0; 
	      int indexK = 0;

	      int tempI = temp_indI + i;
	      if(tempI>=BoxCell.LXcell) 
		{
		  tempI = tempI-BoxCell.LXcell;
		  indexI = 1;
		}
	      else if(tempI<0) 
		{
		  tempI = tempI+BoxCell.LXcell;
		  indexI = -1;
		}

	      int tempJ = temp_indJ + j;
	      if(tempJ>=BoxCell.LYcell) 
		{
		  tempJ = tempJ-BoxCell.LYcell;
		  indexJ = 1;
		}
	      else if(tempJ<0) 
		{
		  tempJ = tempJ+BoxCell.LYcell; 
		  indexJ = -1;
		}

	      int tempK = temp_indK + k;
	      if(tempK>=BoxCell.LZcell) 
		{ 
		  tempK = tempK-BoxCell.LZcell; 
		  indexK = 1;
		}
	      else if(tempK<0) 
		{ 
		  tempK = tempK+BoxCell.LZcell; 
		  indexK = -1;
		}

	      int cell_index = BoxCell.LXcell*BoxCell.LYcell*tempK + BoxCell.LXcell*tempJ + tempI;

	      node* pt = BoxCell.CellList[cell_index]; //node is defined in const.h

	      while(pt!=NULL)
		{
		  int temp_P = pt->index_poly;

		  //the particles overlap
		  //the position of m has been updated, have to exclude the particle when checking
		  if(temp_P!=m && PK.CheckOverlap(PK.PackPoly[temp_P], PK.PackPoly[m], PK.CenterL[temp_P], TCenterL, indexI, indexJ, indexK, m, temp_P)==1) 
		    {		     
		      return 0; //indicating an unsuccessful trial...
		    }

		  pt = pt->next;
		}
	    }//end of loop over all neighboring cells...

      //if we reach here, then means no overlap...
      //only need to check and update the position and the celllist...
         
      //for the old position
      int temp2_indI = (int)floor(PK.CenterL[m][0]*BoxCell.LXcell);
      int temp2_indJ = (int)floor(PK.CenterL[m][1]*BoxCell.LYcell);
      int temp2_indK = (int)floor(PK.CenterL[m][2]*BoxCell.LZcell);

      if(temp2_indI==temp_indI && temp2_indJ==temp_indJ && temp2_indK==temp_indK)
	{
	  //no need to update the CellList
	  UpdateTranslation(PK.CenterL[m]);
	  return 1;
	}
      else
	{
	  BoxCell.UpdateCellList(m, PK.CenterL[m], TCenterL);
	  UpdateTranslation(PK.CenterL[m]);
	  return 1;
	}
        //indicating a successful trail move....

    } //ENDIF
  else //rotational trail move...
    {
      //make a rotation...
      GetRotation(RotMod, PK.PackPoly[m]);

      //printf("Rotation is performed!\n");

      //need to check overlapping...
      int temp_indI = (int)floor(PK.CenterL[m][0]*BoxCell.LXcell);
      int temp_indJ = (int)floor(PK.CenterL[m][1]*BoxCell.LYcell);
      int temp_indK = (int)floor(PK.CenterL[m][2]*BoxCell.LZcell);

      for(int i=-1; i<=1; i++)
	for(int j=-1; j<=1; j++)
	  for(int k=-1; k<=1; k++)
	    {
	      int indexI = 0;
	      int indexJ = 0;
	      int indexK = 0;
	      
	      int tempI = temp_indI + i;
	      if(tempI>=BoxCell.LXcell) 
		{
		  tempI = tempI-BoxCell.LXcell;
		  indexI = 1;
		}
	      else if(tempI<0) 
		{
		  tempI = tempI+BoxCell.LXcell;
		  indexI = -1;
		}

	      int tempJ = temp_indJ + j;
	      if(tempJ>=BoxCell.LYcell) 
		{
		  tempJ = tempJ-BoxCell.LYcell;
		  indexJ = 1;
		}
	      else if(tempJ<0) 
		{
		  tempJ = tempJ+BoxCell.LYcell; 
		  indexJ = -1;
		}

	      int tempK = temp_indK + k;
	      if(tempK>=BoxCell.LZcell) 
		{ 
		  tempK = tempK-BoxCell.LZcell; 
		  indexK = 1;
		}
	      else if(tempK<0) 
		{ 
		  tempK = tempK+BoxCell.LZcell; 
		  indexK = -1;
		}
	      
	      int cell_index = BoxCell.LXcell*BoxCell.LYcell*tempK + BoxCell.LXcell*tempJ + tempI;

	      node* pt = BoxCell.CellList[cell_index];
	     
	      while(pt!=NULL)
		{
		  int temp_P = pt->index_poly;

		  //the particles overlap
		  //the particle configuration is updated
		  if(temp_P!=m && PK.CheckOverlap(PK.PackPoly[temp_P], TParticle, PK.CenterL[temp_P], PK.CenterL[m], indexI, indexJ, indexK, m, temp_P)==1)
		    {
		      return 0; //indicating an unsuccessful trial...
		    }

		  pt = pt->next;
		}
	    }//end of loop over all neighboring cells...
      

      UpdateRotation(PK.PackPoly[m]);
      return 2; //indicating a successful rotation trail 


    }//ENDELSE

}





 

//each move could be a translation or a rotation with specific probability
//return 1 if the trial traslation is successful
//retrun 2 if the trial rotation is sucessful
//return 0 if fails....
 /*
//this uses NNL and SEP_List...
int Moves::ParticleMoveII(double p_trans, double TransMod, double RotMod, int m, Packing &PK, Cells &BoxCell)
{
  double pm = (double)(rand()%MAXX)/(double)MAXX;

  if(pm<=p_trans) //translational trail move
    {
      //translate the particle
      GetTranslation(TransMod, PK.CenterL[m]);

      //printf("Translation is performed!\n");

      for(int i=0; i<PK.NNL_counter[m]; i++)
	{
	  int temp_P = PK.NNL[m][i][0];

	  int indexI = PK.NNL[m][i][1];
	  int indexJ = PK.NNL[m][i][2]; 
	  int indexK = PK.NNL[m][i][3];

	  if(temp_P!=m && PK.CheckOverlap(PK.PackPoly[temp_P], PK.PackPoly[m], PK.CenterL[temp_P], TCenterL, indexI, indexJ, indexK, m, temp_P)==1) 
	    {		     
	      return 0; //indicating an unsuccessful trial...
	    }
	  
	}

      //if we reach here, then means no overlap...
      //only need to check and update the position and the celllist...

      int temp_indI = (int)floor(TCenterL[0]*BoxCell.LXcell);
      int temp_indJ = (int)floor(TCenterL[1]*BoxCell.LYcell);
      int temp_indK = (int)floor(TCenterL[2]*BoxCell.LZcell);
         
      //for the old position
      int temp2_indI = (int)floor(PK.CenterL[m][0]*BoxCell.LXcell);
      int temp2_indJ = (int)floor(PK.CenterL[m][1]*BoxCell.LYcell);
      int temp2_indK = (int)floor(PK.CenterL[m][2]*BoxCell.LZcell);

      if(temp2_indI==temp_indI && temp2_indJ==temp_indJ && temp2_indK==temp_indK)
	{
	  //no need to update the CellList
	  UpdateTranslation(PK.CenterL[m]);
	  return 1;
	}
      else
	{
	  BoxCell.UpdateCellList(m, PK.CenterL[m], TCenterL);
	  UpdateTranslation(PK.CenterL[m]);
	  return 1;
	}
        //indicating a successful trail move....

    } //ENDIF
  else //rotational trail move...
    {
      //make a rotation...
      GetRotation(RotMod, PK.PackPoly[m]);

      //printf("Rotation is performed!\n");
      for(int i=0; i<PK.NNL_counter[m]; i++)
	{
	  int temp_P = PK.NNL[m][i][0];
	  
	  int indexI = PK.NNL[m][i][1];
	  int indexJ = PK.NNL[m][i][2]; 
	  int indexK = PK.NNL[m][i][3];
	  
	  if(temp_P!=m && PK.CheckOverlap(PK.PackPoly[temp_P], TParticle, PK.CenterL[temp_P], PK.CenterL[m], indexI, indexJ, indexK, m, temp_P)==1)
    	    {		     
	      return 0; //indicating an unsuccessful trial...
	    }
	  
	}
      
      
      UpdateRotation(PK.PackPoly[m]);
      return 2; //indicating a successful rotation trail 
      
      
    }//ENDELSE
  
}
 */



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
int Moves::BoundaryMove(int Ndeform, double p_uphill, double Strain_Rescale, double Global_StrainMod, double Shear_StrainMod, Packing &PK, Cells &BoxCell)
{
  //Now we deform the boundary and report density...
  GetStrainIso(p_uphill, Global_StrainMod, Shear_StrainMod); //each time we use the orginally specified strain rate
  
  BoundaryDeform(PK.Lambda);
  
  int overlap_flag = PK.GlobalOverlapCheck(BoxCell);
  int limit = 0;
  
  while(overlap_flag==1 && limit<Ndeform)
    {
      BoundaryRetain(PK.Lambda); //retain the old boundary...
      
      ResetStrain(p_uphill, Strain_Rescale, limit, Global_StrainMod, Shear_StrainMod);
      
      //GetStrainIso();
      
      BoundaryDeform(PK.Lambda);
      
      overlap_flag = PK.GlobalOverlapCheck(BoxCell);
      
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


