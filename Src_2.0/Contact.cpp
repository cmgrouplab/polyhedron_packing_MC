

using namespace std;

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "const.h"
#include "Geometry.h"
#include "Polyhedron.h"
//#include "Cells.h"
#include "Packing.h"
#include "Contact.h"



Contact::Contact(char* temp_name, Packing& PK)
{
  //~~~~~~~~~~~~~~~~~~~~~~~~
  //first, initialize all the corresponding variables...
  N = PK.N;
  d0 = PK.d0;

  Lambda = new double*[SD];
  for(int i=0; i<SD; i++)
    Lambda[i] = new double[SD];

  for(int i=0; i<SD; i++) // read in the lattice vectors...
    for(int j=0; j<SD; j++)
      {
	Lambda[i][j] = PK.Lambda[i][j];
      }

  CenterE = new double*[N];
  CenterL = new double*[N];

  PackPoly = new Polyhedron[N]; //new [] cannot call constructors with parameters
      
  for(int i=0; i<N; i++)
    {
      CenterE[i] = new double[SD];	  
      CenterL[i] = new double[SD];
      
      PackPoly[i].PolyConstr(temp_name);
      
    }
  
  //now get the values from Packing PK...
  for(int i=0; i<N; i++) // read the posistions of the center and the vertex
    {
      //read in the center of mass
      for(int j=0; j<SD; j++)
	{
	  CenterL[i][j] = PK.CenterL[i][j]; 
	  if(CenterL[i][j]>=1.0) CenterL[i][j] = CenterL[i][j]-1.0;
	  else if(CenterL[i][j]<0) CenterL[i][j] = CenterL[i][j]+1.0;
	}
    
	  
      
      //read in vertex with O = 0; and to make sure, shift the vertices again for each particle
      
      for(int v=0; v<PackPoly[i].n_vertex; v++)
	for(int j=0; j<SD; j++)
	  {
	    PackPoly[i].Vertex[v][j] = PK.PackPoly[i].Vertex[v][j];
	  }
      
	 
      PackPoly[i].ShiftVert(); //shift vertices to make the centroid at orgin
      
      PackPoly[i].GetNormal_All(); //get the face normals
            
      PackPoly[i].GetCentDist(); //get d_cf and d_cv
      
      
      
      //get the maximum characteristic length of the particle... 
      double temp_d = 2*GetMax(PackPoly[i].d_cv, PackPoly[i].n_vertex);
      if(temp_d>d0) d0 = temp_d;
      
    }//finishing loop over all particles

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //now we deal with other variable initializations...
  
  Contact_Array = new int*[N];
  
  NP_f2f = new int[N]; NP_e2f = new int[N]; NP_p2f = new int[N]; NP_e2e = new int[N];
  
  n_f2f = 0; n_e2f = 0; n_p2f = 0; n_e2e = 0;

  for(int i=0; i<N; i++)
    {
      Contact_Array[i] = new int[N];

      NP_f2f[i] = 0;  NP_e2f[i] = 0; NP_p2f[i] = 0; NP_e2e[i] = 0;

      for(int j=0; j<N; j++)
	Contact_Array[i][j] = -1;
    }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //now ask for the tolerance 
  TOL = 0.01;
  cout<<"Specify TOL for neighbor cut-off <give a value> OR use default "<<TOL<<" <give -1>?"<<endl;
  double temp_val;
  cin>>temp_val;
  if(temp_val>0)
    {
      TOL = temp_val;
      cout<<"Specified TOL = "<<TOL<<" *d0 for cut-off"<<endl;
    }
  else
    {
      cout<<"Default TOL = "<<TOL<<" *d0 for cut-off"<<endl;
    }
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

  radial_cutoff = 1.0; //see Contact.h file for explanation 
  normal_cutoff = 0.97;
  f2f_cutoff = 0.5;
  e2f_cutoff = 0.5;
  v2f_cutoff = 0.5;
  edge_cutoff = 0.325;

  double temp_cutoff;
  cout<<"Specify radial_cutoff (for determing radial neighbors) <give a number> or default "<<radial_cutoff<<" <-1> "<<endl;
  cin>>temp_cutoff;
  if(temp_cutoff == -1)
    {
      cout<<"Using default radial_cutoff = "<<radial_cutoff<<endl;
    }
  else
    {
      radial_cutoff = temp_cutoff;
      cout<<"The specified radial_cutoff = "<<radial_cutoff<<endl;
    }
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

  cout<<"Specify normal_cutoff (for determing f2f cont. based on normal inner product.) <give a number> or default "<<normal_cutoff<<" <-1> "<<endl;
  cin>>temp_cutoff;
  if(temp_cutoff == -1)
    {
      cout<<"Using default normal_cutoff = "<<normal_cutoff<<endl;
    }
  else
    {
      normal_cutoff = temp_cutoff;
      cout<<"The specified radial_cutoff = "<<normal_cutoff<<endl;
    }
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

  cout<<"Specify f2f_cutoff (for determing face cont. based on ct sep dist.) <give a number> or default "<<f2f_cutoff<<" <-1> "<<endl;
  cin>>temp_cutoff;
  if(temp_cutoff == -1)
    {
      cout<<"Using default f2f_cutoff = "<<f2f_cutoff<<endl;
    }
  else
    {
      f2f_cutoff = temp_cutoff;
      cout<<"The specified f2f_cutoff = "<<f2f_cutoff<<endl;
    }
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

  cout<<"Specify e2f_cutoff (for determing face cont. based on ct sep dist.) <give a number> or default "<<e2f_cutoff<<" <-1> "<<endl;
  cin>>temp_cutoff;
  if(temp_cutoff == -1)
    {
      cout<<"Using default e2f_cutoff = "<<e2f_cutoff<<endl;
    }
  else
    {
      e2f_cutoff = temp_cutoff;
      cout<<"The specified e2f_cutoff = "<<e2f_cutoff<<endl;
    }
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

  cout<<"Specify v2f_cutoff (for determing face cont. based on ct sep dist.) <give a number> or default "<<v2f_cutoff<<" <-1> "<<endl;
  cin>>temp_cutoff;
  if(temp_cutoff == -1)
    {
      cout<<"Using default v2f_cutoff = "<<v2f_cutoff<<endl;
    }
  else
    {
      v2f_cutoff = temp_cutoff;
      cout<<"The specified v2f_cutoff = "<<v2f_cutoff<<endl;
    }
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

  cout<<"Specify edge_cutoff (for determing e2e cont based on edge ct dist.) <give a number> or default "<<edge_cutoff<<" <-1> "<<endl;
  cin>>temp_cutoff;
  if(temp_cutoff == -1)
    {
      cout<<"Using default edge_cutoff = "<<edge_cutoff<<endl;
    }
  else
    {
      edge_cutoff = temp_cutoff;
      cout<<"The specified edge_cutoff = "<<edge_cutoff<<endl;
    }
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

  
}


void Contact::GetGlobalPosition()
{
  //compute the global posistions of the center of mass
  for(int n=0; n<N; n++)
    {
      for(int i=0; i<SD; i++)
	{
	  
	  CenterE[n][i] = Lambda[i][0]*CenterL[n][0]+Lambda[i][1]*CenterL[n][1]+Lambda[i][2]*CenterL[n][2];
	}
      
    }
}

double Contact::MinDis(double CenterLm[SD], double CenterLn[SD], double temp_disE[SD], int indexI, int indexJ, int indexK)
{
  //find the minimal distance between the centers of two polyhedra in Ecludean space...
  //by checking all the images of Polyhedron[m], while keeping Polyhedron[n] in the central box
  //record the index of the box which Tetrah[m] is in...
  //the order (m,n) is important, especially when getting the NNL.... 
  
  double tempLd[SD];

  tempLd[0] = (CenterLm[0]-CenterLn[0])+indexI;
  tempLd[1] = (CenterLm[1]-CenterLn[1])+indexJ;
  tempLd[2] = (CenterLm[2]-CenterLn[2])+indexK;
  
  for(int q=0; q<SD; q++)
    temp_disE[q] = 0;
  
  for(int q=0; q<SD; q++)
    for(int h=0; h<SD; h++)
      temp_disE[q] += Lambda[q][h]*tempLd[h];
  
  double dist = temp_disE[0]*temp_disE[0]+temp_disE[1]*temp_disE[1]+temp_disE[2]*temp_disE[2]; 

  //printf("Dist = %f\n", dist);

  return sqrt(dist); //this is center-to-center distance...
}

double Contact::GetMin(double Arr[], int Size)
{
  double temp_min = 1000000000.0;

  for(int i=0; i<Size; i++)
    {
      if(Arr[i]<temp_min) temp_min = Arr[i];
    }

  return temp_min;
}

double Contact::GetMax(double Arr[], int Size)
{
  double temp_max = -1000000000.0;

  for(int i=0; i<Size; i++)
    {
      if(Arr[i]>temp_max) temp_max = Arr[i];
    }

  return temp_max;
}

//try to generalize this to check all shapes
//make sure the particle do not check its old position for overlap
//return: -1 if particles do not contact
//        -2 if overlap
//         1 for f2f
//         2 for e2f
//         3 for p2f
//         4 for e2e


int Contact::GetContact(Polyhedron &Pm, Polyhedron &Pn, double CenterLm[SD], double CenterLn[SD], int indexI, int indexJ, int indexK)
{
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  double temp_VdisE[SD]; //the relative and Ecludean separation distance...
  
  //first get the center-to-center distance of the particles...
  double temp_dis = MinDis(CenterLm, CenterLn, temp_VdisE, indexI, indexJ, indexK);
  //recall n is fixed in the central box, the images of m are checked...

  //the mininmal and maximal possible separations
  double d_cvmax = radial_cutoff*(GetMax(Pm.d_cv, Pm.n_vertex)+GetMax(Pn.d_cv, Pn.n_vertex));
  //!!!!NOTE: the relaxation here <0.75> is empirical: for two particles contact, their centroids have to be close enough

  double d_cfmin = (GetMin(Pm.d_cf, Pm.n_face)+GetMin(Pn.d_cf, Pn.n_face));

   double Delta = 0.0000001; //a relax variable for floating point number comparision
  
  if(temp_dis>d_cvmax || temp_dis < Delta) return -1; //greater than largest possible separation, definitely not overlapping or checking the particle itself, i.e., m = n; this exclude the possible p2p contact
  else if(temp_dis<d_cfmin &&  temp_dis>=Delta) return -2; //smaller than the smallest possible separation, definitely overlapping


  //if in the intermediate distance, need to check with the separation axis theorem
  //need to do the complete check if overlapping: the most inefficient cases...

  int overlap_flag = -2; //-2 is the indicator for overlap
  double min_gap = 1000000.0; //to save the minimal gap...

  int face_flag = -1; //indicate whether it's a face contact or not; and the nature of the contact
  int face_normal_index; //just specify the index of the face which gives the smallest gap
                   //positive value for Pm, negative for Pn (just take -face_normal for Pn face index)

  //*******************************************************************
 

  //first, check the faces of Pn, i.e., loop over all faces of N and vertex of M...
  for(int i=0; i<Pn.n_face; i++) //the face normal
    {
      //the efficiency for face check can be further improved by 
      //first computing the inner product of face normal and centroid displacement
      //only when the inner product is positive that we need to check for projections...
      
      double temp_dist[Pm.n_vertex];
     
      for(int j=0; j<Pm.n_vertex; j++) //the vertex
	{
	  double temp_vertex[SD];	 

	  for(int k=0; k<SD; k++)
	    temp_vertex[k] = Pm.Vertex[j][k] + temp_VdisE[k];

	  temp_dist[j] = Comput.GetInnerProduct(Pn.Normal[i], temp_vertex);
	}

      //each face may have a different center2face distance, make sure they match
      if(GetMin(temp_dist, Pm.n_vertex)>Pn.d_cf[i]) //the particles do not overlap...
	{
	  overlap_flag = 0;
	  //return overlap_flag;
	  
	  face_flag = 1;
	  double temp_d = fabs(GetMin(temp_dist, Pm.n_vertex)-Pn.d_cf[i]);
	  if(temp_d < min_gap)
	    {
	      min_gap = temp_d;
	      face_normal_index = i;
	    }
	}
    }//end of loop for the faces...

 

  //second, check the faces of Pm  
  for(int i=0; i<Pm.n_face; i++)
   {
     double temp_dist[Pn.n_vertex];
   
      for(int j=0; j<Pn.n_vertex; j++) //the vertex
	{
	  double temp_vertex[SD];

	  for(int k=0; k<SD; k++)
	    temp_vertex[k] = Pn.Vertex[j][k] - temp_VdisE[k];

	  temp_dist[j] = Comput.GetInnerProduct(Pm.Normal[i], temp_vertex);
	}

      if(GetMin(temp_dist, Pn.n_vertex)>Pm.d_cf[i]) //the particles do not overlap...
	{
	  overlap_flag = 0;
	  //return overlap_flag;

	  face_flag = 1;
	  double temp_d = fabs(GetMin(temp_dist, Pn.n_vertex)-Pm.d_cf[i]);
	  if(temp_d < min_gap)
	    {
	      min_gap = temp_d;
	      face_normal_index = -i;
	    }
	}

   }//end of loop for the faces...
  //if a face contact, check it's nature...
 
  if(face_flag == 1)
    {
      int temp_n_vert = 0; //for the number of points...
      double ip_normal =10000.0; //this is to check whether f2f contact
                        //a more sensitive indicator than the num of projected points
    
      int dist_flag = 0; //flag indicating the read distance for possible each contact...
                      //this is important because although the projected distance could be small, 
                      //   the actuall separation could still be very large...
       
      if(face_normal_index>0)
	{
	  //for face of particle n
	  temp_n_vert = 0;
	  ip_normal =10000.0;
	
	  dist_flag = 0;


	  //a characterisitc length to determine whether really contact, we used the edge length here...
	  double temp_edge_len = sqrt((Pn.Vertex[Pn.Face[face_normal_index][0]][0]-Pn.Vertex[Pn.Face[face_normal_index][1]][0])*(Pn.Vertex[Pn.Face[face_normal_index][0]][0]-Pn.Vertex[Pn.Face[face_normal_index][1]][0])+(Pn.Vertex[Pn.Face[face_normal_index][0]][1]-Pn.Vertex[Pn.Face[face_normal_index][1]][1])*(Pn.Vertex[Pn.Face[face_normal_index][0]][1]-Pn.Vertex[Pn.Face[face_normal_index][1]][1])+(Pn.Vertex[Pn.Face[face_normal_index][0]][2]-Pn.Vertex[Pn.Face[face_normal_index][1]][2])*(Pn.Vertex[Pn.Face[face_normal_index][0]][2]-Pn.Vertex[Pn.Face[face_normal_index][1]][2]));

	  
	  //the center coordinates of Pn.Face[face_normal_index]
	  double Pn_face_coord[SD] = {0, 0, 0};
	  for(int k=0; k<SD; k++)
	    {
	      for(int j=0; j<Pn.vert_num_fspecie[Pn.face_type[face_normal_index]]; j++)
		Pn_face_coord[k] += Pn.Vertex[Pn.Face[face_normal_index][j]][k];

	      Pn_face_coord[k] = Pn_face_coord[k]/(double)Pn.vert_num_fspecie[Pn.face_type[face_normal_index]];
	    }
	  
	  double Pm_face_coord[SD] = {0, 0, 0};


	  //array saving the index of Pm vertices contribute to temp_n_vert
	  int Temp_Vert_Ind[Pm.n_vertex];


	  //for vertices of particle m...
	  double temp_dist[Pm.n_vertex];
     
	  for(int j=0; j<Pm.n_vertex; j++) 
	    {
	      double temp_vertex[SD];	 
	      
	      for(int k=0; k<SD; k++)
		temp_vertex[k] = Pm.Vertex[j][k] + temp_VdisE[k];
	      
	      temp_dist[j] = Comput.GetInnerProduct(Pn.Normal[face_normal_index], temp_vertex) - Pn.d_cf[face_normal_index];

	      if(fabs(temp_dist[j])<TOL*d0)
		{
		  Temp_Vert_Ind[temp_n_vert] = j;
		  temp_n_vert++;
		}
	    }

	  //for the face normal of particle m, check for f2f...
	  double temp_ip_normal;

	  for(int j=0; j<Pm.n_face; j++)
	    {
	      temp_ip_normal = Comput.GetInnerProduct(Pn.Normal[face_normal_index], Pm.Normal[j]);

	      if(fabs(temp_ip_normal) >1.0) cout<<"The normals are not properly normalized...!"<<endl;
		 
	      if(temp_ip_normal<ip_normal) 
		{ 
		  ip_normal = temp_ip_normal;

		  for(int k=0; k<SD; k++)
		    {
		      Pm_face_coord[k] = 0;
		      
		      for(int l=0; l<Pm.vert_num_fspecie[Pm.face_type[j]]; l++)
			Pm_face_coord[k] += Pm.Vertex[Pm.Face[j][l]][k];
		      
		      Pm_face_coord[k] = Pm_face_coord[k]/(double)Pm.vert_num_fspecie[Pm.face_type[j]];
		    }
		}
	      
	    }

	  //now check the nature of the contact...
	  //******************************************************
	  if(fabs(ip_normal)>normal_cutoff) 
	    {

	      double temp_dist = 0;
	      
	      //compute the face center2cetner distance
	      for(int k=0; k<SD; k++)
		temp_dist += (Pm_face_coord[k]+temp_VdisE[k]-Pn_face_coord[k])*(Pm_face_coord[k]+temp_VdisE[k]-Pn_face_coord[k]);
	      temp_dist = sqrt(temp_dist);

	      if(temp_dist < f2f_cutoff*temp_edge_len)
		{
		  face_flag = 1;
		  dist_flag = 1;
		  return face_flag;
		}	  
	    }
	  else if(temp_n_vert == 2) 
	    {
	      double temp_dist = 0;
	      
	      //compute the face center2cetner distance
	      for(int k=0; k<SD; k++)
		temp_dist += ((Pm.Vertex[Temp_Vert_Ind[0]][k]+Pm.Vertex[Temp_Vert_Ind[1]][k])/2.0+temp_VdisE[k]-Pn_face_coord[k])*((Pm.Vertex[Temp_Vert_Ind[0]][k]+Pm.Vertex[Temp_Vert_Ind[1]][k])/2.0+temp_VdisE[k]-Pn_face_coord[k]);
	      temp_dist = sqrt(temp_dist);
	      
	      if(temp_dist < e2f_cutoff*temp_edge_len)
		{
		  face_flag = 2;
		  dist_flag = 1;
		  return face_flag;
		}
	    }
	  else if(temp_n_vert == 1) 
	    {
	      double temp_dist = 0;
	      
	      //compute the face center2cetner distance
	      for(int k=0; k<SD; k++)
		temp_dist += (Pm.Vertex[Temp_Vert_Ind[0]][k]+temp_VdisE[k]-Pn_face_coord[k])*(Pm.Vertex[Temp_Vert_Ind[0]][k]+temp_VdisE[k]-Pn_face_coord[k]);
	      temp_dist = sqrt(temp_dist);
	      
	      if(temp_dist < v2f_cutoff*temp_edge_len)
		{ 
		  face_flag = 3;
		  dist_flag = 1;
		  return face_flag;
		}
	    }
	  else if(dist_flag == 0)
	    {
	      //this means the particles are still too far away from each other ...
	      face_flag =-1;
	      return -1; 
	    }
	  
	  //&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%

	}
      else
	{
	  //for face of particle m
	  temp_n_vert = 0;
	  ip_normal =10000.0;
	  dist_flag = 0;

	  face_normal_index = -face_normal_index;

	  //a characterisitc length to determine whether really contact, we used the edge length here...
	  double temp_edge_len = sqrt((Pm.Vertex[Pm.Face[face_normal_index][0]][0]-Pm.Vertex[Pm.Face[face_normal_index][1]][0])*(Pm.Vertex[Pm.Face[face_normal_index][0]][0]-Pm.Vertex[Pm.Face[face_normal_index][1]][0])+(Pm.Vertex[Pm.Face[face_normal_index][0]][1]-Pm.Vertex[Pm.Face[face_normal_index][1]][1])*(Pm.Vertex[Pm.Face[face_normal_index][0]][1]-Pm.Vertex[Pm.Face[face_normal_index][1]][1])+(Pm.Vertex[Pm.Face[face_normal_index][0]][2]-Pm.Vertex[Pm.Face[face_normal_index][1]][2])*(Pm.Vertex[Pm.Face[face_normal_index][0]][2]-Pm.Vertex[Pm.Face[face_normal_index][1]][2]));

	  
	  //the center coordinates of Pm.Face[face_normal_index]
	  double Pm_face_coord[SD] = {0, 0, 0};
	  for(int k=0; k<SD; k++)
	    {
	      for(int j=0; j<Pm.vert_num_fspecie[Pm.face_type[face_normal_index]]; j++)
		Pm_face_coord[k] += Pm.Vertex[Pm.Face[face_normal_index][j]][k];

	      Pm_face_coord[k] = Pm_face_coord[k]/(double)Pm.vert_num_fspecie[Pm.face_type[face_normal_index]];
	    }
	  
	  double Pn_face_coord[SD] = {0, 0, 0};


	  //array saving the index of Pn vertices contribute to temp_n_vert
	  int Temp_Vert_Ind[Pn.n_vertex];

	  
	  //for vertex of particel n...
	  double temp_dist[Pn.n_vertex];
     
	  for(int j=0; j<Pn.n_vertex; j++) 
	    {
	      double temp_vertex[SD];	 
	      
	      for(int k=0; k<SD; k++)
		temp_vertex[k] = Pn.Vertex[j][k] - temp_VdisE[k];
	      
	      temp_dist[j] = Comput.GetInnerProduct(Pm.Normal[face_normal_index], temp_vertex) - Pm.d_cf[face_normal_index];

	      if(fabs(temp_dist[j])<TOL*d0)
		{
		  Temp_Vert_Ind[temp_n_vert] = j;
		  temp_n_vert++;
		}
	    }

	  //for the face normal of particle n, check for f2f...
	  double temp_ip_normal;

	  for(int j=0; j<Pn.n_face; j++)
	    {
	      temp_ip_normal = Comput.GetInnerProduct(Pm.Normal[face_normal_index], Pn.Normal[j]);

	      if(fabs(temp_ip_normal) >1.0) cout<<"The normals are not properly normalized...!"<<endl;
		 
	      if(temp_ip_normal<ip_normal) 
		{
		  ip_normal = temp_ip_normal;

		  for(int k=0; k<SD; k++)
		    {
		      Pn_face_coord[k] = 0;
		      
		      for(int l=0; l<Pn.vert_num_fspecie[Pn.face_type[j]]; l++)
			Pn_face_coord[k] += Pn.Vertex[Pn.Face[j][l]][k];
		      
		      Pn_face_coord[k] = Pn_face_coord[k]/(double)Pn.vert_num_fspecie[Pn.face_type[j]];
		    }
		}
		 
	    }

	  //now check the nature of the contact...
	  if(fabs(ip_normal)>normal_cutoff) 
	    {
	      
	      double temp_dist = 0;
	      
	      //compute the face center2cetner distance
	      for(int k=0; k<SD; k++)
		temp_dist += (Pm_face_coord[k]+temp_VdisE[k]-Pn_face_coord[k])*(Pm_face_coord[k]+temp_VdisE[k]-Pn_face_coord[k]);
	      temp_dist = sqrt(temp_dist);
	      
	      if(temp_dist < f2f_cutoff*temp_edge_len)
		{
		  face_flag = 1;
		  dist_flag = 1;
		  return face_flag;
		}	  
	    }
	  else if(temp_n_vert == 2) 
	    {
	      double temp_dist = 0;
	      
	      //compute the face center2cetner distance
	      for(int k=0; k<SD; k++)
		temp_dist += ((Pn.Vertex[Temp_Vert_Ind[0]][k]+Pn.Vertex[Temp_Vert_Ind[1]][k])/2.0-temp_VdisE[k]-Pm_face_coord[k])*((Pn.Vertex[Temp_Vert_Ind[0]][k]+Pn.Vertex[Temp_Vert_Ind[1]][k])/2.0-temp_VdisE[k]-Pm_face_coord[k]);
	      temp_dist = sqrt(temp_dist);
	      
	      if(temp_dist < e2f_cutoff*temp_edge_len)
		{
		  face_flag = 2;
		  dist_flag = 1;
		  return face_flag;
		}
	    }
	  else if(temp_n_vert == 1) 
	    {
	      double temp_dist = 0;
	      
	      //compute the face center2cetner distance
	      for(int k=0; k<SD; k++)
		temp_dist += (Pn.Vertex[Temp_Vert_Ind[0]][k]-temp_VdisE[k]-Pm_face_coord[k])*(Pn.Vertex[Temp_Vert_Ind[0]][k]-temp_VdisE[k]-Pm_face_coord[k]);
	      temp_dist = sqrt(temp_dist);
	      
	      if(temp_dist < v2f_cutoff*temp_edge_len)
		{ 
		  face_flag = 3;
		  dist_flag = 1;
		  return face_flag;
		}
	    }
	  else if(dist_flag == 0)
	    {
	      //this means the particles are still too far away from each other ...
	      face_flag =-1;
	      return -1; 
	    }
	}
      
    }
  

  //*************************************************************
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //next we need to check all pair of edges to find the separating axis...

  min_gap = 1000000.0; //reset the minimal gap...
  //int edge_dist_flag = 0;
  //also check the distance between the edges

  double edgeM[Pm.n_edge][SD], edgeN[Pn.n_edge][SD];
  double eM_center[Pm.n_edge][SD], eN_center[Pn.n_edge][SD];
  double eM_len[Pm.n_edge], eN_len[Pn.n_edge];

  for(int j=0; j<SD; j++)
    {
      for(int i=0; i<Pm.n_edge; i++)
	{
	  edgeM[i][j] = Pm.Vertex[Pm.Edge[i][0]][j] - Pm.Vertex[Pm.Edge[i][1]][j];
	  eM_center[i][j] = (Pm.Vertex[Pm.Edge[i][0]][j] + Pm.Vertex[Pm.Edge[i][1]][j])/2.0;
	}
      
      for(int i=0; i<Pn.n_edge; i++)
	{
	  edgeN[i][j] = Pn.Vertex[Pn.Edge[i][0]][j] - Pn.Vertex[Pn.Edge[i][1]][j];
	  eN_center[i][j] = (Pn.Vertex[Pn.Edge[i][0]][j] + Pn.Vertex[Pn.Edge[i][1]][j])/2.0;
	}
    }

  for(int i=0; i<Pm.n_edge; i++)
    {
      eM_len[i] = 0;
      for(int j=0; j<SD; j++)
	eM_len[i] += edgeM[i][j]*edgeM[i][j];
      eM_len[i] = sqrt(eM_len[i]);
    }
  for(int i=0; i<Pn.n_edge; i++)
    {
      eN_len[i] = 0;
      for(int j=0; j<SD; j++)
	eN_len[i] += edgeN[i][j]*edgeN[i][j];
      eN_len[i] = sqrt(eN_len[i]);
    }


  for(int i=0; i<Pm.n_edge; i++)
    for(int j=0; j<Pn.n_edge; j++)
      {
	double temp_axis[SD];
	
	Comput.GetCrossProduct(edgeM[i], edgeN[j], temp_axis);

	//dealing with the parallel case
	if(Comput.GetLength(temp_axis)==0)
	  {
	    double temp_vect1[SD], temp_vect2[SD];

	    //either end points of the edge is fine
	    for(int v=0; v<SD; v++)
	      temp_vect1[v] = Pm.Vertex[Pm.Edge[i][0]][v]+temp_VdisE[v]; 
	    for(int v=0; v<SD; v++)
	      temp_vect2[v] = Pn.Vertex[Pn.Edge[j][0]][v]; 

	    for(int v=0; v<SD; v++)
	      temp_vect1[v] = temp_vect1[v] - temp_vect2[v];

	    Comput.GetPLine(edgeM[i], temp_vect1, temp_axis);
	  }


	//the separating axis has been obtained, now check the vertices...

	double temp_disM[Pm.n_vertex], temp_disN[Pn.n_vertex];

	for(int k=0; k<Pm.n_vertex; k++)
	  {
	    temp_disM[k] = 0;
    	  }
	for(int k=0; k<Pn.n_vertex; k++)
	  {
	    temp_disN[k] = 0;
    	  }

	for(int u=0; u<Pm.n_vertex; u++)
	  for(int k=0; k<SD; k++)
	    temp_disM[u] += temp_axis[k]*(Pm.Vertex[u][k] + temp_VdisE[k]);
	
	for(int u=0; u<Pn.n_vertex; u++)
	  for(int k=0; k<SD; k++)
	    temp_disN[u] += temp_axis[k]*Pn.Vertex[u][k];


	//now we comput the inner product of two vectors connecting the ends of an edge to the center of another edge
	//if this is >0, meaning the edges are not crossing
	double ec_inner = 0;
	//this is the distance between the edge cetners
	double ec_dist = 0;
	

	for(int k=0; k<SD; k++)
	  {
	    ec_inner += (Pm.Vertex[Pm.Edge[i][0]][k]+temp_VdisE[k]-eN_center[j][k])*(Pm.Vertex[Pm.Edge[i][1]][k]+temp_VdisE[k]-eN_center[j][k]);
	    ec_dist += (eM_center[i][k]+temp_VdisE[k]-eN_center[j][k])*(eM_center[i][k]+temp_VdisE[k]-eN_center[j][k]);
	  }
	ec_dist = sqrt(ec_dist);

	    
	if(GetMax(temp_disM, Pm.n_vertex)<GetMin(temp_disN, Pn.n_vertex) || GetMax(temp_disN, Pn.n_vertex)<GetMin(temp_disM, Pm.n_vertex))
	  {
	    overlap_flag = 0;
	    
	    //return overlap_flag;

	    face_flag = 4;

	    double temp_d1 = fabs(GetMax(temp_disM, Pm.n_vertex)-GetMin(temp_disN, Pn.n_vertex));
	    if(temp_d1<min_gap && ec_inner<0 && ec_dist<edge_cutoff*(eM_len[i]+eN_len[j])/2.0) 
	      {
		//this additional conditions makes sure the distance between particles are not large
		min_gap = temp_d1;
	      }

	    double temp_d2 = fabs(GetMax(temp_disN, Pn.n_vertex)-GetMin(temp_disM, Pm.n_vertex));
	    if(temp_d2<min_gap && ec_inner<0 && ec_dist<edge_cutoff*(eM_len[i]+eN_len[j])/2.0) 
	      {
		min_gap = temp_d2;
	      }
	  }
      }//end of the loop for the edge pairs...
  
  //meaning there is possible e2e contact
  if(face_flag == 4)
    {
      if(min_gap<TOL*d0) return face_flag;
      else return -1;
    }

  //if we get here, then the M and N overlap...
  
  return -2; //-2 for overlapping

 

}

void Contact::GlobalGetContact()
{
  //loop over all particle pairs 
  for(int m=0; m<N; m++)
    for(int n=m; n<N; n++)
      {
	//cout<<"In GOC "<<m<<" "<<n<<endl;
	
	for(int i=-1; i<=1; i++)
	  for(int j=-1; j<=1; j++)
	    for(int k=-1; k<=1; k++)
	      {
		int temp_ca = GetContact(PackPoly[n], PackPoly[m], CenterL[n], CenterL[m], i, j, k);

		if(temp_ca>0)
		  Contact_Array[m][n] = Contact_Array[n][m] = temp_ca;
	      }
      }
  

  //now check the value of contact array and get statistics 
  for(int i=0; i<N; i++)
    {
      int temp_ca;

      for(int j=0; j<N; j++)
	{
	  temp_ca = Contact_Array[i][j];

	  if(temp_ca == 1)
	    {
	      NP_f2f[i]++;
	      n_f2f++;
	    }
	  else if(temp_ca == 2)
	    {
	      NP_e2f[i]++;
	      n_e2f++;
	    }
	  else if(temp_ca == 3)
	    {
	      NP_p2f[i]++;
	      n_p2f++;
	    }
	  else if(temp_ca == 4)
	    {
	      NP_e2e[i]++;
	      n_e2e++;
	    }
	}
    }

  //we double count the global contact
  n_f2f = int(n_f2f);
  n_e2f = int(n_e2f);
  n_p2f = int(n_p2f);
  n_e2e = int(n_e2e);
}

void Contact::Print_Stat()
{
  
  ofstream fout;

  fout.open("Contact_Neighbors.txt");
  fout<<" index of particle, and the index of its contacting neighbors "<<endl;
  for(int i=0; i<N; i++)
    {
      fout<<i<<"\t\t";
      for(int j=0; j<N; j++)
	{
	  if(Contact_Array[i][j]>0)
	    fout<<j<<" ";
	}
      fout<<endl;
    }
  fout.close();

  fout.open("Contact_Num.txt");
  fout<<" index of particle, n_f2f    n_e2f    n_p2f    n_e2e "<<endl;
  for(int i=0; i<N; i++)
    {
      fout<<i<<"\t"<<NP_f2f[i]<<" "<<NP_e2f[i]<<" "<<NP_p2f[i]<<" "<<NP_e2e[i]<<endl;
    }
  fout.close();


  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
  cout<<"Average contact n_ave = n_f2f + n_e2f + n_p2f + n_e2e = "<<(double)n_f2f/(double)N<<" + "<<(double)n_e2f/(double)N<<" + "<<(double)n_p2f/(double)N<<" + "<<(double)n_e2e/(double)N<<" = "<<(double)(n_f2f+n_e2f+n_p2f+n_e2e)/(double)N<<endl;

  cout<<"Total D.O.F = "<<N*6<<endl;
  cout<<"Total Constraints = "<<(n_f2f*3+n_e2f*2+n_p2f+n_e2e)/2.0<<endl;
  
  
}


