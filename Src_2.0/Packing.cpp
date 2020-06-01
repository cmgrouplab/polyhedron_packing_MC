


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


void Packing::Rescale(double relax_ratio) //if packing containing overlapping pairs, rescale the packing
{
  for(int i=0; i<SD; i++)
    for(int j=0; j<SD; j++)
      Lambda[i][j] = relax_ratio*Lambda[i][j];
}

void Packing::Check_Iconfig(double temp_option, double relax_ratio)
{
  //cout<<"In Check_Iconfig 1"<<endl;

  int temp_flag = GlobalOverlapCheck(1);
  
  if(temp_flag == 1)
    {
  
      if(temp_option == 0)
	{
	  while(GlobalOverlapCheck(1)==1)
	    {
	      //printf("The Initial Configuration Contains Overlapping Pairs!\n");
	      Rescale(relax_ratio);
	      //exit(1);
	    }
	}
      else if(temp_option == 1)
	{ 
	  exit(1);
	}
    }

}

Packing::Packing(char* temp_name)
{
  FILE* fp;

  if((fp=fopen("Iconfig.txt","r"))==NULL)
    {
      printf("Cannot open file Iconfig.txt to initialize Packing!\n");
      exit(1);
    }
  else //read the initial configuration...
    {
      double val;
      fscanf(fp, "%d", &N); // read in the number of particles...
      
      fscanf(fp, "%lf", &d0); // read in the length...
      //for monodisperse packings, d0 is the edge length of the shape
      //for polydisperse or random shapes, should be the maximal length for all shapes.

      Lambda = new double*[SD];
      for(int i=0; i<SD; i++)
	Lambda[i] = new double[SD];
      
      for(int i=0; i<SD; i++) // read in the lattice vectors...
	for(int j=0; j<SD; j++)
	  {
	    fscanf(fp, "%lf", &val);
	    Lambda[j][i] = val;
	  }

        //now initalize the variables

      CenterE = new double*[N];
      CenterL = new double*[N];

      PackPoly = new Polyhedron[N]; //new [] cannot call constructors with parameters
      
      for(int i=0; i<N; i++)
	{
	  CenterE[i] = new double[SD];	  
	  CenterL[i] = new double[SD];

	  PackPoly[i].PolyConstr(temp_name);
	  
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

	  
	   
	  //read in vertex with O = 0; and to make sure, shift the vertices again for each particle
	  
	  for(int v=0; v<PackPoly[i].n_vertex; v++)
	    for(int j=0; j<SD; j++)
	      {
		fscanf(fp, "%lf", &val);
		PackPoly[i].Vertex[v][j] = val;
	      }

	 
	  PackPoly[i].ShiftVert(); //shift vertices to make the centroid at orgin

	  //PackPoly[i].GetNormal_All(); //get the face normals
	  //no need to compute all the face normals, only compute those used for overlap check

	  PackPoly[i].GetCentDist(); //get d_cf and d_cv

	  

	  //get the maximum characteristic length of the particle... 
	  double temp_d = 2*GetMax(PackPoly[i].d_cv, PackPoly[i].n_vertex);
	  if(temp_d>d0) d0 = temp_d;

	}//finishing loop over all particles
      
      
      
      fclose(fp);
    }

  
  //now we compute the total volume of all particles, which is a constant of the packing

  GetVParticle();

}

void Packing::GetGlobalPosition()
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


double Packing::MinDis(double CenterLm[SD], double CenterLn[SD], double temp_disE[SD], int indexI, int indexJ, int indexK)
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


double Packing::MinDisII(int m, int n, double temp_disE[SD], int &indexI, int &indexJ, int &indexK)
{
  double tempLd[SD];
  double tempEd[SD];
  double min_dist = 10000000000.0;
  double temp_dist;

  for(int i=-1; i<=1; i++)
    for(int j=-1; j<=1; j++)
      for(int k=-1; k<=1; k++)
	{

	  tempLd[0] = (CenterL[m][0]-CenterL[n][0])+i;
	  tempLd[1] = (CenterL[m][1]-CenterL[n][1])+j;
	  tempLd[2] = (CenterL[m][2]-CenterL[n][2])+k;

	  for(int q=0; q<SD; q++)
	    tempEd[q] = 0;
	  
	  for(int q=0; q<SD; q++)
	    for(int h=0; h<SD; h++)
	      tempEd[q] += Lambda[q][h]*tempLd[h];
	  
	  temp_dist = tempEd[0]*tempEd[0]+tempEd[1]*tempEd[1]+tempEd[2]*tempEd[2]; 

	  if(temp_dist < min_dist)
	    {
	      min_dist = temp_dist;

	      for(int q=0; q<SD; q++)
		temp_disE[q] = tempEd[q];

	      indexI = i; indexJ = j; indexK = k;
	    }
	}

  return sqrt(min_dist);
  
}


void Packing::GetNNL(double cut_dist)
{
  NNL = new int**[N];

  for(int i=0; i<N; i++)
    {
      NNL[i] = new int*[N];

      for(int j=0; j<N; j++)
	NNL[i][j] = new int [4];
    }

  NNL_counter = new int [N];

  int temp_ind1, temp_ind2, temp_ind3;
  double temp_dist[SD];
  double d_cvmax;
  
  //loop for all particles for near-neighbors...
  for(int i=0; i<N; i++)
    {
      NNL_counter[i] = 0;
      
      for(int j=0; j<N; j++)
	{
	  d_cvmax = (GetMax(PackPoly[i].d_cv, PackPoly[i].n_vertex)+GetMax(PackPoly[j].d_cv, PackPoly[j].n_vertex));

	  if(MinDisII(i, j, temp_dist, temp_ind1, temp_ind2, temp_ind3) < cut_dist*d_cvmax)
	    {
	      NNL[i][NNL_counter[i]][0] = j;

	      NNL[i][NNL_counter[i]][1] = temp_ind1;
	      NNL[i][NNL_counter[i]][2] = temp_ind2;
	      NNL[i][NNL_counter[i]][3] = temp_ind3;

	      NNL_counter[i]++;
	    }
	}

      //cout<<NNL_counter[i]<<endl;
    }
}


void Packing::UpdateNNL(double cut_dist)
{
  int temp_ind1, temp_ind2, temp_ind3;
  double temp_dist[SD];
  double d_cvmax;
  
  //loop for all particles for near-neighbors...
  for(int i=0; i<N; i++)
    {
      NNL_counter[i] = 0;
      
      for(int j=0; j<N; j++)
	{
	  d_cvmax = (GetMax(PackPoly[i].d_cv, PackPoly[i].n_vertex)+GetMax(PackPoly[j].d_cv, PackPoly[j].n_vertex));

	  if(MinDisII(i, j, temp_dist, temp_ind1, temp_ind2, temp_ind3) < cut_dist*d_cvmax)
	    {
	      NNL[i][NNL_counter[i]][0] = j;

	      NNL[i][NNL_counter[i]][1] = temp_ind1;
	      NNL[i][NNL_counter[i]][2] = temp_ind2;
	      NNL[i][NNL_counter[i]][3] = temp_ind3;
	      
	      NNL_counter[i]++;
	    }
	}

      //cout<<NNL_counter[i]<<endl;
    }
}


double Packing::GetMin(double Arr[], int Size)
{
  double temp_min = 1000000000.0;

  for(int i=0; i<Size; i++)
    {
      if(Arr[i]<temp_min) temp_min = Arr[i];
    }

  return temp_min;
}

double Packing::GetMax(double Arr[], int Size)
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

int Packing::CheckOverlap(Polyhedron &Pm, Polyhedron &Pn, double CenterLm[SD], double CenterLn[SD], int indexI, int indexJ, int indexK)
{
  double temp_VdisE[SD]; //the relative and Ecludean separation distance... the vector from n to m
  
  //first get the center-to-center distance of the particles...
  double temp_dis = MinDis(CenterLm, CenterLn, temp_VdisE, indexI, indexJ, indexK);
  //recall n is fixed in the central box, the images of m are checked...

  //the mininmal and maximal possible separations
  double d_cvmax = (GetMax(Pm.d_cv, Pm.n_vertex)+GetMax(Pn.d_cv, Pn.n_vertex));
  double d_cfmin = (GetMin(Pm.d_cf, Pm.n_face)+GetMin(Pn.d_cf, Pn.n_face));

   double Delta = 0.0000001; //a relax variable for floating point number comparision
  
  if(temp_dis>d_cvmax || temp_dis < Delta) return 0; //greater than largest possible separation, definitely not overlapping or checking the particle itself, i.e., m = n
  else if(temp_dis<d_cfmin &&  temp_dis>=Delta) return 1; //smaller than the smallest possible separation, definitely overlapping


  //if in the intermediate distance, need to check with the separation axis theorem
  //need to do the complete check if overlapping: the most inefficient cases...


  double temp_VdistEII[SD]; //the vector from m to n
  for(int k=0; k<SD; k++)
    temp_VdistEII[k] = -temp_VdisE[k];

  //*******************************************************************
  int overlap_flag = 1;

  //first, check the faces of Pn, i.e., loop over all faces of N and vertex of M...
  for(int i=0; i<Pn.n_face; i++) //the face normal
    {
      //the efficiency for face check can be further improved by 
      //first computing the inner product of face normal and centroid displacement
      //only when the inner product is positive that we need to check for projections...

      //compute the normal 
      Pn.GetNormal_F(i);
      
      if(Comput.GetInnerProduct(Pn.Normal[i], temp_VdisE) <= 0.0)
	continue;

      
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
	  return overlap_flag;
	}


    }//end of loop for the faces...

 

  //second, check the faces of Pm  
  for(int i=0; i<Pm.n_face; i++)
   {

     //similar to the previous case, first check the direction of the normal

     //compute the normal to be checked
     Pm.GetNormal_F(i);
   

     if(Comput.GetInnerProduct(Pm.Normal[i], temp_VdistEII) <= 0.0)
       continue;


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
	  return overlap_flag;
	}

   }//end of loop for the faces...


  

  //*************************************************************
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //next we need to check all pair of edges to find the separating axis...
  //this can also been improved by first check ''edge normal'', i.e., by 
  // by computing the inner product of the edge normal and the centroid vector 
  // the edge normal is better to computed as an intrinsic geometric property of the particle

  double edgeM[Pm.n_edge][SD], edgeN[Pn.n_edge][SD];
  double edge_ctM[Pm.n_edge][SD], edge_ctN[Pn.n_edge][SD];
  //the center of the edge, for fast check, if point to the opposite direction from
  // the vector connecting particle centers, directly skip the pair

  for(int j=0; j<SD; j++)
    {
      for(int i=0; i<Pm.n_edge; i++)
	{
	  edgeM[i][j] = Pm.Vertex[Pm.Edge[i][0]][j] - Pm.Vertex[Pm.Edge[i][1]][j];
	  edge_ctM[i][j] = (Pm.Vertex[Pm.Edge[i][0]][j] + Pm.Vertex[Pm.Edge[i][1]][j])/2.0;
	}
      
      for(int i=0; i<Pn.n_edge; i++)
	{
	  edgeN[i][j] = Pn.Vertex[Pn.Edge[i][0]][j] - Pn.Vertex[Pn.Edge[i][1]][j];
	  edge_ctN[i][j] = (Pn.Vertex[Pn.Edge[i][0]][j] + Pn.Vertex[Pn.Edge[i][1]][j])/2.0;
	}
    }


  for(int i=0; i<Pm.n_edge; i++)
    for(int j=0; j<Pn.n_edge; j++)
      {

	//first check the edge center vector and particle centroid vector
	if(Comput.GetInnerProduct(edge_ctN[i], temp_VdisE)<=0.0 &&  Comput.GetInnerProduct(edge_ctM[i], temp_VdistEII)<=0.0)
	  continue;

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
	    
	if(GetMax(temp_disM, Pm.n_vertex)<GetMin(temp_disN, Pn.n_vertex) || GetMax(temp_disN, Pn.n_vertex)<GetMin(temp_disM, Pm.n_vertex))
	  {
	    overlap_flag = 0;
	    
	    return overlap_flag;
	  }
      }//end of the loop for the edge pairs...
  
  
  //if we get here, then the M and N overlap...
  
  return overlap_flag;

}


void Packing::GetSEP_List(char* temp_name)
{
  
  //first, initialize the data
  SEP_List = new int**[N];

  for(int i=0; i<N; i++)
    {
      SEP_List[i] = new int*[N];

      for(int j=0; j<N; j++)
	{
	  SEP_List[i][j] = new int [3];

	  for(int k=0; k<3; k++)
	    SEP_List[i][j][k] = -1;
	}
    }

  
  //now set up the intial SEP_List
  //loop over all particle for overlap check...

  
  Polyhedron Pm(temp_name);
  Polyhedron Pn(temp_name);
  int indexI, indexJ, indexK;

  for(int m=0; m<N; m++)
    for(int n=0; n<N; n++)
      {
  
	Pm.CopyPoly(PackPoly[m]);
	Pn.CopyPoly(PackPoly[n]);

	double temp_VdisE[SD]; //the relative and Ecludean separation distance... the vector from n to m
  
	//first get the center-to-center distance of the particles...
	double temp_dis = MinDisII(m, n, temp_VdisE, indexI, indexJ, indexK);
	//recall n is fixed in the central box, the images of m are checked...

	//the mininmal and maximal possible separations
	double d_cvmax = (GetMax(Pm.d_cv, Pm.n_vertex)+GetMax(Pn.d_cv, Pn.n_vertex));
	double d_cfmin = (GetMin(Pm.d_cf, Pm.n_face)+GetMin(Pn.d_cf, Pn.n_face));

	double Delta = 0.0000001; //a relax variable for floating point number comparision
  
	if(temp_dis>d_cvmax || temp_dis < Delta)
	  {
	    SEP_List[m][n][0] = -1;
	    SEP_List[m][n][1] = -1;
	    SEP_List[m][n][2] = -1;
	    continue; //greater than largest possible separation, definitely not overlapping or checking the particle itself, i.e., m = n
	  }
	else if(temp_dis<d_cfmin &&  temp_dis>=Delta) 
	  {
	    printf("Unexpected overlap is detected when setting up SEP_List! Recheck! \n");
	    exit(1);
	    //smaller than the smallest possible separation, definitely overlapping
	  }


	//if in the intermediate distance, need to check with the separation axis theorem
	//need to do the complete check if overlapping: the most inefficient cases...


	double temp_VdistEII[SD]; //the vector from m to n
	for(int k=0; k<SD; k++)
	  temp_VdistEII[k] = -temp_VdisE[k];

	//*******************************************************************
	int overlap_flag = 1;

	//first, check the faces of Pn, i.e., loop over all faces of N and vertex of M...
	for(int i=0; i<Pn.n_face; i++) //the face normal
	  {
	    //the efficiency for face check can be further improved by 
	    //first computing the inner product of face normal and centroid displacement
	    //only when the inner product is positive that we need to check for projections...

	    //compute the normal 
	    Pn.GetNormal_F(i);
      
	    if(Comput.GetInnerProduct(Pn.Normal[i], temp_VdisE) <= 0.0)
	      continue;

      
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
		
		SEP_List[m][n][0] = 0;
		SEP_List[m][n][1] = 1;
		SEP_List[m][n][2] = i;
		
		goto L1;
	      }
	    

	  }//end of loop for the faces...

 

	//second, check the faces of Pm  
	for(int i=0; i<Pm.n_face; i++)
	  {

	    //similar to the previous case, first check the direction of the normal

	    //compute the normal to be checked
	    Pm.GetNormal_F(i);
   

	    if(Comput.GetInnerProduct(Pm.Normal[i], temp_VdistEII) <= 0.0)
	      continue;


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
	     
		SEP_List[m][n][0] = 0;
		SEP_List[m][n][1] = 0;
		SEP_List[m][n][2] = i;
		
		goto L1;
	      }
	    
	  }//end of loop for the faces...


	//*************************************************************
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//next we need to check all pair of edges to find the separating axis...
	//this can also been improved by first check ''edge normal'', i.e., by 
	// by computing the inner product of the edge normal and the centroid vector 
	// the edge normal is better to computed as an intrinsic geometric property of the particle
	
	double edgeM[Pm.n_edge][SD], edgeN[Pn.n_edge][SD];
	double edge_ctM[Pm.n_edge][SD], edge_ctN[Pn.n_edge][SD];
	//the center of the edge, for fast check, if point to the opposite direction from
	// the vector connecting particle centers, directly skip the pair
	
	for(int j=0; j<SD; j++)
	  {
	    for(int i=0; i<Pm.n_edge; i++)
	      {
		edgeM[i][j] = Pm.Vertex[Pm.Edge[i][0]][j] - Pm.Vertex[Pm.Edge[i][1]][j];
		edge_ctM[i][j] = (Pm.Vertex[Pm.Edge[i][0]][j] + Pm.Vertex[Pm.Edge[i][1]][j])/2.0;
	      }
	    
	    for(int i=0; i<Pn.n_edge; i++)
	      {
		edgeN[i][j] = Pn.Vertex[Pn.Edge[i][0]][j] - Pn.Vertex[Pn.Edge[i][1]][j];
		edge_ctN[i][j] = (Pn.Vertex[Pn.Edge[i][0]][j] + Pn.Vertex[Pn.Edge[i][1]][j])/2.0;
	      }
	  }
	
	
	for(int i=0; i<Pm.n_edge; i++)
	  for(int j=0; j<Pn.n_edge; j++)
	    {
	      
	      //first check the edge center vector and particle centroid vector
	      if(Comput.GetInnerProduct(edge_ctN[i], temp_VdisE)<=0.0 &&  Comput.GetInnerProduct(edge_ctM[i], temp_VdistEII)<=0.0)
		continue;
	      
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
	      
	      if(GetMax(temp_disM, Pm.n_vertex)<GetMin(temp_disN, Pn.n_vertex) || GetMax(temp_disN, Pn.n_vertex)<GetMin(temp_disM, Pm.n_vertex))
		{
		  overlap_flag = 0;
		  
		  SEP_List[m][n][0] = 1;
		  SEP_List[m][n][1] = i;
		  SEP_List[m][n][2] = j;

		  goto L1;
		  

		}
	    }//end of the loop for the edge pairs...
	
	
	//if we get here, then the M and N overlap...
       
	
      L1: if(overlap_flag == 1)
	  {
	    printf("Unexpected overlap is detected when setting up SEP_List! Recheck! \n");
	    exit(1);
	  }
      }
}



//this is the overlap check using SEP_List, the last two are the index of the particle
//if the particles where separated from faces, we restrict that all the following moves has to maintain this separation
//if the particles were separated by edges, it can be violated....

int Packing::CheckOverlap(Polyhedron &Pm, Polyhedron &Pn, double CenterLm[SD], double CenterLn[SD], int indexI, int indexJ, int indexK, int m, int n)
{

  int overlap_flag = 1;

  double temp_VdisE[SD]; //the relative and Ecludean separation distance... the vector from n to m
  
  //first get the center-to-center distance of the particles...
  double temp_dis = MinDis(CenterLm, CenterLn, temp_VdisE, indexI, indexJ, indexK);
  //recall n is fixed in the central box, the images of m are checked...

  double temp_VdistEII[SD]; //the vector from m to n
  for(int k=0; k<SD; k++)
    temp_VdistEII[k] = -temp_VdisE[k];
  
  if(SEP_List[m][n][0] == 0) //face-separate
    {
      int face_index;
      
      if(SEP_List[m][n][1] == 0) //face of m
	{
	  face_index = SEP_List[m][n][2];

	  //compute the normal to be checked
	  Pm.GetNormal_F(face_index);
	  	  
	  double temp_dist[Pn.n_vertex];
	  
	  for(int j=0; j<Pn.n_vertex; j++) //the vertex
	    {
	      double temp_vertex[SD];
	      
	      for(int k=0; k<SD; k++)
		temp_vertex[k] = Pn.Vertex[j][k] - temp_VdisE[k];
	      
	      temp_dist[j] = Comput.GetInnerProduct(Pm.Normal[face_index], temp_vertex);
	    }
	  
	  if(GetMin(temp_dist, Pn.n_vertex)>Pm.d_cf[face_index]) //the particles do not overlap...
	    {
	      overlap_flag = 0;
	      return overlap_flag;
	    }
	  //else return 1;
	}
      else //face of n
	{
	  face_index = SEP_List[m][n][2];

	  Pn.GetNormal_F(face_index);
	 
	  double temp_dist[Pm.n_vertex];
     
	  for(int j=0; j<Pm.n_vertex; j++) //the vertex
	    {
	      double temp_vertex[SD];	 
	      
	      for(int k=0; k<SD; k++)
		temp_vertex[k] = Pm.Vertex[j][k] + temp_VdisE[k];

	      temp_dist[j] = Comput.GetInnerProduct(Pn.Normal[face_index], temp_vertex);
	    }
	  
	  //each face may have a different center2face distance, make sure they match
	  if(GetMin(temp_dist, Pm.n_vertex)>Pn.d_cf[face_index]) //the particles do not overlap...
	    {
	      overlap_flag = 0;
	      return overlap_flag;
	    }
	  //else return 1;
	}
      
    }
  else if(SEP_List[m][n][0] == 1) //edge-separate
    {
      int edge_indm = SEP_List[m][n][1];
      int edge_indn = SEP_List[m][n][2];

      double edge_coordM[SD];
      double edge_coordN[SD];

      for(int j=0; j<SD; j++)
	{	 
	  edge_coordM[j] = Pm.Vertex[Pm.Edge[edge_indm][0]][j] - Pm.Vertex[Pm.Edge[edge_indm][1]][j];
	  edge_coordN[j] = Pn.Vertex[Pn.Edge[edge_indn][0]][j] - Pn.Vertex[Pn.Edge[edge_indn][1]][j];
	}
	

      double temp_axis[SD];
	      
      Comput.GetCrossProduct(edge_coordM, edge_coordN, temp_axis);
	      
      //dealing with the parallel case
      if(Comput.GetLength(temp_axis)==0)
	{
	  double temp_vect1[SD], temp_vect2[SD];
	  
	  //either end points of the edge is fine
	  for(int v=0; v<SD; v++)
	    temp_vect1[v] = Pm.Vertex[Pm.Edge[edge_indm][0]][v]+temp_VdisE[v]; 
	  for(int v=0; v<SD; v++)
	    temp_vect2[v] = Pn.Vertex[Pn.Edge[edge_indn][0]][v]; 
	  
	  for(int v=0; v<SD; v++)
	    temp_vect1[v] = temp_vect1[v] - temp_vect2[v];
	  
	  Comput.GetPLine(edge_coordM, temp_vect1, temp_axis);
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
	      
      if(GetMax(temp_disM, Pm.n_vertex)<GetMin(temp_disN, Pn.n_vertex) || GetMax(temp_disN, Pn.n_vertex)<GetMin(temp_disM, Pm.n_vertex))
	{
	  overlap_flag = 0;
		  
	  return overlap_flag;
	  
	}
      //else return 1; 
      
      
    }
  //if neither is true, it means need to check overlap from scratch...
  

  
  //double temp_VdisE[SD]; //the relative and Ecludean separation distance... the vector from n to m
  
  //first get the center-to-center distance of the particles...
  //double temp_dis = MinDis(CenterLm, CenterLn, temp_VdisE, indexI, indexJ, indexK);
  //recall n is fixed in the central box, the images of m are checked...

  //the mininmal and maximal possible separations
  double d_cvmax = (GetMax(Pm.d_cv, Pm.n_vertex)+GetMax(Pn.d_cv, Pn.n_vertex));
  double d_cfmin = (GetMin(Pm.d_cf, Pm.n_face)+GetMin(Pn.d_cf, Pn.n_face));

   double Delta = 0.0000001; //a relax variable for floating point number comparision
  
  if(temp_dis>d_cvmax || temp_dis < Delta) 
    {
      SEP_List[m][n][0] = -1;
      SEP_List[m][n][1] = -1;
      SEP_List[m][n][2] = -1;
      
      return 0; //greater than largest possible separation, definitely not overlapping or checking the particle itself, i.e., m = n
    }
  else if(temp_dis<d_cfmin &&  temp_dis>=Delta) 
    {
      SEP_List[m][n][0] = -1;
      SEP_List[m][n][1] = -1;
      SEP_List[m][n][2] = -1;

      return 1; //smaller than the smallest possible separation, definitely overlapping
    }


  //if in the intermediate distance, need to check with the separation axis theorem
  //need to do the complete check if overlapping: the most inefficient cases...


  /*
  double temp_VdistEII[SD]; //the vector from m to n
  for(int k=0; k<SD; k++)
    temp_VdistEII[k] = -temp_VdisE[k];
  */
  //*******************************************************************

  //first, check the faces of Pn, i.e., loop over all faces of N and vertex of M...
  for(int i=0; i<Pn.n_face; i++) //the face normal
    {
      //the efficiency for face check can be further improved by 
      //first computing the inner product of face normal and centroid displacement
      //only when the inner product is positive that we need to check for projections...

      //compute the normal 
      Pn.GetNormal_F(i);
      
      if(Comput.GetInnerProduct(Pn.Normal[i], temp_VdisE) <= 0.0)
	continue;

      
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

	  SEP_List[m][n][0] = 0;
	  SEP_List[m][n][1] = 0;
	  SEP_List[m][n][2] = i;

	  return overlap_flag;
	}


    }//end of loop for the faces...

 

  //second, check the faces of Pm  
  for(int i=0; i<Pm.n_face; i++)
   {

     //similar to the previous case, first check the direction of the normal

     //compute the normal to be checked
     Pm.GetNormal_F(i);
   

     if(Comput.GetInnerProduct(Pm.Normal[i], temp_VdistEII) <= 0.0)
       continue;


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

	  SEP_List[m][n][0] = 0;
	  SEP_List[m][n][1] = 1;
	  SEP_List[m][n][2] = i;
	  
	  
	  return overlap_flag;
	}

   }//end of loop for the faces...


  

  //*************************************************************
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //next we need to check all pair of edges to find the separating axis...
  //this can also been improved by first check ''edge normal'', i.e., by 
  // by computing the inner product of the edge normal and the centroid vector 
  // the edge normal is better to computed as an intrinsic geometric property of the particle

  double edgeM[Pm.n_edge][SD], edgeN[Pn.n_edge][SD];
  double edge_ctM[Pm.n_edge][SD], edge_ctN[Pn.n_edge][SD];
  //the center of the edge, for fast check, if point to the opposite direction from
  // the vector connecting particle centers, directly skip the pair

  for(int j=0; j<SD; j++)
    {
      for(int i=0; i<Pm.n_edge; i++)
	{
	  edgeM[i][j] = Pm.Vertex[Pm.Edge[i][0]][j] - Pm.Vertex[Pm.Edge[i][1]][j];
	  edge_ctM[i][j] = (Pm.Vertex[Pm.Edge[i][0]][j] + Pm.Vertex[Pm.Edge[i][1]][j])/2.0;
	}
      
      for(int i=0; i<Pn.n_edge; i++)
	{
	  edgeN[i][j] = Pn.Vertex[Pn.Edge[i][0]][j] - Pn.Vertex[Pn.Edge[i][1]][j];
	  edge_ctN[i][j] = (Pn.Vertex[Pn.Edge[i][0]][j] + Pn.Vertex[Pn.Edge[i][1]][j])/2.0;
	}
    }


  for(int i=0; i<Pm.n_edge; i++)
    for(int j=0; j<Pn.n_edge; j++)
      {

	//first check the edge center vector and particle centroid vector
	if(Comput.GetInnerProduct(edge_ctN[i], temp_VdisE)<=0.0 &&  Comput.GetInnerProduct(edge_ctM[i], temp_VdistEII)<=0.0)
	  continue;

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
	    
	if(GetMax(temp_disM, Pm.n_vertex)<GetMin(temp_disN, Pn.n_vertex) || GetMax(temp_disN, Pn.n_vertex)<GetMin(temp_disM, Pm.n_vertex))
	  {
	    overlap_flag = 0;

	    SEP_List[m][n][0] = 1;
	    SEP_List[m][n][1] = i;
	    SEP_List[m][n][2] = j;
	    
	    return overlap_flag;
	  }
      }//end of the loop for the edge pairs...
  
  
  //if we get here, then the M and N overlap...
  
  return overlap_flag;

}




int Packing::GlobalOverlapCheck(int a)
{
  int overlap_flag = 0;
  
  for(int m=0; m<N; m++)
    for(int n=m; n<N; n++)
      {
	//cout<<"In GOC "<<m<<" "<<n<<endl;
	
	for(int i=-1; i<=1; i++)
	  for(int j=-1; j<=1; j++)
	    for(int k=-1; k<=1; k++)
	      {
		
		
		if(CheckOverlap(PackPoly[n], PackPoly[m], CenterL[n], CenterL[m], i, j, k) == 1) 
		  {
		    printf("Polyhedron %d <%d, %d, %d> and Polyhedron %d overlap!\n", n, i, j, k, m);
		    //PrintOverlapPair(n, m);
		    overlap_flag = 1;
		  }
	      }
      }
  

  return overlap_flag;
}

//this is used in Move(), don't print out overlap pairs...
int Packing::GlobalOverlapCheck(double a)
{
  int overlap_flag = 0;
  
  for(int m=0; m<N; m++)
    for(int n=m; n<N; n++)
      {
	//cout<<"In GOC "<<m<<" "<<n<<endl;
	
	for(int i=-1; i<=1; i++)
	  for(int j=-1; j<=1; j++)
	    for(int k=-1; k<=1; k++)
	      {
		
		
		if(CheckOverlap(PackPoly[n], PackPoly[m], CenterL[n], CenterL[m], i, j, k) == 1) 
		  {
		    //printf("Polyhedron %d <%d, %d, %d> and Polyhedron %d overlap!\n", n, i, j, k, m);
		    //PrintOverlapPair(n, m);
		    overlap_flag = 1;
		  }
	      }
      }
  

  return overlap_flag;
}

//assuming the cell list is already there...
int Packing::GlobalOverlapCheck(Cells &BoxCell)
{
  int overlap_flag = 0; //assuming no overlapping...
  
  for(int m=0; m<N; m++)
    {
      
      int temp_indI = (int)floor(CenterL[m][0]*BoxCell.LXcell);
      int temp_indJ = (int)floor(CenterL[m][1]*BoxCell.LYcell);
      int temp_indK = (int)floor(CenterL[m][2]*BoxCell.LZcell);
      
      
      //cout<<" Neighbor of ["<<m<<"] = ( ";
      
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
		  
		  //cout<<temp_P<<"<"<<indexI<<","<<indexJ<<","<<indexK<<">"<<"\t";

		  if(CheckOverlap(PackPoly[temp_P], PackPoly[m], CenterL[temp_P], CenterL[m], indexI, indexJ, indexK)==1) //the particles overlap
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

void Packing::GetVParticle()
{
  VParticle = 0.0;

  for(int i=0; i<N; i++)
    {      
      VParticle += PackPoly[i].GetVol();      
    }
}

double Packing::GetDensity()
{
  double VLambda;
  VLambda = Lambda[0][2]*(Lambda[1][0]*Lambda[2][1]-Lambda[1][1]*Lambda[2][0]) + Lambda[1][2]*(Lambda[2][0]*Lambda[0][1] - Lambda[0][0]*Lambda[2][1]) + Lambda[2][2]*(Lambda[0][0]*Lambda[1][1] - Lambda[1][0]*Lambda[0][1]);
 
  return VParticle/VLambda;

}

void Packing::PrintOverlapPair(int m, int n)
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
  
  for(int j=0; j<PackPoly[m].n_vertex; j++)
    fprintf(fp, "%.12f\t%.12f\t%.12f\n", PackPoly[m].Vertex[j][0], PackPoly[m].Vertex[j][1], PackPoly[m].Vertex[j][2]);
  
  fprintf(fp, "\n");
  
  fprintf(fp, "%.12f\t%.12f\t%.12f\n", CenterL[n][0], CenterL[n][1], CenterL[n][2]);
  
  for(int j=0; j<PackPoly[n].n_vertex; j++)
    fprintf(fp, "%.12f\t%.12f\t%.12f\n", PackPoly[n].Vertex[j][0], PackPoly[n].Vertex[j][1], PackPoly[n].Vertex[j][2]);
  
  fprintf(fp, "\n");
 
  fclose(fp);
}


void Packing::PrintPacking()
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
      
      for(int j=0; j<PackPoly[i].n_vertex; j++)
	fprintf(fp, "%.12f\t%.12f\t%.12f\n", PackPoly[i].Vertex[j][0], PackPoly[i].Vertex[j][1], PackPoly[i].Vertex[j][2]);
      
      fprintf(fp, "\n");
    }
  
  fclose(fp);
}


void Packing::PrintDensity()
{
  FILE* fp;
  if((fp=fopen("rho.txt","a"))==NULL)
    {
      fp = fopen("rho.txt","w");
      fprintf(fp, "%f\n", GetDensity());
      fclose(fp);
    }
  else
    {
      fp=fopen("rho.txt","a");
      fprintf(fp, "%f\n", GetDensity());
      fclose(fp);
    }
}

void Packing::PrintDensity(double rho_t)
{
  FILE* fp;
  if((fp=fopen("rho.txt","a"))==NULL)
    {
      fp = fopen("rho.txt","w");
      fprintf(fp, "%f\n", rho_t);
      fclose(fp);
    }
  else
    {
      fp=fopen("rho.txt","a");
      fprintf(fp, "%f\n", rho_t);
      fclose(fp);
    }
}
