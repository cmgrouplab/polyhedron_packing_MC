//CPP file for obtaining statistics of the packing, including:
//     (1) pair correlation functions
//     (2) allignment correlation functions 
//     (3) face-normal correlation functions

//NOTE: the maximum sample length is relaxed from NLcounter/2.0 to MaxNL

using namespace std;

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


//the order is important!!!!
#include "const.h"
#include "Geometry.h"
#include "Polyhedron.h"
#include "Packing.h"
#include "Stat.h"


Stat::Stat()
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout<<"Object for class Stat is constructed!"<<endl;
  
}


Stat::Stat(int choice_flag, Packing & PK)
{
  //Now specify the parameters for statistics...
  
  if(choice_flag == 0)
    {
      cout<<" The number of bins per d_0: N_bin = "<<endl;
      cin >> N_bin; cout <<endl;
      cout<<" The maximal sample distance in terms of per bin_width: MaxNL = "<<endl;
      cin >> MaxNL; cout <<endl;
      
      bin = PK.d0/(double)N_bin;	  
      
    }
  else
    {
      N_bin = 25;
      MaxNL = 20*N_bin; //the range of correlation function is 20*d_0
      bin = PK.d0/(double)N_bin;
      
      cout<<" The number of bins per d_0: N_bin = "<< N_bin <<endl;
      cout<<" The maximal sample distance in terms of per bin_width: MaxNL = "<< MaxNL <<endl;
    }
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //initialize the correlation functions...
  /*
  G = new double**[MaxNL];
  for(int i=0; i<MaxNL; i++)
    {
      G[i] = new double*[MaxNL];
      
      for(int j=0; j<MaxNL; j++)
	G[i][j] = new double[MaxNL];
    }
  */

  g = new double[MaxNL];
  g_Num = new int[MaxNL];
  
  g_Normal = new double[MaxNL];
  g_Allign = new double[MaxNL];


  LocalOrder = new double[PK.N]; //face-normal order for each particle averaged over neigbors
  LocalAllign = new double[PK.N]; //allignment order for each particle averaged over neighbors

  NormalOP = new double*[PK.N];
  AllignOP = new double*[PK.N];

  for(int i=0; i<PK.N; i++)
    {
      NormalOP[i] = new double[PK.N];
      AllignOP[i] = new double[PK.N];

      LocalOrder[i] = 0.0; //contributed by nearest neighbors...
      LocalAllign[i] = 0.0;
    }
  
  
  //hopefully everything is properly initialized... if segmentation fault, go back and check these codes..

  FILE* fp;
  //fp = fopen("G2.txt","w");
  //fclose(fp);

  fp = fopen("g2.txt","w");
  fclose(fp);

  fp  = fopen("g_Normal.txt","w");
  fclose(fp);

  fp = fopen("Local_Normal.txt","w");
  fclose(fp);

  fp  = fopen("g_Allign.txt","w");
  fclose(fp);

  fp = fopen("Local_Allign.txt","w");
  fclose(fp);
  
}





double Stat::MinDis(int m, int n, Packing &PK, int IndexI, int IndexJ, int IndexK)
{
  //find the minimal distance between the centers of two polyhedra in Ecludean space...
  //by checking all the images of Poly[m], while keeping Poly[n] in the central box
  //record the index of the box which Poly[m] is in...
  //the order (m,n) is important, especially when getting the NNL.... 
  double dx = (PK.CenterL[m][0]-PK.CenterL[n][0]);
  double dy = (PK.CenterL[m][1]-PK.CenterL[n][1]);
  double dz = (PK.CenterL[m][2]-PK.CenterL[n][2]);

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
	      tempGd[q] += PK.Lambda[q][h]*tempLd[h];
	  
	  
	  //including the translational displacement...
	  double tempdist = tempGd[0]*tempGd[0]+tempGd[1]*tempGd[1]+tempGd[2]*tempGd[2]; 
	  
	  //printf("TempDist[%d][%d][%d] = %f\n", i, j, k, tempdist);
	  
	  if(tempdist<dist) // store the smallest distance...
	    {
	      dist = tempdist;
	      
     	      IndexI = i; IndexJ = j; IndexK = k;
	      
	    }
	  
	}
  
  // printf("IMAmndx = %f\t IMAmndy = %f\t IMAmndz = %f\n", IMAmn_Ldx, IMAmn_Ldy, IMAmn_Ldz);
  
  //printf("Dist = %f\n", dist);
  
  return sqrt(dist); //this is center-to-center distance...
}








double Stat::Get_NDensity(Packing &PK) // the number density, for computing g2...
{
  double VLambda;
  VLambda = PK.Lambda[0][2]*(PK.Lambda[1][0]*PK.Lambda[2][1]-PK.Lambda[1][1]*PK.Lambda[2][0]) + PK.Lambda[1][2]*(PK.Lambda[2][0]*PK.Lambda[0][1] - PK.Lambda[0][0]*PK.Lambda[2][1]) + PK.Lambda[2][2]*(PK.Lambda[0][0]*PK.Lambda[1][1] - PK.Lambda[1][0]*PK.Lambda[0][1]);

  return (double)PK.N/VLambda;
}








void Stat::Get_PairCorr(Packing &PK)
{

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //compute the minimum length of the lattice vectors
  double temp_d;
  double temp_dmin = 10000000000000.0;

  for(int i=0; i<SD; i++)
    {
      temp_d = 0;

      for(int j=0; j<SD; j++)
	temp_d += PK.Lambda[i][j]*PK.Lambda[i][j];

      if(temp_d < temp_dmin) temp_dmin = temp_d;
    }

  temp_dmin = sqrt(temp_dmin);

  //obtain the sample length....
  //note it's no problem if this is a little larger, the tail can always be cut-off
  NLcounter = 2*(int)floor(temp_dmin/bin);

  if(NLcounter>MaxNL) NLcounter = MaxNL;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /*
  for(int i=0; i<MaxNL; i++)
    for(int j=0; j<MaxNL; j++)
      for(int k=0; k<MaxNL; k++)
	G[i][j][k] = 0;
  */

  for(int i=0; i<MaxNL; i++)
    g[i] = 0;

  printf("Computing G2 and g2 now....\n");

  int IndexI, IndexJ, IndexK; //to indicate the index of the minimal images

  //loop over all particles
  for(int i=0; i<PK.N; i++)
    {
      for(int j=0; j<PK.N; j++)
	{
	  
	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  //for the radiual g2....
	  double temp_dis = MinDis(j, i, PK, IndexI, IndexJ, IndexK)/bin;
	
	  //the upper limit should be NLcounter/2.0, here relax this to MaxNL
	  // && temp_dis<(MaxNL-1)
	  if(temp_dis>0)
	    {
	      int int_dis = (int)floor(temp_dis);
	      double delta_dis = temp_dis - int_dis;
	      if(delta_dis>0.5) int_dis++;
	      
	      g[int_dis]++;
	      g_Num[int_dis]++;
	    }

	  /*
	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  //for the vector G2....
	  double temp_dL[SD];
	  temp_dL[0] = PK.CenterL[j][0] - PK.CenterL[i][0] + IndexI;
	  temp_dL[1] = PK.CenterL[j][1] - PK.CenterL[i][1] + IndexJ;
	  temp_dL[2] = PK.CenterL[j][2] - PK.CenterL[i][2] + IndexK;

	  double temp_dE[SD];
	  for(int u=0; u<SD; u++)
	    temp_dE[u] = 0;
	  
	  for(int u=0; u<SD; u++)
	    for(int v=0; v<SD; v++)
	      temp_dE[u] += PK.Lambda[u][v]*temp_dL[v];

	  double temp_dx = temp_dE[0]/bin + NLcounter/2.0;
	  double temp_dy = temp_dE[1]/bin + NLcounter/2.0;
	  double temp_dz = temp_dE[2]/bin + NLcounter/2.0;


	  if(temp_dx >=0 && temp_dx < NLcounter && temp_dy >=0 && temp_dy < NLcounter && temp_dz >=0 && temp_dz < NLcounter)
	    {
	      int int_dx = floor(temp_dx);
	      double delta_dx = temp_dx - int_dx;

	      if(delta_dx>0.5) int_dx++;

	      int int_dy = floor(temp_dy);
	      double delta_dy = temp_dy - int_dy;

	      if(delta_dy>0.5) int_dy++;

	      int int_dz = floor(temp_dz);
	      double delta_dz = temp_dz - int_dz;

	      if(delta_dz>0.5) int_dz++;

	      G[int_dx][int_dy][int_dz]++;
	    }
	  */

	}//end of the inner loop
    }//end of the outer loop

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //now we normalize g2 and print the results....
  double rho_n = Get_NDensity(PK);
  //double pi = 3.141592654;

  for(int r=0; r<NLcounter/2.0; r++)
    {
      g[r] = g[r]/(PK.N*rho_n*4*pi*(r+0.5)*(r+0.5)*bin*bin*bin);
    }

  /*
  FILE* fp;
  fp = fopen("G2.txt","a");
  for(int i=0; i<NLcounter; i++)
    for(int j=0; j<NLcounter; j++)
      for(int k=0; k<NLcounter; k++)
	fprintf(fp, "%f  %f  %f\t%f\n", bin*(i-NLcounter/2.0), bin*(j-NLcounter/2.0), bin*(k-NLcounter/2.0), G[i][j][k]);
  fprintf(fp, "******************************************************************\n");
  fclose(fp);
  */

  FILE* fp = fopen("g2.txt","a");
  for(int r=0; r<NLcounter/2.0; r++)
    fprintf(fp, "%f\t%f\n", (r)*bin, g[r]);
  fprintf(fp, "#******************************************************************\n");
  fclose(fp);

}








//NormalOP is simply the minus of the smallest face normal inner products of any two paritcles
//Here, we want to characterize the degree of f2f contact and do not distinguish faces with different types
void Stat::Get_NormalOP(Packing &PK) 
{  
  //loop for all the polyhedron pairs...
  for(int i=0; i<PK.N; i++)
    for(int j=0; j<=i; j++)
      {
	double Normal_InnerProduct[PK.PackPoly[i].n_face][PK.PackPoly[j].n_face]; //the inner product of the normals ...

    	//obtain the normal inner product....
	for(int p=0; p<PK.PackPoly[i].n_face; p++)
	  for(int q=0; q<PK.PackPoly[j].n_face; q++)
	    {
	      double temp_normip  = 0;
	      double temp_l1 = 0;
	      double temp_l2 = 0;

	      for(int k=0; k<SD; k++)
		{
	  	   temp_normip += PK.PackPoly[i].Normal[p][k]*PK.PackPoly[j].Normal[q][k];
		   temp_l1 += PK.PackPoly[i].Normal[p][k]*PK.PackPoly[i].Normal[p][k];
		   temp_l2 += PK.PackPoly[j].Normal[q][k]*PK.PackPoly[j].Normal[q][k];
		}


	      Normal_InnerProduct[p][q] = temp_normip/sqrt(temp_l1*temp_l2);
              //printf("Normal_IP[%d][%d] = %f\n", p, q, Normal_InnerProduct[p][q]);

	    }
	
	//compute the smallest of those inner product...

	double min_norm = 10000000000.0;

	for(int p=0; p<PK.PackPoly[i].n_face; p++)
	  for(int q=0; q<PK.PackPoly[j].n_face; q++)
	    {
	      if(Normal_InnerProduct[p][q]<min_norm) min_norm = Normal_InnerProduct[p][q];
	    }

	NormalOP[i][j] = NormalOP[j][i] = -min_norm; 
        //(fabs(min_norm)-sqrt(2.0/3.0))/(1-sqrt(2.0/3.0));
        //sqrt(2/3) is smallest value of cos(theta), theta is the angle between arbitary direction and the face normals...
        //-min_norm; //this could make the quantity a little more sensitive...hopefully
        // = fabs(min_norm-1.0)/2.0;
        //printf("NormalOP = %f, %f\n", NormalOP[i][j], NormalOP[j][i]);
      }
}







//compute the face-normal correlation distribution
//get the local order face-normal parameter for each particle, contributed from nearest neighbors
//then get the average over all particles
void Stat::Get_LocalNormOrder(Packing &PK)
{
  Get_NormalOP(PK);
  
  int neig_ct[PK.N]; //the number of neighbors for each particle..
  for(int i=0; i<PK.N; i++)
    neig_ct[i] = 0;
 

  //compute the distribution, add values to different bin and average over the number in that bin
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(int i=0; i<MaxNL; i++)
    {
      g_Normal[i] = 0;
    }
  //NOTE: g_Num[i] is already available...

  printf("We ASSUME the array g_Num[] is already set up...\n");
  printf("Computing Face Normal Correlation function now....\n");

  int tempI, tempJ, tempK; //for the function MinDis
 
  //loop over all particles
  for(int i=0; i<PK.N; i++)
    {
      for(int j=0; j<PK.N; j++)
	{
	  //for the radiual function....
	  double temp_dis = MinDis(j, i, PK, tempI, tempJ, tempK)/bin;

          //the upper limit should be NLcounter/2.0, here relax this to MaxNL
	  if(temp_dis>0 && temp_dis<(MaxNL-1))
	    {
	      int int_dis = (int)floor(temp_dis);
	      double delta_dis = temp_dis - int_dis;
	      if(delta_dis>0.5) int_dis++;

	      g_Normal[int_dis] += NormalOP[i][j];

	      //now check for local contribution
	      if(temp_dis<PK.d0/bin)
		{
		  LocalOrder[i] += NormalOP[i][j];
		  neig_ct[i]++;
		}
	    }

	}//end of the inner loop
    }//end of the outer loop

  //now we normalize g_Normal (divided by the number of pairs in that bin) and print the results....
  for(int r=0; r<NLcounter/2.0; r++)
    {
      if(g_Num[r]>0)
	{
	  g_Normal[r] = g_Normal[r]/(double)g_Num[r];
	}
      else
	{
	  g_Normal[r] = 0;
	}
    }
  
  FILE* fp  = fopen("g_Normal.txt","a");
  for(int r=0; r<NLcounter/2.0; r++)
    fprintf(fp, "%f\t%f\n", (r)*bin, g_Normal[r]);
  fprintf(fp, "*******************************************************************\n");
  fclose(fp);

   //compute the local normal order: first average over contact neighbors, then all particles
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Ave_LocalOrder = 0.0; //the average local order...
  

  //loop over all the particles...
  for(int i=0; i<PK.N; i++)
    {
      LocalOrder[i] = LocalOrder[i]/(double)neig_ct[i]; //get the averge LocalOrder for Poly[i] 
      Ave_LocalOrder += LocalOrder[i];  
    }
  
  //printf("AVELocalOrder = %f, Face = %f, Edge = %f\n", Ave_LocalOrder, Ave_FaceNum, Ave_EdgeNum);
  
  Ave_LocalOrder = Ave_LocalOrder/(double)PK.N;

 
  fp = fopen("Local_Normal.txt","a");
  fprintf(fp, "*******************************************************************\n");
  printf("The Averaged Local Normal Order Parameter: Ave_LOP = %f .\n", Ave_LocalOrder);
  fprintf(fp, "The Averaged Local Normal Order Parameter: Ave_LOP = %f .\n", Ave_LocalOrder);
  fprintf(fp, "#Index\tLocal_Normal_Order\n");
  for(int i=0; i<PK.N; i++)
    fprintf(fp, "%d\t%f\n", i, LocalOrder[i]);
  fprintf(fp, "*******************************************************************\n");
  fclose(fp);
  

}








//AllignOP is already the order parameter measuring the allignment of particle
//the current way allignOp is computed has given normal allignement more weight
//when any two faces of the same kind perfectly allign with each other, the two particles also allign (rigorous)
//we simply consider the faces  
void Stat::Get_AllignOP(Packing &PK) 
{

  //loop for all the polyhedron pairs...
  for(int i=0; i<PK.N; i++)
    for(int j=0; j<=i; j++)
      {
	double Normal_InnerProduct[PK.PackPoly[i].n_face][PK.PackPoly[j].n_face]; //the inner product of the normals ...
	double Vector_InnerProduct[PK.PackPoly[i].max_face_vert][PK.PackPoly[j].max_face_vert]; 
	                  //the inner product of center to vertex vector, once normal is fixed...

	//obtain the normal inner product....
	for(int p=0; p<PK.PackPoly[i].n_face; p++)
	  for(int q=0; q<PK.PackPoly[j].n_face; q++)
	    {
	     
	      double temp_normip  = 0;
	      double temp_l1 = 0;
	      double temp_l2 = 0;

	      for(int k=0; k<SD; k++)
		{
		  temp_normip += PK.PackPoly[i].Normal[p][k]*PK.PackPoly[j].Normal[q][k];
		  temp_l1 += PK.PackPoly[i].Normal[p][k]*PK.PackPoly[i].Normal[p][k];
		  temp_l2 += PK.PackPoly[j].Normal[q][k]*PK.PackPoly[j].Normal[q][k];
		}
	      

	      Normal_InnerProduct[p][q] = temp_normip/sqrt(temp_l1*temp_l2);
              //printf("Normal_IP[%d][%d] = %f\n", p, q, Normal_InnerProduct[p][q]);

	    }
	
	//compute the max of those inner product...
	
	double max_norm = -10000000.0;
	double max_vect = -10000000.0;
	int index_p, index_q;
	
	for(int p=0; p<PK.PackPoly[i].n_face; p++)
	  for(int q=0; q<PK.PackPoly[j].n_face; q++)
	    {
	      //make sure the face is of the same type...
	      if(Normal_InnerProduct[p][q]>max_norm && PK.PackPoly[i].face_type[p]==PK.PackPoly[j].face_type[q]) 
		{
		  max_norm = Normal_InnerProduct[p][q];
		  index_p = p;
		  index_q = q;
		}
	    }

	//now get
	double Temp_Vertex_I[PK.PackPoly[i].max_face_vert][SD];
	double Temp_Vertex_J[PK.PackPoly[j].max_face_vert][SD];
	int ct_vert_I = 0;
	int ct_vert_J = 0;

	for(int m=0; m<PK.PackPoly[i].vert_num_fspecie[PK.PackPoly[i].face_type[index_p]]; m++)
	  {
	    for(int n=0; n<SD; n++)
	      Temp_Vertex_I[ct_vert_I][n] = PK.PackPoly[i].Vertex[PK.PackPoly[i].Face[index_p][m]][n];
	    
	    ct_vert_I++;
	  }

	for(int m=0; m<PK.PackPoly[j].vert_num_fspecie[PK.PackPoly[j].face_type[index_q]]; m++)
	  {
	    for(int n=0; n<SD; n++)
	      Temp_Vertex_J[ct_vert_J][n] = PK.PackPoly[j].Vertex[PK.PackPoly[j].Face[index_q][m]][n];
	    
	    ct_vert_J++;
	  }
      

	for(int m=0; m<ct_vert_I; m++)
	  for(int n=0; n<ct_vert_J; n++)
	    {
	      double temp_vect_inner = 0.0;
	      double temp_l1 = 0.0;
	      double temp_l2 = 0.0;

	      for(int k=0; k<SD; k++)
		{
		  temp_vect_inner += Temp_Vertex_I[m][k]*Temp_Vertex_J[n][k];
		  temp_l1 += Temp_Vertex_I[m][k]*Temp_Vertex_I[m][k];
		  temp_l2 += Temp_Vertex_J[n][k]*Temp_Vertex_J[n][k];
		}
	      
	      
	      temp_vect_inner = temp_vect_inner/sqrt(temp_l1*temp_l2);

	      if(temp_vect_inner > max_vect)
		{
		  max_vect = temp_vect_inner;
		}
	    }
	
	AllignOP[i][j] = AllignOP[j][i] = (max_norm+max_vect)/2.0;
	
	
      }
}









//compute the allignment correlation distribution
//get the local allignment order parameter for each particle, contributed from nearest neighbors
//then get the average over all particles
void Stat::Get_LocalAllign(Packing &PK)
{

  Get_AllignOP(PK); 

  int neig_ct[PK.N]; //the number of neighbors for each particle..
  for(int i=0; i<PK.N; i++)
    neig_ct[i] = 0;
 

  //compute the distribution, add values to different bin and average over the number in that bin
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  for(int i=0; i<MaxNL; i++)
    {
      g_Allign[i] = 0;
    }
  //NOTE: g_Num[i] is already available...

  printf("We ASSUME the array g_Num[] is already set up...\n");
  printf("Computing Allignment Correlation function now....\n");

  int tempI, tempJ, tempK;

  //loop over all particles
  for(int i=0; i<PK.N; i++)
    {
      for(int j=0; j<PK.N; j++)
	{
	  //for the radiual function....
	  double temp_dis = MinDis(j, i, PK, tempI, tempJ, tempK)/bin;

          //the upper limit should be NLcounter/2.0, here relax this to MaxNL
	  if(temp_dis>0 && temp_dis<(MaxNL-1))
	    {
	      int int_dis = (int)floor(temp_dis);
	      double delta_dis = temp_dis - int_dis;
	      if(delta_dis>0.5) int_dis++;

	      g_Allign[int_dis] += AllignOP[i][j];

	      //now check for local contribution
	      if(temp_dis<PK.d0/bin)
		{
		  LocalAllign[i] += AllignOP[i][j];
		  neig_ct[i]++;
		}
	    }

	}//end of the inner loop
    }//end of the outer loop
  
  //now we normalize g_Allign (divided by the number of pairs in that bin) and print the results....
  
  for(int r=0; r<NLcounter/2.0; r++)
    {
      if(g_Num[r]>0)
	{
	  g_Allign[r] = g_Allign[r]/(double)g_Num[r];
	}
      else
	{
	  g_Allign[r] = 0;
	}
    }
  
  FILE* fp  = fopen("g_Allign.txt","a");
  for(int r=0; r<NLcounter/2.0; r++)
    fprintf(fp, "%f\t%f\n", (r)*bin, g_Allign[r]);
  fprintf(fp, "*******************************************************************\n");
  fclose(fp);


  //compute the local allignment order: first average over contact neighbors, then all particles
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Ave_LocalAllign = 0; //the average local order...

  //loop over all the particles...
  for(int i=0; i<PK.N; i++)
    {
      LocalAllign[i] = LocalAllign[i]/(double)neig_ct[i]; //get the averge LocalOrder for Poly[i] 
      Ave_LocalAllign += LocalAllign[i];  
    }
  
  //printf("AVELocalAllign = %f, Face = %f, Edge = %f\n", Ave_LocalAllign, Ave_FaceNum, Ave_EdgeNum);
  
  Ave_LocalAllign = Ave_LocalAllign/(double)PK.N;

  fp = fopen("Local_Allign.txt","a");
  fprintf(fp, "*******************************************************************\n");
  printf("The Averaged Local Allignment Order Parameter: Ave_LAllign = %f .\n", Ave_LocalAllign);
  fprintf(fp, "The Averaged Local Allignment Order Parameter: Ave_LAllign = %f .\n", Ave_LocalAllign);
  fprintf(fp, "#Index\tLocal_Normal_Order\n");
  for(int i=0; i<PK.N; i++)
    fprintf(fp, "%d\t%f\n", i, LocalAllign[i]);
  fprintf(fp, "*******************************************************************\n");
  fclose(fp);


 
}
