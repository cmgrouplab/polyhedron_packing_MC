//compute the number variance as a function of R, which has to be smaller than L/2.0
//   using the analytical Sal and Frank derived in Phys. Rev. E 68, 041113 (2003), Eq.(58)

//author: Y. Jiao, yjiao@princeton.edu
//started: Oct. 21, 2010

//modified: Dec. 14, 2010
//now this works for three dimensional point configurations only, 
//BUT now we allow a non-cubic fundamental cell....

using namespace std;

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

const double Pi = 3.141592654;

int dim = 3; //spatial dimension, though only programmed the 2D formula here, make the code flexible to adpat...
int N; //number of particles
double L; //the box length, 435
double** lat;
double** coord; //the coordinates
double** dist; //the distance of all particle pairs
int mesh; //the increamental of R is L/mesh, 100

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double Get_Min(double Array[], int Size)
{
  double temp_min = 1000000.0;
  
  for(int i=0; i<Size; i++)
    {
      if(Array[i] < temp_min)
	temp_min = Array[i];
    }

  return temp_min;
}


double Get_Dist()
{
  for(int i=0; i<N; i++)
    for(int j=0; j<=i; j++)
      {
	double* dx = new double[dim];
	
	dist[i][j] = 0;

	double image_dist[27];
	int image_dist_ct = 0;
	
	//this only works for d=3...
	for(int n1=-1; n1<=1; n1++)
	  for(int n2=-1; n2<=1; n2++)
	    for(int n3=-1; n3<=1; n3++)
	      {
		//assume i is in the central box, j is in other box

		image_dist[image_dist_ct] = 0;

		for(int k=0; k<dim; k++)
		  {
		    dx[k] = coord[j][k]+(n1*lat[0][k]+n2*lat[1][k]+n3*lat[2][k])-coord[i][k] ;		   
		    
		    image_dist[image_dist_ct] += dx[k]*dx[k];
		  }

		image_dist_ct++;
	      }
	
	dist[i][j] = Get_Min(image_dist, image_dist_ct);

	dist[i][j] = sqrt(dist[i][j]);
	dist[j][i] = dist[i][j];
	
      }
}


double Get_Vol(double R)
{
  if(dim == 2)
    {
      return Pi*R*R;
    }
  else if(dim == 3)
    {
      return 4.0*Pi*R*R*R/3.0;
    }
  else
    {
      cout<<"This dimension is not programmed..."<<endl;
      exit(1);
    }
}

double Get_Alpha(int index_I, int index_J, double R)
{
  if(dim == 2)
    {
      double r_ij = dist[index_I][index_J];

      if(r_ij >= 2*R)
	return 0;
      else
	{
	  double t_val = acos(r_ij/(2*R))-(r_ij/(2*R))*sqrt((1-r_ij*r_ij/(4*R*R)));
	  return 2.0*t_val/Pi;
	}
    }
  else if(dim == 3)
    {
      double r_ij = dist[index_I][index_J];

      if(r_ij >= 2*R)
	return 0;
      else
	{
	  double t_val = 1-3*r_ij/(4*R)+(r_ij/R)*(r_ij/R)*(r_ij/R)/16.0;
	  return t_val;
	}
    }
  else
    {
      cout<<"This dimension is not programmed..."<<endl;
      exit(1);
    }
}

double Get_BoxVol()
{
  double Vlat;
  Vlat = lat[0][2]*(lat[1][0]*lat[2][1]-lat[1][1]*lat[2][0]) + lat[1][2]*(lat[2][0]*lat[0][1] - lat[0][0]*lat[2][1]) + lat[2][2]*(lat[0][0]*lat[1][1] - lat[1][0]*lat[0][1]);

  return Vlat;
}

double Get_Variance(double R)
{
  
  double rho = (double)N/Get_BoxVol(); //the number density
  double vol = Get_Vol(R); //the volume of the observation window

  double sum_alpha = 0.0;

  for(int i=0; i<N; i++)
    for(int j=0; j<i; j++)
      sum_alpha += 2*Get_Alpha(i, j, R);

  return rho*vol*(1-rho*vol+sum_alpha/(double)N);
  
}

main()
{
  //ask input for dim, L...
  //cout<<"Please specify the dimension, dim = "; cin>>dim;
  cout<<"Please specify the linear size, L = "; cin>>L;
  cout<<"Please specify the increamental of R = L/mesh, mesh = "; cin>>mesh;
  
  

  ifstream fin;
  fin.open("coords");
  if(!fin)
    {
      cout<<"Can not open target file coords! Abort!"<<endl;
      exit(1);
    }

  cout<<"read the number of points N..."<<endl;
  fin>>N;


  cout<<"read in the lattice vectors"<<endl;
  lat = new double*[dim];
  for(int i=0; i<dim; i++)
    lat[i] = new double[dim];

  //read in the lattice vectors...
  for(int i=0; i<dim; i++)
    {
      for(int j=0; j<dim; j++)
	fin>>lat[i][j];
    }
      
  cout<<"read in the point positions..."<<endl;
  //read in the position of points...
  coord = new double*[N];
  dist = new double*[N];

  for(int i=0; i<N; i++)
    {
      coord[i] = new double[dim];
      dist[i] = new double[N];

      for(int j=0; j<dim; j++)
	fin>>coord[i][j];
    }
  fin.close();

  Get_Dist();

  ofstream fout;
  fout.open("Num_Var.txt");

  double dR = L/(double)mesh;
  int Rmax = (int)floor(mesh/2.0);

  cout<<"Computing the number variance now..."<<endl;
  for(int r=1; r<Rmax; r++)
    {
      if(r*dR > L/4.0) break;

      double t_var = Get_Variance(r*dR);
      fout<<r*dR<<"   "<<t_var<<endl;
      cout<<r*dR<<"   "<<t_var<<endl;
      
    }

  fout.close();
}
