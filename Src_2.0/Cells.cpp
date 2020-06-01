

using namespace std;

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "const.h"
#include "Cells.h"


Cells::Cells()
{
  //these are the default values
  LXcell = 3; 
  LYcell = 3;
  LZcell = 3;

  Ncell = LXcell*LYcell*LZcell;
}


int Cells::NumMax(int a, int b)
{
  if(a>b) return a;
  else return b;
}

void Cells::DelCellList()
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


void Cells::GetCellList(int flag, double &d0, double** &Lambda, double** &CenterL, int &N, double &temp_rho)
{
  
  if(flag != 0) //not the first time to establish CellList, so need to clean the list first
    {
      DelCellList();
    }

  //now need to compute the number of cell along each direction
  double LatX, LatY, LatZ; //the length of lattice vectors along each direction
  
  LatX = sqrt(Lambda[0][0]*Lambda[0][0]+Lambda[1][0]*Lambda[1][0]+Lambda[2][0]*Lambda[2][0]);
  LatY = sqrt(Lambda[0][1]*Lambda[0][1]+Lambda[1][1]*Lambda[1][1]+Lambda[2][1]*Lambda[2][1]);
  LatZ = sqrt(Lambda[0][2]*Lambda[0][2]+Lambda[1][2]*Lambda[1][2]+Lambda[2][2]*Lambda[2][2]);

  double linear_size = pow((LatX*LatY*LatZ)/(double)N, 1.0/3.0);

  LXcell = (int)NumMax(3, (int)floor(LatX/(1.3*linear_size)));
  LYcell = (int)NumMax(3, (int)floor(LatY/(1.3*linear_size)));
  LZcell = (int)NumMax(3, (int)floor(LatZ/(1.3*linear_size)));
  //the size of cell has to be large enough such that next neighbor cell can never contain possibly overlapped particles

  Ncell = LXcell*LYcell*LZcell;
  
  cout<<"The total number of cells Ncell = "<<Ncell<<endl;

  CellList = new node*[Ncell]; //re-setup the cell list
  for(int i=0; i<Ncell; i++)
    CellList[i] = NULL;

  for(int n=0; n<N; n++)
    {
      int temp_indexI = (int)floor(CenterL[n][0]*LXcell);
      int temp_indexJ = (int)floor(CenterL[n][1]*LYcell);
      int temp_indexK = (int)floor(CenterL[n][2]*LZcell);

      if(temp_indexI>=LXcell || temp_indexJ>=LYcell || temp_indexK>=LZcell)
	{
	  printf("CenterL includes invalid/un-scaled coordinates! Recheck!\n");
	  exit(1);
	}

      int cell_index = LXcell*LYcell*temp_indexK + LXcell*temp_indexJ + temp_indexI;

      //cout<<"Tetrah"<<n<<" cell_index = "<<cell_index<<endl;

      node* pt = new node;
      pt->index_poly = n;
      pt->next = CellList[cell_index];
      CellList[cell_index] = pt;

    }

  //record the density of the most recent update
  rho_cell = temp_rho;
  
}


void Cells::UpdateCellList(int m, double Old_CenterL[SD], double CenterLm[SD])
{
  //delete the particle from the old list...
  int temp_indI = (int)floor(Old_CenterL[0]*LXcell);
  int temp_indJ = (int)floor(Old_CenterL[1]*LYcell);
  int temp_indK = (int)floor(Old_CenterL[2]*LZcell);

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
	  if((pt1->index_poly)==m)
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
  temp_indI = (int)floor(CenterLm[0]*LXcell);
  temp_indJ = (int)floor(CenterLm[1]*LYcell);
  temp_indK = (int)floor(CenterLm[2]*LZcell);

  if(temp_indI>=LXcell || temp_indJ>=LYcell || temp_indK>=LZcell)
	{
	  printf("CenterL includes invalid/un-scaled coordinates! Recheck!\n");
	  exit(1);
	}

  int new_cell_index = LXcell*LYcell*temp_indK + LXcell*temp_indJ + temp_indI;

  node* pt = (node *)malloc(sizeof(node));
  pt->index_poly = m;
  pt->next = CellList[new_cell_index];

  CellList[new_cell_index] = pt;

}

void Cells::PrintCellList()
{
   for(int i=0; i<Ncell; i++)
     {
       node* pt = CellList[i];

       int temp_indexK = (int)floor((double)i/(double)(LXcell*LYcell));
       int temp_indexJ = (int)floor((double)(i%(LXcell*LYcell))/(double)LXcell);
       int temp_indexI = (i%(LXcell*LYcell))%LXcell;

       printf("CellList(%d,%d,%d) = ", temp_indexI, temp_indexJ, temp_indexK);
       
       while(pt!=NULL)
	 {
	   printf(" %d ", pt->index_poly);
	   pt = pt->next;
	 }
       
       printf("\n");
       
     }
   
}

Cells::~Cells()
{
  DelCellList();

  delete [] CellList;

  //cout<<"Destructor for Cells is called!"<<endl;
}
