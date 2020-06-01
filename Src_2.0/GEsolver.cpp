//function definitions

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include "const.h"
#include "GEsolver.h"

GESOL::GESOL()
{
  moreeq = 1;
}

GESOL::~GESOL()
{
  int i,j; 
  for(i=0; i<neq; i++)
    free(*(sys+i));
  free(sys);
  free(back);

}

double GESOL::evaluate(int g) // get the final solutions
{
 int u;
 double rt=0;
 rt=sys[g][nun];
 for(u=0;u<nun;u++)
  {
   if(g==u) continue;
   rt-=sys[g][u]*back[u];
  }
 rt/=sys[g][g];
 return rt;
}

int GESOL::consistent() // whether the system is well defined
{
 int i,j,n=0;
 for(i=0;i<neq;i++)
  {
   for(j=0;j<nun;j++)
    {
     if(sys[i][j]!=0)
     {
      n++;
      break;
     }
    }
  }
  if(n==neq)
   return 1;
   return 0;
}

void GESOL::result()
{
int i,j;
 if(!consistent())
 {
  printf("\n\n Given system is inconsistent system...\n Solution can't be found\n");
  return;
 }
 else
 {
   for(i=nun-1;i>=0;i--)
    {
     back[i]=evaluate(i);
    }
   return;
 }
}

void GESOL::solution()
{
  int i,j,k,r;
  float c;
  for(r=0;r<neq;r++)
   {
    c=sys[r][r];

      for(k=0;k<=nun;k++)
       {
	sys[r][k]/=c;
       }

     for(i=r+1;i<neq;i++)
      {
       c=sys[i][r];
	for(j=0;j<=nun;j++)
	 {
	  sys[i][j]-=c*sys[r][j];
	 }
      }
     
     /*
    for(i=0;i<neq;i++)
      {
       printf("\n\n");
	 for(j=0;j<=nun;j++)
	  {
	   printf("%15.4f",sys[i][j]);
	  }
      }
       printf("\n\n");
     */

       }
  result();
}


void GESOL::GEsolver(int Neq, int Nun, double Sys[SD][SD], double RB[SD])
{
 int i,j;
 neq = Neq;
 nun = Nun;

	 if(nun>neq)
	  {
	   printf("Data inadequate to find the solution.....\n");
	   exit(0);
	  }
	 if(neq>nun)
	  moreeq=1;
	  sys=(double **)malloc((neq)*sizeof(double *));
	  back=(double *)malloc((nun)*sizeof(double ));
	 for(i=0;i<neq;i++)
	 {
	   *(sys+i)=(double *)malloc((nun+1)*sizeof(double));
	 }


	 //printf("INITIALIZING THE COEFFICIENTS OF EQUATIONS...\n");
	 for(i=0;i<neq;i++)
	  {
	   for(j=0;j<nun;j++)
	   {
	    *(*(sys+i)+j) = Sys[i][j];
	    back[i]=0;
	   }
	   //right hand side: b
	   *(*(sys+i)+j) = RB[i];
	  }
	 

	 //deal with the case when we have zero diagonal elements...
	 //this simple version only works for 3*3 system...
	 double TempSys[SD+1];
	 //double TempRB;
	 int order_index[SD];
	 int temp_index;
	 order_index[0] = 0; order_index[1] = 1; order_index[2] = 2; 
	 double zero = 0.0000001;

	 if(fabs(sys[0][0])<=zero)
	   {
	     if(fabs(sys[1][0])>zero && fabs(sys[0][1])>zero) 
	       {
		 for(i=0; i<SD+1; i++)
		   TempSys[i] = sys[0][i];
		 for(i=0; i<SD+1; i++)
		   sys[0][i] = sys[1][i];
		 for(i=0; i<SD+1; i++)
		   sys[1][i] = TempSys[i];

		 temp_index = order_index[0];
		 order_index[0] = order_index[1];
		 order_index[1] = temp_index;
	       }
	     else
	       {
                 for(i=0; i<SD+1; i++)
		   TempSys[i] = sys[0][i];
		 for(i=0; i<SD+1; i++)
		   sys[0][i] = sys[2][i];
		 for(i=0; i<SD+1; i++)
		   sys[2][i] = TempSys[i];

		 temp_index = order_index[0];
		 order_index[0] = order_index[2];
		 order_index[2] = temp_index;
	       }
	   }


	 if(fabs(sys[1][1])<=zero)
	   {
	     if(fabs(sys[0][1])>zero && fabs(sys[1][0])>zero) 
	       {
		 for(i=0; i<SD+1; i++)
		   TempSys[i] = sys[1][i];
		 for(i=0; i<SD+1; i++)
		   sys[1][i] = sys[0][i];
		 for(i=0; i<SD+1; i++)
		   sys[0][i] = TempSys[i];

		 temp_index = order_index[1];
		 order_index[1] = order_index[0];
		 order_index[0] = temp_index;
	       }
	     else
	       {
                 for(i=0; i<SD+1; i++)
		   TempSys[i] = sys[1][i];
		 for(i=0; i<SD+1; i++)
		   sys[1][i] = sys[2][i];
		 for(i=0; i<SD+1; i++)
		   sys[2][i] = TempSys[i];

		 temp_index = order_index[1];
		 order_index[1] = order_index[2];
		 order_index[2] = temp_index;
	       }
	   }

          if(fabs(sys[2][2])<=zero)
	   {
	     if(fabs(sys[0][2])>zero && fabs(sys[2][0])>zero) 
	       {
		 for(i=0; i<SD+1; i++)
		   TempSys[i] = sys[2][i];
		 for(i=0; i<SD+1; i++)
		   sys[2][i] = sys[0][i];
		 for(i=0; i<SD+1; i++)
		   sys[0][i] = TempSys[i];

		 temp_index = order_index[2];
		 order_index[2] = order_index[0];
		 order_index[0] = temp_index;
	       }
	     else
	       {
                 for(i=0; i<SD+1; i++)
		   TempSys[i] = sys[2][i];
		 for(i=0; i<SD+1; i++)
		   sys[2][i] = sys[1][i];
		 for(i=0; i<SD+1; i++)
		   sys[1][i] = TempSys[i];

		 temp_index = order_index[2];
		 order_index[2] = order_index[1];
		 order_index[1] = temp_index;
	       }
	   }




	  /*
	  printf("GIVEN SYSTEM OF EQUATION\n\n");
	 for(i=0;i<neq;i++)
	  {
	   printf("\t");
	   for(j=0;j<nun;j++)
	    {
	     if(*(*(sys+i)+j)<0)
	     {
	     printf("\b");
	     printf("\b");
	     }
	     printf("%f",*(*(sys+i)+j));
	     printf(" x");
	     printf("%d",j+1);
	     if(!(j==nun-1))
	     printf(" + ");
	    }
	   printf(" = ");
	   printf("%f",*(*(sys+i)+j));
	   printf("\n");
	   }
	  */

     solution();
     
     
     //printf("\nSOLUTION FOR THE GIVEN SYSTEM OF LINEAR SYSTEM OF EQUATIONS\n");
     
     /*
     for(i=0;i<nun;i++)
      {
	//printf("\t");
       printf("x");
       printf("%d", i+1);
       printf(" = ");
       printf("%f \n", back[i]);

      }
     */
     
         
}
