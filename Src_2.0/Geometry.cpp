using namespace std;

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "const.h" //this only includes spatial dimension
#include "Geometry.h"

void Geometry::GetCrossProduct(double Vect1[SD], double Vect2[SD], double Product[SD])
{
  Product[0] = Vect1[1]*Vect2[2] - Vect1[2]*Vect2[1];
  Product[1] = Vect1[2]*Vect2[0] - Vect1[0]*Vect2[2];
  Product[2] = Vect1[0]*Vect2[1] - Vect1[1]*Vect2[0];
}

double Geometry::GetInnerProduct(double Vect1[SD], double Vect2[SD])
{
  double sum = 0;

  for(int i=0; i<SD; i++)
    sum += Vect1[i]*Vect2[i];

  return sum;
}

//get the common 
void Geometry::GetPLine(double Direct[SD], double Vect[SD], double PLine[SD])
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


double Geometry::GetLength(double Vect[SD])
{
  double temp_length = 0;
  
  for(int j=0; j<SD; j++)
    temp_length += Vect[j]*Vect[j];
  
  temp_length = sqrt(temp_length);

  return temp_length;
}

void Geometry::Normalize(double Vect[SD])
{
  double temp_length = 0;
  
  for(int j=0; j<SD; j++)
    temp_length += Vect[j]*Vect[j];
  
  temp_length = sqrt(temp_length);
  
  for(int j=0; j<SD; j++)
    Vect[j] = Vect[j]/temp_length;
}

double Geometry::Normalize(double Vect[SD], int a)
{
  double temp_length = 0;
  
  for(int j=0; j<SD; j++)
    temp_length += Vect[j]*Vect[j];
  
  temp_length = sqrt(temp_length);
  
  for(int j=0; j<SD; j++)
    Vect[j] = Vect[j]/temp_length;

  return temp_length;
}
