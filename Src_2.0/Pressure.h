


using namespace std;

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class Pressure{

 public:
  Pressure();
  void GetPressure();

};


Pressure::Pressure()
{
  cout<<"Object for class Pressure.h is constructed!"<<endl;
}

void Pressure::GetPressure()
{
  cout<<"Pressure will be sampled!"<<endl;
}
