

using namespace std;

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class OrderMetric{

 public:
  OrderMetric();
  void GetOrderMetric();

};


OrderMetric::OrderMetric()
{
  cout<<"Object for class OrderMetric.h is constructed!"<<endl;
}

void OrderMetric::GetOrderMetric()
{
  cout<<"OrderMetric will be sampled!"<<endl;
  cout<<"I decided to convert the packing coordinates to sphere packings and use Aleks' exe to compute order metrics... "<<endl;
}
