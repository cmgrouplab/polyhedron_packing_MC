//Read in lattice vector from standard input 
//for each edge generate a hollow square tube to show the simulation box


using namespace std;

#include <iostream>
#include <fstream>
#include <math.h>



double vertex[8][3]; //the vertex of the box
double node[96][3]; //the total number of nodes for the 12 tubes
double face[48][4]; //the four nodes for each tube

double temp_node[8][3];

ifstream fin;
ofstream fout;

//read in from the lattice vectors
void Get_Vertex()
{
  cout<<"read lattice vectors from standard input"<<endl;

  for(int i=0; i<3; i++)
    cin>>vertex[1][i];
  for(int i=0; i<3; i++)
    cin>>vertex[3][i];
  for(int i=0; i<3; i++)
    cin>>vertex[4][i];
  

  for(int i=0; i<3; i++)
    vertex[0][i] = 0;
  
  for(int i=0; i<3; i++)
    vertex[2][i] = vertex[1][i]+vertex[3][i];
  for(int i=0; i<3; i++)
    vertex[5][i] = vertex[1][i]+vertex[4][i];
  for(int i=0; i<3; i++)
    vertex[7][i] = vertex[3][i]+vertex[4][i];
  for(int i=0; i<3; i++)
    vertex[6][i] = vertex[2][i]+vertex[4][i];
}

double Get_Length(double vect[3])
{
  double sum = 0.0;

  for(int i=0; i<3; i++)
    sum += vect[i]*vect[i];

  return sqrt(sum);
}


void Make_Tubes(double len0)
{
  int node_ct = 0;
  int face_ct = 0;

  double temp_len1, temp_len2;
  double temp_vect1[3];
  double temp_vect2[3];

  //****************************************************************
  //make the tube with 0-1
  for(int i=0; i<3; i++)
    {
      temp_vect1[i] = vertex[2][i]-vertex[1][i];
      temp_vect2[i] = vertex[5][i]-vertex[1][i];
    }
  temp_len1 = Get_Length(temp_vect1);
  temp_len2 = Get_Length(temp_vect2);

  for(int i=0; i<3; i++)
    {
      temp_node[0][i] = vertex[0][i];
      temp_node[1][i] = vertex[1][i];

      temp_node[2][i] = temp_node[1][i] + temp_vect1[i]*len0/temp_len1;
      temp_node[3][i] = temp_node[0][i] + temp_vect1[i]*len0/temp_len1;

      temp_node[5][i] = temp_node[1][i] + temp_vect2[i]*len0/temp_len2;
      temp_node[4][i] = temp_node[0][i] + temp_vect2[i]*len0/temp_len2;

      temp_node[6][i] = temp_node[5][i] + temp_vect1[i]*len0/temp_len1;
      temp_node[7][i] = temp_node[4][i] + temp_vect1[i]*len0/temp_len1;
    }

  //write the nodes according to their global index
  for(int i=0; i<8; i++)
    {
      for(int j=0; j<3; j++)
	{
	  node[i+node_ct][j] = temp_node[i][j];
	}
    }
  
  //write the faces according to their gloabl index
  face[face_ct+0][0] = node_ct+0; face[face_ct+0][1] = node_ct+1; face[face_ct+0][2] = node_ct+2; face[face_ct+0][3] = node_ct+3;   
  face[face_ct+1][0] = node_ct+0; face[face_ct+1][1] = node_ct+1; face[face_ct+1][2] = node_ct+5; face[face_ct+1][3] = node_ct+4;   
  face[face_ct+2][0] = node_ct+4; face[face_ct+2][1] = node_ct+5; face[face_ct+2][2] = node_ct+6; face[face_ct+2][3] = node_ct+7;   
  face[face_ct+3][0] = node_ct+3; face[face_ct+3][1] = node_ct+2; face[face_ct+3][2] = node_ct+6; face[face_ct+3][3] = node_ct+7;   

  node_ct += 8;
  face_ct += 4;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //now translate this edge to vertex 3, 4 and 7
  for(int i=0; i<8; i++)
    {
      for(int j=0; j<3; j++)
	node[i+node_ct][j] = temp_node[i][j] + vertex[3][j];
    }

  face[face_ct+0][0] = node_ct+0; face[face_ct+0][1] = node_ct+1; face[face_ct+0][2] = node_ct+2; face[face_ct+0][3] = node_ct+3;   
  face[face_ct+1][0] = node_ct+0; face[face_ct+1][1] = node_ct+1; face[face_ct+1][2] = node_ct+5; face[face_ct+1][3] = node_ct+4;   
  face[face_ct+2][0] = node_ct+4; face[face_ct+2][1] = node_ct+5; face[face_ct+2][2] = node_ct+6; face[face_ct+2][3] = node_ct+7;   
  face[face_ct+3][0] = node_ct+3; face[face_ct+3][1] = node_ct+2; face[face_ct+3][2] = node_ct+6; face[face_ct+3][3] = node_ct+7;   

  node_ct += 8;
  face_ct += 4;

  for(int i=0; i<8; i++)
    {
      for(int j=0; j<3; j++)
	node[i+node_ct][j] = temp_node[i][j] + vertex[4][j];
    }

  face[face_ct+0][0] = node_ct+0; face[face_ct+0][1] = node_ct+1; face[face_ct+0][2] = node_ct+2; face[face_ct+0][3] = node_ct+3;   
  face[face_ct+1][0] = node_ct+0; face[face_ct+1][1] = node_ct+1; face[face_ct+1][2] = node_ct+5; face[face_ct+1][3] = node_ct+4;   
  face[face_ct+2][0] = node_ct+4; face[face_ct+2][1] = node_ct+5; face[face_ct+2][2] = node_ct+6; face[face_ct+2][3] = node_ct+7;   
  face[face_ct+3][0] = node_ct+3; face[face_ct+3][1] = node_ct+2; face[face_ct+3][2] = node_ct+6; face[face_ct+3][3] = node_ct+7;   

  node_ct += 8;
  face_ct += 4;

  for(int i=0; i<8; i++)
    {
      for(int j=0; j<3; j++)
	node[i+node_ct][j] = temp_node[i][j] + vertex[7][j];
    }

  face[face_ct+0][0] = node_ct+0; face[face_ct+0][1] = node_ct+1; face[face_ct+0][2] = node_ct+2; face[face_ct+0][3] = node_ct+3;   
  face[face_ct+1][0] = node_ct+0; face[face_ct+1][1] = node_ct+1; face[face_ct+1][2] = node_ct+5; face[face_ct+1][3] = node_ct+4;   
  face[face_ct+2][0] = node_ct+4; face[face_ct+2][1] = node_ct+5; face[face_ct+2][2] = node_ct+6; face[face_ct+2][3] = node_ct+7;   
  face[face_ct+3][0] = node_ct+3; face[face_ct+3][1] = node_ct+2; face[face_ct+3][2] = node_ct+6; face[face_ct+3][3] = node_ct+7;   
  
  node_ct += 8;
  face_ct += 4;


  //****************************************************************
  //make the tube with 0-3
  for(int i=0; i<3; i++)
    {
      temp_vect1[i] = vertex[1][i]-vertex[0][i];
      temp_vect2[i] = vertex[4][i]-vertex[0][i];
    }
  temp_len1 = Get_Length(temp_vect1);
  temp_len2 = Get_Length(temp_vect2);

  for(int i=0; i<3; i++)
    {
      temp_node[0][i] = vertex[0][i];
      temp_node[3][i] = vertex[3][i];

      temp_node[1][i] = temp_node[0][i] + temp_vect1[i]*len0/temp_len1;
      temp_node[2][i] = temp_node[3][i] + temp_vect1[i]*len0/temp_len1;

      temp_node[4][i] = temp_node[0][i] + temp_vect2[i]*len0/temp_len2;
      temp_node[7][i] = temp_node[3][i] + temp_vect2[i]*len0/temp_len2;

      temp_node[5][i] = temp_node[4][i] + temp_vect1[i]*len0/temp_len1;
      temp_node[6][i] = temp_node[7][i] + temp_vect1[i]*len0/temp_len1;
    }

  //write the nodes according to their global index
  for(int i=0; i<8; i++)
    {
      for(int j=0; j<3; j++)
	{
	  node[i+node_ct][j] = temp_node[i][j];
	}
    }
  
  //write the faces according to their gloabl index
  face[face_ct+0][0] = node_ct+0; face[face_ct+0][1] = node_ct+1; face[face_ct+0][2] = node_ct+2; face[face_ct+0][3] = node_ct+3;   
  face[face_ct+1][0] = node_ct+1; face[face_ct+1][1] = node_ct+2; face[face_ct+1][2] = node_ct+6; face[face_ct+1][3] = node_ct+5;   
  face[face_ct+2][0] = node_ct+0; face[face_ct+2][1] = node_ct+3; face[face_ct+2][2] = node_ct+7; face[face_ct+2][3] = node_ct+4;   
  face[face_ct+3][0] = node_ct+4; face[face_ct+3][1] = node_ct+5; face[face_ct+3][2] = node_ct+6; face[face_ct+3][3] = node_ct+7;   

  node_ct += 8;
  face_ct += 4;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //now translate this edge to vertex 1, 4 and 5
  for(int i=0; i<8; i++)
    {
      for(int j=0; j<3; j++)
	node[i+node_ct][j] = temp_node[i][j] + vertex[1][j];
    }
  face[face_ct+0][0] = node_ct+0; face[face_ct+0][1] = node_ct+1; face[face_ct+0][2] = node_ct+2; face[face_ct+0][3] = node_ct+3;   
  face[face_ct+1][0] = node_ct+1; face[face_ct+1][1] = node_ct+2; face[face_ct+1][2] = node_ct+6; face[face_ct+1][3] = node_ct+5;   
  face[face_ct+2][0] = node_ct+0; face[face_ct+2][1] = node_ct+3; face[face_ct+2][2] = node_ct+7; face[face_ct+2][3] = node_ct+4;   
  face[face_ct+3][0] = node_ct+4; face[face_ct+3][1] = node_ct+5; face[face_ct+3][2] = node_ct+6; face[face_ct+3][3] = node_ct+7;   
 

  node_ct += 8;
  face_ct += 4;

  for(int i=0; i<8; i++)
    {
      for(int j=0; j<3; j++)
	node[i+node_ct][j] = temp_node[i][j] + vertex[4][j];
    }
  face[face_ct+0][0] = node_ct+0; face[face_ct+0][1] = node_ct+1; face[face_ct+0][2] = node_ct+2; face[face_ct+0][3] = node_ct+3;   
  face[face_ct+1][0] = node_ct+1; face[face_ct+1][1] = node_ct+2; face[face_ct+1][2] = node_ct+6; face[face_ct+1][3] = node_ct+5;   
  face[face_ct+2][0] = node_ct+0; face[face_ct+2][1] = node_ct+3; face[face_ct+2][2] = node_ct+7; face[face_ct+2][3] = node_ct+4;   
  face[face_ct+3][0] = node_ct+4; face[face_ct+3][1] = node_ct+5; face[face_ct+3][2] = node_ct+6; face[face_ct+3][3] = node_ct+7;   


  node_ct += 8;
  face_ct += 4;

  for(int i=0; i<8; i++)
    {
      for(int j=0; j<3; j++)
	node[i+node_ct][j] = temp_node[i][j] + vertex[5][j];
    }
  face[face_ct+0][0] = node_ct+0; face[face_ct+0][1] = node_ct+1; face[face_ct+0][2] = node_ct+2; face[face_ct+0][3] = node_ct+3;   
  face[face_ct+1][0] = node_ct+1; face[face_ct+1][1] = node_ct+2; face[face_ct+1][2] = node_ct+6; face[face_ct+1][3] = node_ct+5;   
  face[face_ct+2][0] = node_ct+0; face[face_ct+2][1] = node_ct+3; face[face_ct+2][2] = node_ct+7; face[face_ct+2][3] = node_ct+4;   
  face[face_ct+3][0] = node_ct+4; face[face_ct+3][1] = node_ct+5; face[face_ct+3][2] = node_ct+6; face[face_ct+3][3] = node_ct+7;   

  
  node_ct += 8;
  face_ct += 4;


  //****************************************************************
  //make the tube with 0-4
  for(int i=0; i<3; i++)
    {
      temp_vect1[i] = vertex[1][i]-vertex[0][i];
      temp_vect2[i] = vertex[3][i]-vertex[0][i];
    }
  temp_len1 = Get_Length(temp_vect1);
  temp_len2 = Get_Length(temp_vect2);

  for(int i=0; i<3; i++)
    {
      temp_node[0][i] = vertex[0][i];
      temp_node[4][i] = vertex[4][i];

      temp_node[1][i] = temp_node[0][i] + temp_vect1[i]*len0/temp_len1;
      temp_node[5][i] = temp_node[4][i] + temp_vect1[i]*len0/temp_len1;

      temp_node[3][i] = temp_node[0][i] + temp_vect2[i]*len0/temp_len2;
      temp_node[7][i] = temp_node[4][i] + temp_vect2[i]*len0/temp_len2;

      temp_node[2][i] = temp_node[3][i] + temp_vect1[i]*len0/temp_len1;
      temp_node[6][i] = temp_node[7][i] + temp_vect1[i]*len0/temp_len1;
    }

  //write the nodes according to their global index
  for(int i=0; i<8; i++)
    {
      for(int j=0; j<3; j++)
	{
	  node[i+node_ct][j] = temp_node[i][j];
	}
    }
  
  //write the faces according to their gloabl index
  face[face_ct+0][0] = node_ct+0; face[face_ct+0][1] = node_ct+1; face[face_ct+0][2] = node_ct+5; face[face_ct+0][3] = node_ct+4;   
  face[face_ct+1][0] = node_ct+1; face[face_ct+1][1] = node_ct+2; face[face_ct+1][2] = node_ct+6; face[face_ct+1][3] = node_ct+5;   
  face[face_ct+2][0] = node_ct+3; face[face_ct+2][1] = node_ct+2; face[face_ct+2][2] = node_ct+6; face[face_ct+2][3] = node_ct+7;   
  face[face_ct+3][0] = node_ct+0; face[face_ct+3][1] = node_ct+3; face[face_ct+3][2] = node_ct+7; face[face_ct+3][3] = node_ct+4;   

  node_ct += 8;
  face_ct += 4;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //now translate this edge to vertex 1, 2 and 3
  for(int i=0; i<8; i++)
    {
      for(int j=0; j<3; j++)
	node[i+node_ct][j] = temp_node[i][j] + vertex[1][j];
    }
  face[face_ct+0][0] = node_ct+0; face[face_ct+0][1] = node_ct+1; face[face_ct+0][2] = node_ct+5; face[face_ct+0][3] = node_ct+4;   
  face[face_ct+1][0] = node_ct+1; face[face_ct+1][1] = node_ct+2; face[face_ct+1][2] = node_ct+6; face[face_ct+1][3] = node_ct+5;   
  face[face_ct+2][0] = node_ct+3; face[face_ct+2][1] = node_ct+2; face[face_ct+2][2] = node_ct+6; face[face_ct+2][3] = node_ct+7;   
  face[face_ct+3][0] = node_ct+0; face[face_ct+3][1] = node_ct+3; face[face_ct+3][2] = node_ct+7; face[face_ct+3][3] = node_ct+4;   

  node_ct += 8;
  face_ct += 4;

  for(int i=0; i<8; i++)
    {
      for(int j=0; j<3; j++)
	node[i+node_ct][j] = temp_node[i][j] + vertex[2][j];
    }
  face[face_ct+0][0] = node_ct+0; face[face_ct+0][1] = node_ct+1; face[face_ct+0][2] = node_ct+5; face[face_ct+0][3] = node_ct+4;   
  face[face_ct+1][0] = node_ct+1; face[face_ct+1][1] = node_ct+2; face[face_ct+1][2] = node_ct+6; face[face_ct+1][3] = node_ct+5;   
  face[face_ct+2][0] = node_ct+3; face[face_ct+2][1] = node_ct+2; face[face_ct+2][2] = node_ct+6; face[face_ct+2][3] = node_ct+7;   
  face[face_ct+3][0] = node_ct+0; face[face_ct+3][1] = node_ct+3; face[face_ct+3][2] = node_ct+7; face[face_ct+3][3] = node_ct+4;   

  node_ct += 8;
  face_ct += 4;

  for(int i=0; i<8; i++)
    {
      for(int j=0; j<3; j++)
	node[i+node_ct][j] = temp_node[i][j] + vertex[3][j];
    }
  face[face_ct+0][0] = node_ct+0; face[face_ct+0][1] = node_ct+1; face[face_ct+0][2] = node_ct+5; face[face_ct+0][3] = node_ct+4;   
  face[face_ct+1][0] = node_ct+1; face[face_ct+1][1] = node_ct+2; face[face_ct+1][2] = node_ct+6; face[face_ct+1][3] = node_ct+5;   
  face[face_ct+2][0] = node_ct+3; face[face_ct+2][1] = node_ct+2; face[face_ct+2][2] = node_ct+6; face[face_ct+2][3] = node_ct+7;   
  face[face_ct+3][0] = node_ct+0; face[face_ct+3][1] = node_ct+3; face[face_ct+3][2] = node_ct+7; face[face_ct+3][3] = node_ct+4;   

  node_ct += 8;
  face_ct += 4;

  cout<<"node_ct = "<<node_ct<<endl;
  cout<<"face_ct = "<<face_ct<<endl;

}

void Print_OFF()
{
  fout.open("Box.off");

  fout<<"OFF"<<endl;
  fout<<"96 48 0"<<endl;
  
  //print out the vertices
  for(int i=0; i<96; i++)
    {
      for(int j=0; j<3; j++)
	fout<<node[i][j]<<"  ";
      fout<<endl;
    }

  //print out the faces
  for(int i=0; i<48; i++)
    {
      fout<<4<<"    ";
      for(int j=0; j<4; j++)
	fout<<face[i][j]<<"  ";
      fout<<"255 0 0"<<endl;
    }

  fout.close();
  
}



main()
{
  Get_Vertex();

  //now make each square tube
  //first provide the length of the short edge of the tube
  double len0 = Get_Length(vertex[1])/50.0;

  Make_Tubes(len0);

  Print_OFF();
   
}
