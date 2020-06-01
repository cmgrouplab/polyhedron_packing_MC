
using namespace std;

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include "const.h" //specify spatial dimension
#include "Geometry.h"
#include "Polyhedron.h"

Polyhedron::Polyhedron()
{
  //this is the default constructor...
}


//IMPORTANT NOTE: for the codes specific for each shape, the vertex order is arbitary since
//                the face configuration is computed by search for nearest neighbors of each vertex
//                Here, we specify the face configuration, therefore vertex order is CRUCIAL!!!!!
Polyhedron::Polyhedron(char* temp_name)
{
  if(strcmp(temp_name,"tetrahedron")==0 || strcmp(temp_name,"P1")==0)
    {

      poly_name = "tetrahedron";
      n_vertex = 4;
      n_edge = 6;
      n_face = 4;
      n_fs = 1; //not really need this one, even for complicated shapes...
      max_face_vert = 3;

      face_num_fspecie = new int[n_fs];
      face_num_fspecie[0] = 4;

      vert_num_fspecie = new int[n_fs];
      vert_num_fspecie[0] = 3;
      
      //now initialize the data type
      Vertex = new double*[n_vertex];
      for(int i=0; i<n_vertex; i++)
	Vertex[i] = new double[SD];
      
      Edge = new int*[n_edge];
      for(int i=0; i<n_edge; i++)
	Edge[i] = new int[2];
      //now specify the index of vertices for each edge
      Edge[0][0] = 0; Edge[0][1] = 1;
      Edge[1][0] = 0; Edge[1][1] = 2;
      Edge[2][0] = 0; Edge[2][1] = 3;
      Edge[3][0] = 1; Edge[3][1] = 2;
      Edge[4][0] = 1; Edge[4][1] = 3;
      Edge[5][0] = 2; Edge[5][1] = 3;
      
      
      Face = new int*[n_face];
      Normal = new double*[n_face];  
      face_type = new int[n_face];    
      for(int i=0; i<n_face; i++)
	{
	  Face[i] = new int[max_face_vert]; //make sure the list is feasible for all faces
	  Normal[i] = new double[SD];
	  face_type[i] = 0;
	}
      Face[0][0] = 0; Face[0][1] = 1; Face[0][2] = 2;
      Face[1][0] = 0; Face[1][1] = 1; Face[1][2] = 3;
      Face[2][0] = 0; Face[2][1] = 2; Face[2][2] = 3;
      Face[3][0] = 1; Face[3][1] = 2; Face[3][2] = 3;
      
      d_cf = new double[n_face];
      d_cv = new double[n_vertex];
      
      
    } //add other polyhedron types here
  else if(strcmp(temp_name,"icosahedron")==0 || strcmp(temp_name,"P2")==0)
    {

      poly_name = "icosahedron";
      n_vertex = 12;
      n_edge = 30;
      n_face = 20;
      n_fs = 1; //only needed when printing out the configuration...
      max_face_vert = 3;

      face_num_fspecie = new int[n_fs];
      face_num_fspecie[0] = 20;

      vert_num_fspecie = new int[n_fs];
      vert_num_fspecie[0] = 3;
      
      //now initialize the data type
      Vertex = new double*[n_vertex];
      for(int i=0; i<n_vertex; i++)
	Vertex[i] = new double[SD];
      
      Edge = new int*[n_edge];
      for(int i=0; i<n_edge; i++)
	Edge[i] = new int[2];
      //now specify the index of vertices for each edge, based on the face configuration
      Edge[0][0] = 0; Edge[0][1] = 1;
      Edge[1][0] = 0; Edge[1][1] = 4;
      Edge[2][0] = 0; Edge[2][1] = 5;
      Edge[3][0] = 0; Edge[3][1] = 8;
      Edge[4][0] = 0; Edge[4][1] = 9;

      Edge[5][0] = 1; Edge[5][1] = 4;
      Edge[6][0] = 1; Edge[6][1] = 5;
      Edge[7][0] = 1; Edge[7][1] = 10;
      Edge[8][0] = 1; Edge[8][1] = 11;
      Edge[9][0] = 2; Edge[9][1] = 3;

      Edge[10][0] = 2; Edge[10][1] = 6;
      Edge[11][0] = 2; Edge[11][1] = 7;
      Edge[12][0] = 2; Edge[12][1] = 8;
      Edge[13][0] = 2; Edge[13][1] = 9;
      Edge[14][0] = 3; Edge[14][1] = 6;
     
      Edge[15][0] = 3; Edge[15][1] = 7;
      Edge[16][0] = 3; Edge[16][1] = 10;
      Edge[17][0] = 3; Edge[17][1] = 11;
      Edge[18][0] = 4; Edge[18][1] = 6;
      Edge[19][0] = 4; Edge[19][1] = 8;

      Edge[20][0] = 4; Edge[20][1] = 10;
      Edge[21][0] = 5; Edge[21][1] = 7;
      Edge[22][0] = 5; Edge[22][1] = 9;
      Edge[23][0] = 5; Edge[23][1] = 11;
      Edge[24][0] = 6; Edge[24][1] = 8;

      Edge[25][0] = 6; Edge[25][1] = 10;
      Edge[26][0] = 7; Edge[26][1] = 9;
      Edge[27][0] = 7; Edge[27][1] = 11;
      Edge[28][0] = 8; Edge[28][1] = 9;
      Edge[29][0] = 10; Edge[29][1] = 11;


      
      
      Face = new int*[n_face];
      Normal = new double*[n_face];   
      face_type = new int[n_face];
      for(int i=0; i<n_face; i++)
	{
	  Face[i] = new int[max_face_vert]; //make sure the list is feasible for all faces
	  Normal[i] = new double[SD];
	  face_type[i] = 0;
	}
      //the following list copied from AntiView output
      Face[0][0] = 0; Face[0][1] = 1; Face[0][2] = 4;
      Face[1][0] = 0; Face[1][1] = 1; Face[1][2] = 5;
      Face[2][0] = 0; Face[2][1] = 4; Face[2][2] = 8;
      Face[3][0] = 0; Face[3][1] = 5; Face[3][2] = 9;
      Face[4][0] = 0; Face[4][1] = 8; Face[4][2] = 9;
     
      Face[5][0] = 1; Face[5][1] = 4; Face[5][2] = 10;
      Face[6][0] = 1; Face[6][1] = 5; Face[6][2] = 11;
      Face[7][0] = 1; Face[7][1] = 10; Face[7][2] = 11;
      Face[8][0] = 2; Face[8][1] = 3; Face[8][2] = 6;
      Face[9][0] = 2; Face[9][1] = 3; Face[9][2] = 7;

      Face[10][0] = 2; Face[10][1] = 6; Face[10][2] = 8;
      Face[11][0] = 2; Face[11][1] = 7; Face[11][2] = 9;
      Face[12][0] = 2; Face[12][1] = 8; Face[12][2] = 9;
      Face[13][0] = 3; Face[13][1] = 6; Face[13][2] = 10;
      Face[14][0] = 3; Face[14][1] = 7; Face[14][2] = 11;

      Face[15][0] = 3; Face[15][1] = 10; Face[15][2] = 11;
      Face[16][0] = 4; Face[16][1] = 6; Face[16][2] = 8;
      Face[17][0] = 4; Face[17][1] = 6; Face[17][2] = 10;
      Face[18][0] = 5; Face[18][1] = 7; Face[18][2] = 9;
      Face[19][0] = 5; Face[19][1] = 7; Face[19][2] = 11;


      //the following is the AntiView output
      //3  0 1 4 
      //3  0 1 5 
      //3  0 4 8 
      //3  0 5 9 
      //3  0 8 9 

      //3  1 4 10 
      //3  1 5 11 
      //3  1 10 11 
      //3  2 3 6 
      //3  2 3 7 

      //3  2 6 8 
      //3  2 7 9 
      //3  2 8 9 
      //3  3 6 10 
      //3  3 7 11
 
      //3  3 10 11 
      //3  4 6 8 
      //3  4 6 10 
      //3  5 7 9 
      //3  5 7 11 
      
      d_cf = new double[n_face];
      d_cv = new double[n_vertex];
      
      
    } //add other polyhedron types here
  else if(strcmp(temp_name,"dodecahedron")==0 || strcmp(temp_name,"P3")==0)
    {

      poly_name = "dodecahedron";
      n_vertex = 20;
      n_edge = 30;
      n_face = 12;
      n_fs = 1; //not really need this one, even for complicated shapes...
      max_face_vert = 5;

      face_num_fspecie = new int[n_fs];
      face_num_fspecie[0] = 12;

      vert_num_fspecie = new int[n_fs];
      vert_num_fspecie[0] = 5;
      
      //now initialize the data type
      Vertex = new double*[n_vertex];
      for(int i=0; i<n_vertex; i++)
	Vertex[i] = new double[SD];
      
      Edge = new int*[n_edge];
      for(int i=0; i<n_edge; i++)
	Edge[i] = new int[2];
      //now specify the index of vertices for each edge     
      Edge[0][0] = 0; Edge[0][1] = 8;
      Edge[1][0] = 0; Edge[1][1] = 12;
      Edge[2][0] = 0; Edge[2][1] = 16;
      Edge[3][0] = 1; Edge[3][1] = 9;
      Edge[4][0] = 1; Edge[4][1] = 12;
      Edge[5][0] = 1; Edge[5][1] = 17;

      Edge[6][0] = 2; Edge[6][1] = 10;
      Edge[7][0] = 2; Edge[7][1] = 13;
      Edge[8][0] = 2; Edge[8][1] = 16;
      Edge[9][0] = 3; Edge[9][1] = 11;
      Edge[10][0] = 3; Edge[10][1] = 13;
      Edge[11][0] = 3; Edge[11][1] = 17;

      Edge[12][0] = 4; Edge[12][1] = 8;
      Edge[13][0] = 4; Edge[13][1] = 14;
      Edge[14][0] = 4; Edge[14][1] = 18;
      Edge[15][0] = 5; Edge[15][1] = 9;
      Edge[16][0] = 5; Edge[16][1] = 14;
      Edge[17][0] = 5; Edge[17][1] = 19;

      Edge[18][0] = 6; Edge[18][1] = 10;
      Edge[19][0] = 6; Edge[19][1] = 15;
      Edge[20][0] = 6; Edge[20][1] = 18;
      Edge[21][0] = 7; Edge[21][1] = 11;
      Edge[22][0] = 7; Edge[22][1] = 15;
      Edge[23][0] = 7; Edge[23][1] = 19;

      Edge[24][0] = 8; Edge[24][1] = 10;
      Edge[25][0] = 9; Edge[25][1] = 11;
      Edge[26][0] = 12; Edge[26][1] = 14;
      Edge[27][0] = 13; Edge[27][1] = 15;
      Edge[28][0] = 16; Edge[28][1] = 17;
      Edge[29][0] = 18; Edge[29][1] = 19;
      
      
      Face = new int*[n_face];
      Normal = new double*[n_face];  
      face_type = new int[n_face];
      for(int i=0; i<n_face; i++)
	{
	  Face[i] = new int[max_face_vert]; //make sure the list is feasible for all faces
	  Normal[i] = new double[SD];
	  face_type[i] = 0;
	}
      Face[0][0] = 3; Face[0][1] = 11; Face[0][2] = 7; Face[0][3] = 15; Face[0][4] = 13;
      Face[1][0] = 6; Face[1][1] = 15; Face[1][2] = 7; Face[1][3] = 19; Face[1][4] = 18;
      Face[2][0] = 5; Face[2][1] = 9; Face[2][2] = 11; Face[2][3] = 7; Face[2][4] = 19;
      Face[3][0] = 1; Face[3][1] = 9; Face[3][2] = 5; Face[3][3] = 14; Face[3][4] = 12;
      Face[4][0] = 0; Face[4][1] = 12; Face[4][2] = 1; Face[4][3] = 17; Face[4][4] = 16;
      Face[5][0] = 1; Face[5][1] = 9; Face[5][2] = 11; Face[5][3] = 3; Face[5][4] = 17;

      Face[6][0] = 2; Face[6][1] = 10; Face[6][2] = 6; Face[6][3] = 15; Face[6][4] = 13;
      Face[7][0] = 2; Face[7][1] = 13; Face[7][2] = 3; Face[7][3] = 17; Face[7][4] = 16;
      Face[8][0] = 0; Face[8][1] = 8; Face[8][2] = 10; Face[8][3] = 2; Face[8][4] = 16;
      Face[9][0] = 0; Face[9][1] = 8; Face[9][2] = 4; Face[9][3] = 14; Face[9][4] = 12;
      Face[10][0] = 4; Face[10][1] = 14; Face[10][2] = 5; Face[10][3] = 19; Face[10][4] = 18;
      Face[11][0] = 4; Face[11][1] = 8; Face[11][2] = 10; Face[11][3] = 6; Face[11][4] = 18;

      //the following is the face configuration from AntiView
      /*
	5  3 11 7 15 13
	5  6 15 7 19 18
	5  5 9 11 7 19
	5  1 9 5 14 12
	5  0 12 1 17 16
	5  1 9 11 3 17

	5  2 10 6 15 13
	5  2 13 3 17 16
	5  0 8 10 2 16
	5  0 8 4 14 12
	5  4 14 5 19 18
	5  4 8 10 6 18
      */



      d_cf = new double[n_face];
      d_cv = new double[n_vertex];
      
      
    } //add other polyhedron types here
  else if(strcmp(temp_name,"octahedron")==0 || strcmp(temp_name,"P4")==0)
    {

      poly_name = "octahedron";
      n_vertex = 6;
      n_edge = 12;
      n_face = 8;
      n_fs = 1; //not really need this one, even for complicated shapes...
      max_face_vert = 3;

      face_num_fspecie = new int[n_fs];
      face_num_fspecie[0] = 8;

      vert_num_fspecie = new int[n_fs];
      vert_num_fspecie[0] = 3;
      
      //now initialize the data type
      Vertex = new double*[n_vertex];
      for(int i=0; i<n_vertex; i++)
	Vertex[i] = new double[SD];
      
      Edge = new int*[n_edge];
      for(int i=0; i<n_edge; i++)
	Edge[i] = new int[2];
      //now specify the index of vertices for each edge
      Edge[0][0] = 0; Edge[0][1] = 2;
      Edge[1][0] = 0; Edge[1][1] = 3;
      Edge[2][0] = 0; Edge[2][1] = 4;
      Edge[3][0] = 0; Edge[3][1] = 5;

      Edge[4][0] = 1; Edge[4][1] = 2;
      Edge[5][0] = 1; Edge[5][1] = 3;
      Edge[6][0] = 1; Edge[6][1] = 4;
      Edge[7][0] = 1; Edge[7][1] = 5;

      Edge[8][0] = 2; Edge[8][1] = 4;
      Edge[9][0] = 2; Edge[9][1] = 5;
      Edge[10][0] = 3; Edge[10][1] = 4;
      Edge[11][0] = 3; Edge[11][1] = 5;
      
      
      Face = new int*[n_face];
      Normal = new double*[n_face]; 
      face_type = new int[n_face];
      for(int i=0; i<n_face; i++)
	{
	  Face[i] = new int[max_face_vert]; //make sure the list is feasible for all faces
	  Normal[i] = new double[SD];
	  face_type[i] = 0;
	}
      Face[0][0] = 0; Face[0][1] = 2; Face[0][2] = 4;
      Face[1][0] = 0; Face[1][1] = 2; Face[1][2] = 5;
      Face[2][0] = 0; Face[2][1] = 3; Face[2][2] = 4;
      Face[3][0] = 0; Face[3][1] = 3; Face[3][2] = 5;

      Face[4][0] = 1; Face[4][1] = 2; Face[4][2] = 4;
      Face[5][0] = 1; Face[5][1] = 2; Face[5][2] = 5;
      Face[6][0] = 1; Face[6][1] = 3; Face[6][2] = 4;
      Face[7][0] = 1; Face[7][1] = 3; Face[7][2] = 5;


      //face configuration taken from AntiView
      /*
      3  0 2 4 
      3  0 2 5 
      3  0 3 4 
      3  0 3 5 
      3  1 2 4 
      3  1 2 5 
      3  1 3 4 
      3  1 3 5
      */ 

      
      d_cf = new double[n_face];
      d_cv = new double[n_vertex];
      
      
    } //add other polyhedron types here
  else if(strcmp(temp_name,"cube")==0 || strcmp(temp_name,"P5")==0)
    {

      poly_name = "cube";
      n_vertex = 8;
      n_edge = 12;
      n_face = 6;
      n_fs = 1; //not really need this one, even for complicated shapes...
      max_face_vert = 4;

      face_num_fspecie = new int[n_fs];
      face_num_fspecie[0] = 6;

      vert_num_fspecie = new int[n_fs];
      vert_num_fspecie[0] = 4;
      
      //now initialize the data type
      Vertex = new double*[n_vertex];
      for(int i=0; i<n_vertex; i++)
	Vertex[i] = new double[SD];
      
      Edge = new int*[n_edge];
      for(int i=0; i<n_edge; i++)
	Edge[i] = new int[2];
      //now specify the index of vertices for each edge
      Edge[0][0] = 0; Edge[0][1] = 1;
      Edge[1][0] = 0; Edge[1][1] = 2;
      Edge[2][0] = 0; Edge[2][1] = 4;
      Edge[3][0] = 1; Edge[3][1] = 3;
      Edge[4][0] = 1; Edge[4][1] = 5;
      Edge[5][0] = 2; Edge[5][1] = 3;

      Edge[6][0] = 2; Edge[6][1] = 6;
      Edge[7][0] = 3; Edge[7][1] = 7;
      Edge[8][0] = 4; Edge[8][1] = 5;
      Edge[9][0] = 4; Edge[9][1] = 6;
      Edge[10][0] = 5; Edge[10][1] = 7;
      Edge[11][0] = 6; Edge[11][1] = 7;
      
      
      Face = new int*[n_face];
      Normal = new double*[n_face];  
      face_type = new int[n_face];
      for(int i=0; i<n_face; i++)
	{
	  Face[i] = new int[max_face_vert]; //make sure the list is feasible for all faces
	  Normal[i] = new double[SD];
	  face_type[i] = 0;
	}
      Face[0][0] = 0; Face[0][1] = 1; Face[0][2] = 3; Face[0][3] = 2;
      Face[1][0] = 0; Face[1][1] = 1; Face[1][2] = 5; Face[1][3] = 4;
      Face[2][0] = 0; Face[2][1] = 4; Face[2][2] = 6; Face[2][3] = 2;

      Face[3][0] = 1; Face[3][1] = 3; Face[3][2] = 7; Face[3][3] = 5;
      Face[4][0] = 3; Face[4][1] = 2; Face[4][2] = 6; Face[4][3] = 7;
      Face[5][0] = 4; Face[5][1] = 5; Face[5][2] = 7; Face[5][3] = 6;
      

      //the following face configuration is taken from AntiView
      //we re-order the points such as edged can be easily specified...
      /*
      4  4 5 6 7
      4  1 3 5 7
      4  2 3 6 7
      4  0 1 2 3
      4  0 1 4 5
      4  0 2 4 6
      */
      
      d_cf = new double[n_face];
      d_cv = new double[n_vertex];
      
      
    } //add other polyhedron types here
  else if(strcmp(temp_name,"user")==0 || strcmp(temp_name,"U")==0)
    {
      //for general shapes, need to play with the face_type

      cout<<" This USER SPECIFIED polyhedron shape: "<<endl;
      cout<<" Need to provide basic geometric parameters, edge-config, face-config in required format"<<endl;
      //cout<<" Still under construction..."<<endl;
      //exit(1);

       poly_name = "user";

      ifstream fin;
      fin.open("shape_param.txt");
      if(!fin)
	{
	  cout<<"Can not open file shape_param.txt! Re-check! "<<endl;
	  exit(1);
	}
      
      char dummy[100]; //for polyhedron name
      fin>>dummy; cout<<dummy<<endl;

      fin>>n_vertex; //cout<<"n_vertex = "<<n_vertex<<endl;
      fin>>n_edge; //cout<<"n_edge = "<<n_edge<<endl;
      fin>>n_face; //cout<<"n_face = "<<n_face<<endl;
      fin>>n_fs; //cout<<"n_fs = "<<n_fs<<endl;
      fin>>max_face_vert; //cout<<"max_face_vert = "<<max_face_vert<<endl;
      

      face_num_fspecie = new int[n_fs];
      for(int i=0; i<n_fs; i++)
	{
	  fin>>face_num_fspecie[i];
	  //cout<<"face_num_fspecie_"<<i<<" = "<<face_num_fspecie[i]<<endl;
	}
      
      
      vert_num_fspecie = new int[n_fs];
      for(int i=0; i<n_fs; i++)
	{
	  fin>>vert_num_fspecie[i];
	  //cout<<"vert_num_fspecie_"<<i<<" = "<<vert_num_fspecie[i]<<endl;
	}
      
      
      Vertex = new double*[n_vertex];
      for(int i=0; i<n_vertex; i++)
	Vertex[i] = new double[SD];
      
      Edge = new int*[n_edge];
      for(int i=0; i<n_edge; i++)
	Edge[i] = new int[2];
      
      
      Face = new int*[n_face];
      Normal = new double*[n_face];      
      for(int i=0; i<n_face; i++)
	{
	  Face[i] = new int[max_face_vert]; //make sure the list is feasible for all faces
	  Normal[i] = new double[SD];
	}


      d_cf = new double[n_face];
      d_cv = new double[n_vertex];

      face_type = new int[n_face];

      //now we read in the edge pairs
      for(int i=0; i<n_edge; i++)
	{
	  fin>>Edge[i][0]; fin>>Edge[i][1];
	  //cout<<"Edge_"<<i<<" = ("<<Edge[i][0]<<","<<Edge[i][1]<<")"<<endl;
	}
      
      //now we read in the faces...
      //IMPORTANT: THIS IS DIFFERENT THAN P1-P5, THE FIRST Number IS <<<<FACE TYPE>>>>!!!!!!
      for(int i=0; i<n_face; i++)
	{
	  fin>>face_type[i];
	  for(int j=0; j<vert_num_fspecie[face_type[i]]; j++)
	    fin>>Face[i][j];

	  /*
	  cout<<"Face_"<<i<<" = (";
	  for(int j=0; j<vert_num_fspecie[face_type[i]]; j++)
	    cout<<Face[i][j]<<" ";
	  cout<<")"<<endl;
	  */
	  
	}

      fin.close();
      
    }
  else
    {
      cout<<" Polyhedron "<<temp_name<<" has not been programmed yet!"<<endl;
      exit(1);
    }
  //here just prepare the data type, the values to be read in later
}



//exactly the same as Polyhedron(char*)
//for "new []" which cannot call constructor with parameters
void Polyhedron::PolyConstr(char* temp_name)
{
  if(strcmp(temp_name,"tetrahedron")==0 || strcmp(temp_name,"P1")==0)
    {

      poly_name = "tetrahedron";
      n_vertex = 4;
      n_edge = 6;
      n_face = 4;
      n_fs = 1; //not really need this one, even for complicated shapes...
      max_face_vert = 3;

      face_num_fspecie = new int[n_fs];
      face_num_fspecie[0] = 4;

      vert_num_fspecie = new int[n_fs];
      vert_num_fspecie[0] = 3;
      
      //now initialize the data type
      Vertex = new double*[n_vertex];
      for(int i=0; i<n_vertex; i++)
	Vertex[i] = new double[SD];
      
      Edge = new int*[n_edge];
      for(int i=0; i<n_edge; i++)
	Edge[i] = new int[2];
      //now specify the index of vertices for each edge
      Edge[0][0] = 0; Edge[0][1] = 1;
      Edge[1][0] = 0; Edge[1][1] = 2;
      Edge[2][0] = 0; Edge[2][1] = 3;
      Edge[3][0] = 1; Edge[3][1] = 2;
      Edge[4][0] = 1; Edge[4][1] = 3;
      Edge[5][0] = 2; Edge[5][1] = 3;
      
      
      Face = new int*[n_face];
      Normal = new double*[n_face]; 
      face_type = new int[n_face];
      for(int i=0; i<n_face; i++)
	{
	  Face[i] = new int[max_face_vert]; //make sure the list is feasible for all faces
	  Normal[i] = new double[SD];
	  face_type[i] = 0;
	}
      Face[0][0] = 0; Face[0][1] = 1; Face[0][2] = 2;
      Face[1][0] = 0; Face[1][1] = 1; Face[1][2] = 3;
      Face[2][0] = 0; Face[2][1] = 2; Face[2][2] = 3;
      Face[3][0] = 1; Face[3][1] = 2; Face[3][2] = 3;
      
      d_cf = new double[n_face];
      d_cv = new double[n_vertex];
      
      
    } //add other polyhedron types here
  else if(strcmp(temp_name,"icosahedron")==0 || strcmp(temp_name,"P2")==0)
    {

      poly_name = "icosahedron";
      n_vertex = 12;
      n_edge = 30;
      n_face = 20;
      n_fs = 1; //only needed when printing out the configuration...
      max_face_vert = 3;

      face_num_fspecie = new int[n_fs];
      face_num_fspecie[0] = 20;

      vert_num_fspecie = new int[n_fs];
      vert_num_fspecie[0] = 3;
      
      //now initialize the data type
      Vertex = new double*[n_vertex];
      for(int i=0; i<n_vertex; i++)
	Vertex[i] = new double[SD];
      
      Edge = new int*[n_edge];
      for(int i=0; i<n_edge; i++)
	Edge[i] = new int[2];
      //now specify the index of vertices for each edge, based on the face configuration
      Edge[0][0] = 0; Edge[0][1] = 1;
      Edge[1][0] = 0; Edge[1][1] = 4;
      Edge[2][0] = 0; Edge[2][1] = 5;
      Edge[3][0] = 0; Edge[3][1] = 8;
      Edge[4][0] = 0; Edge[4][1] = 9;

      Edge[5][0] = 1; Edge[5][1] = 4;
      Edge[6][0] = 1; Edge[6][1] = 5;
      Edge[7][0] = 1; Edge[7][1] = 10;
      Edge[8][0] = 1; Edge[8][1] = 11;
      Edge[9][0] = 2; Edge[9][1] = 3;

      Edge[10][0] = 2; Edge[10][1] = 6;
      Edge[11][0] = 2; Edge[11][1] = 7;
      Edge[12][0] = 2; Edge[12][1] = 8;
      Edge[13][0] = 2; Edge[13][1] = 9;
      Edge[14][0] = 3; Edge[14][1] = 6;
     
      Edge[15][0] = 3; Edge[15][1] = 7;
      Edge[16][0] = 3; Edge[16][1] = 10;
      Edge[17][0] = 3; Edge[17][1] = 11;
      Edge[18][0] = 4; Edge[18][1] = 6;
      Edge[19][0] = 4; Edge[19][1] = 8;

      Edge[20][0] = 4; Edge[20][1] = 10;
      Edge[21][0] = 5; Edge[21][1] = 7;
      Edge[22][0] = 5; Edge[22][1] = 9;
      Edge[23][0] = 5; Edge[23][1] = 11;
      Edge[24][0] = 6; Edge[24][1] = 8;

      Edge[25][0] = 6; Edge[25][1] = 10;
      Edge[26][0] = 7; Edge[26][1] = 9;
      Edge[27][0] = 7; Edge[27][1] = 11;
      Edge[28][0] = 8; Edge[28][1] = 9;
      Edge[29][0] = 10; Edge[29][1] = 11;


      
      
      Face = new int*[n_face];
      Normal = new double*[n_face]; 
      face_type = new int[n_face];
      for(int i=0; i<n_face; i++)
	{
	  Face[i] = new int[max_face_vert]; //make sure the list is feasible for all faces
	  Normal[i] = new double[SD];
	  face_type[i] = 0;
	}
      //the following list copied from AntiView output
      Face[0][0] = 0; Face[0][1] = 1; Face[0][2] = 4;
      Face[1][0] = 0; Face[1][1] = 1; Face[1][2] = 5;
      Face[2][0] = 0; Face[2][1] = 4; Face[2][2] = 8;
      Face[3][0] = 0; Face[3][1] = 5; Face[3][2] = 9;
      Face[4][0] = 0; Face[4][1] = 8; Face[4][2] = 9;
     
      Face[5][0] = 1; Face[5][1] = 4; Face[5][2] = 10;
      Face[6][0] = 1; Face[6][1] = 5; Face[6][2] = 11;
      Face[7][0] = 1; Face[7][1] = 10; Face[7][2] = 11;
      Face[8][0] = 2; Face[8][1] = 3; Face[8][2] = 6;
      Face[9][0] = 2; Face[9][1] = 3; Face[9][2] = 7;

      Face[10][0] = 2; Face[10][1] = 6; Face[10][2] = 8;
      Face[11][0] = 2; Face[11][1] = 7; Face[11][2] = 9;
      Face[12][0] = 2; Face[12][1] = 8; Face[12][2] = 9;
      Face[13][0] = 3; Face[13][1] = 6; Face[13][2] = 10;
      Face[14][0] = 3; Face[14][1] = 7; Face[14][2] = 11;

      Face[15][0] = 3; Face[15][1] = 10; Face[15][2] = 11;
      Face[16][0] = 4; Face[16][1] = 6; Face[16][2] = 8;
      Face[17][0] = 4; Face[17][1] = 6; Face[17][2] = 10;
      Face[18][0] = 5; Face[18][1] = 7; Face[18][2] = 9;
      Face[19][0] = 5; Face[19][1] = 7; Face[19][2] = 11;


      //the following is the AntiView output
      //3  0 1 4 
      //3  0 1 5 
      //3  0 4 8 
      //3  0 5 9 
      //3  0 8 9 

      //3  1 4 10 
      //3  1 5 11 
      //3  1 10 11 
      //3  2 3 6 
      //3  2 3 7 

      //3  2 6 8 
      //3  2 7 9 
      //3  2 8 9 
      //3  3 6 10 
      //3  3 7 11
 
      //3  3 10 11 
      //3  4 6 8 
      //3  4 6 10 
      //3  5 7 9 
      //3  5 7 11 
      
      d_cf = new double[n_face];
      d_cv = new double[n_vertex];
      
      
    } //add other polyhedron types here
  else if(strcmp(temp_name,"dodecahedron")==0 || strcmp(temp_name,"P3")==0)
    {

      poly_name = "dodecahedron";
      n_vertex = 20;
      n_edge = 30;
      n_face = 12;
      n_fs = 1; //not really need this one, even for complicated shapes...
      max_face_vert = 5;

      face_num_fspecie = new int[n_fs];
      face_num_fspecie[0] = 12;

      vert_num_fspecie = new int[n_fs];
      vert_num_fspecie[0] = 5;
      
      //now initialize the data type
      Vertex = new double*[n_vertex];
      for(int i=0; i<n_vertex; i++)
	Vertex[i] = new double[SD];
      
      Edge = new int*[n_edge];
      for(int i=0; i<n_edge; i++)
	Edge[i] = new int[2];
      //now specify the index of vertices for each edge     
      Edge[0][0] = 0; Edge[0][1] = 8;
      Edge[1][0] = 0; Edge[1][1] = 12;
      Edge[2][0] = 0; Edge[2][1] = 16;
      Edge[3][0] = 1; Edge[3][1] = 9;
      Edge[4][0] = 1; Edge[4][1] = 12;
      Edge[5][0] = 1; Edge[5][1] = 17;

      Edge[6][0] = 2; Edge[6][1] = 10;
      Edge[7][0] = 2; Edge[7][1] = 13;
      Edge[8][0] = 2; Edge[8][1] = 16;
      Edge[9][0] = 3; Edge[9][1] = 11;
      Edge[10][0] = 3; Edge[10][1] = 13;
      Edge[11][0] = 3; Edge[11][1] = 17;

      Edge[12][0] = 4; Edge[12][1] = 8;
      Edge[13][0] = 4; Edge[13][1] = 14;
      Edge[14][0] = 4; Edge[14][1] = 18;
      Edge[15][0] = 5; Edge[15][1] = 9;
      Edge[16][0] = 5; Edge[16][1] = 14;
      Edge[17][0] = 5; Edge[17][1] = 19;

      Edge[18][0] = 6; Edge[18][1] = 10;
      Edge[19][0] = 6; Edge[19][1] = 15;
      Edge[20][0] = 6; Edge[20][1] = 18;
      Edge[21][0] = 7; Edge[21][1] = 11;
      Edge[22][0] = 7; Edge[22][1] = 15;
      Edge[23][0] = 7; Edge[23][1] = 19;

      Edge[24][0] = 8; Edge[24][1] = 10;
      Edge[25][0] = 9; Edge[25][1] = 11;
      Edge[26][0] = 12; Edge[26][1] = 14;
      Edge[27][0] = 13; Edge[27][1] = 15;
      Edge[28][0] = 16; Edge[28][1] = 17;
      Edge[29][0] = 18; Edge[29][1] = 19;
      
      
      Face = new int*[n_face];
      Normal = new double*[n_face];  
      face_type = new int[n_face];    
      for(int i=0; i<n_face; i++)
	{
	  Face[i] = new int[max_face_vert]; //make sure the list is feasible for all faces
	  Normal[i] = new double[SD];
	  face_type[i] = 0;
	}
      Face[0][0] = 3; Face[0][1] = 11; Face[0][2] = 7; Face[0][3] = 15; Face[0][4] = 13;
      Face[1][0] = 6; Face[1][1] = 15; Face[1][2] = 7; Face[1][3] = 19; Face[1][4] = 18;
      Face[2][0] = 5; Face[2][1] = 9; Face[2][2] = 11; Face[2][3] = 7; Face[2][4] = 19;
      Face[3][0] = 1; Face[3][1] = 9; Face[3][2] = 5; Face[3][3] = 14; Face[3][4] = 12;
      Face[4][0] = 0; Face[4][1] = 12; Face[4][2] = 1; Face[4][3] = 17; Face[4][4] = 16;
      Face[5][0] = 1; Face[5][1] = 9; Face[5][2] = 11; Face[5][3] = 3; Face[5][4] = 17;

      Face[6][0] = 2; Face[6][1] = 10; Face[6][2] = 6; Face[6][3] = 15; Face[6][4] = 13;
      Face[7][0] = 2; Face[7][1] = 13; Face[7][2] = 3; Face[7][3] = 17; Face[7][4] = 16;
      Face[8][0] = 0; Face[8][1] = 8; Face[8][2] = 10; Face[8][3] = 2; Face[8][4] = 16;
      Face[9][0] = 0; Face[9][1] = 8; Face[9][2] = 4; Face[9][3] = 14; Face[9][4] = 12;
      Face[10][0] = 4; Face[10][1] = 14; Face[10][2] = 5; Face[10][3] = 19; Face[10][4] = 18;
      Face[11][0] = 4; Face[11][1] = 8; Face[11][2] = 10; Face[11][3] = 6; Face[11][4] = 18;

      //the following is the face configuration from AntiView
      /*
	5  3 11 7 15 13
	5  6 15 7 19 18
	5  5 9 11 7 19
	5  1 9 5 14 12
	5  0 12 1 17 16
	5  1 9 11 3 17

	5  2 10 6 15 13
	5  2 13 3 17 16
	5  0 8 10 2 16
	5  0 8 4 14 12
	5  4 14 5 19 18
	5  4 8 10 6 18
      */



      d_cf = new double[n_face];
      d_cv = new double[n_vertex];
      
      
    } //add other polyhedron types here
  else if(strcmp(temp_name,"octahedron")==0 || strcmp(temp_name,"P4")==0)
    {

      poly_name = "octahedron";
      n_vertex = 6;
      n_edge = 12;
      n_face = 8;
      n_fs = 1; //not really need this one, even for complicated shapes...
      max_face_vert = 3;

      face_num_fspecie = new int[n_fs];
      face_num_fspecie[0] = 8;

      vert_num_fspecie = new int[n_fs];
      vert_num_fspecie[0] = 3;
      
      //now initialize the data type
      Vertex = new double*[n_vertex];
      for(int i=0; i<n_vertex; i++)
	Vertex[i] = new double[SD];
      
      Edge = new int*[n_edge];
      for(int i=0; i<n_edge; i++)
	Edge[i] = new int[2];
      //now specify the index of vertices for each edge
      Edge[0][0] = 0; Edge[0][1] = 2;
      Edge[1][0] = 0; Edge[1][1] = 3;
      Edge[2][0] = 0; Edge[2][1] = 4;
      Edge[3][0] = 0; Edge[3][1] = 5;

      Edge[4][0] = 1; Edge[4][1] = 2;
      Edge[5][0] = 1; Edge[5][1] = 3;
      Edge[6][0] = 1; Edge[6][1] = 4;
      Edge[7][0] = 1; Edge[7][1] = 5;

      Edge[8][0] = 2; Edge[8][1] = 4;
      Edge[9][0] = 2; Edge[9][1] = 5;
      Edge[10][0] = 3; Edge[10][1] = 4;
      Edge[11][0] = 3; Edge[11][1] = 5;
      
      
      Face = new int*[n_face];
      Normal = new double*[n_face]; 
      face_type = new int[n_face];
      for(int i=0; i<n_face; i++)
	{
	  Face[i] = new int[max_face_vert]; //make sure the list is feasible for all faces
	  Normal[i] = new double[SD];
	  face_type[i] = 0;
	}
      Face[0][0] = 0; Face[0][1] = 2; Face[0][2] = 4;
      Face[1][0] = 0; Face[1][1] = 2; Face[1][2] = 5;
      Face[2][0] = 0; Face[2][1] = 3; Face[2][2] = 4;
      Face[3][0] = 0; Face[3][1] = 3; Face[3][2] = 5;

      Face[4][0] = 1; Face[4][1] = 2; Face[4][2] = 4;
      Face[5][0] = 1; Face[5][1] = 2; Face[5][2] = 5;
      Face[6][0] = 1; Face[6][1] = 3; Face[6][2] = 4;
      Face[7][0] = 1; Face[7][1] = 3; Face[7][2] = 5;


      //face configuration taken from AntiView
      /*
      3  0 2 4 
      3  0 2 5 
      3  0 3 4 
      3  0 3 5 
      3  1 2 4 
      3  1 2 5 
      3  1 3 4 
      3  1 3 5
      */ 

      
      d_cf = new double[n_face];
      d_cv = new double[n_vertex];
      
      
    } //add other polyhedron types here
  else if(strcmp(temp_name,"cube")==0 || strcmp(temp_name,"P5")==0)
    {

      poly_name = "cube";
      n_vertex = 8;
      n_edge = 12;
      n_face = 6;
      n_fs = 1; //not really need this one, even for complicated shapes...
      max_face_vert = 4;

      face_num_fspecie = new int[n_fs];
      face_num_fspecie[0] = 6;

      vert_num_fspecie = new int[n_fs];
      vert_num_fspecie[0] = 4;
      
      //now initialize the data type
      Vertex = new double*[n_vertex];
      for(int i=0; i<n_vertex; i++)
	Vertex[i] = new double[SD];
      
      Edge = new int*[n_edge];
      for(int i=0; i<n_edge; i++)
	Edge[i] = new int[2];
      //now specify the index of vertices for each edge
      Edge[0][0] = 0; Edge[0][1] = 1;
      Edge[1][0] = 0; Edge[1][1] = 2;
      Edge[2][0] = 0; Edge[2][1] = 4;
      Edge[3][0] = 1; Edge[3][1] = 3;
      Edge[4][0] = 1; Edge[4][1] = 5;
      Edge[5][0] = 2; Edge[5][1] = 3;

      Edge[6][0] = 2; Edge[6][1] = 6;
      Edge[7][0] = 3; Edge[7][1] = 7;
      Edge[8][0] = 4; Edge[8][1] = 5;
      Edge[9][0] = 4; Edge[9][1] = 6;
      Edge[10][0] = 5; Edge[10][1] = 7;
      Edge[11][0] = 6; Edge[11][1] = 7;
      
      
      Face = new int*[n_face];
      Normal = new double*[n_face]; 
      face_type = new int[n_face];
      for(int i=0; i<n_face; i++)
	{
	  Face[i] = new int[max_face_vert]; //make sure the list is feasible for all faces
	  Normal[i] = new double[SD];
	  face_type[i] = 0;
	}
      Face[0][0] = 0; Face[0][1] = 1; Face[0][2] = 3; Face[0][3] = 2;
      Face[1][0] = 0; Face[1][1] = 1; Face[1][2] = 5; Face[1][3] = 4;
      Face[2][0] = 0; Face[2][1] = 4; Face[2][2] = 6; Face[2][3] = 2;

      Face[3][0] = 1; Face[3][1] = 3; Face[3][2] = 7; Face[3][3] = 5;
      Face[4][0] = 3; Face[4][1] = 2; Face[4][2] = 6; Face[4][3] = 7;
      Face[5][0] = 4; Face[5][1] = 5; Face[5][2] = 7; Face[5][3] = 6;
      

      //the following face configuration is taken from AntiView
      //we re-order the points such as edged can be easily specified...
      /*
      4  4 5 6 7
      4  1 3 5 7
      4  2 3 6 7
      4  0 1 2 3
      4  0 1 4 5
      4  0 2 4 6
      */
      
      d_cf = new double[n_face];
      d_cv = new double[n_vertex];
      
      
    } //add other polyhedron types here
  else if(strcmp(temp_name,"user")==0 || strcmp(temp_name,"U")==0)
    {
      //for general shapes, need to play with the face_type

      //cout<<" This USER SPECIFIED polyhedron shape: "<<endl;
      //cout<<" Need to provide basic geometric parameters, edge-config, face-config in required format"<<endl;
      //cout<<" Still under construction..."<<endl;
      //exit(1);

      poly_name = "user";

      ifstream fin;
      fin.open("shape_param.txt");
      if(!fin)
	{
	  cout<<"Can not open file shape_param.txt! Re-check! "<<endl;
	  exit(1);
	}
      
      char dummy[100]; //for polyhedron name
      fin>>dummy; cout<<dummy<<endl;

      fin>>n_vertex; //cout<<"n_vertex = "<<n_vertex<<endl;
      fin>>n_edge; //cout<<"n_edge = "<<n_edge<<endl;
      fin>>n_face; //cout<<"n_face = "<<n_face<<endl;
      fin>>n_fs; //cout<<"n_fs = "<<n_fs<<endl;
      fin>>max_face_vert; //cout<<"max_face_vert = "<<max_face_vert<<endl;
      

      face_num_fspecie = new int[n_fs];
      for(int i=0; i<n_fs; i++)
	{
	  fin>>face_num_fspecie[i];
	  //cout<<"face_num_fspecie_"<<i<<" = "<<face_num_fspecie[i]<<endl;
	}
      
      
      vert_num_fspecie = new int[n_fs];
      for(int i=0; i<n_fs; i++)
	{
	  fin>>vert_num_fspecie[i];
	  //cout<<"vert_num_fspecie_"<<i<<" = "<<vert_num_fspecie[i]<<endl;
	}
      //cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
      
      Vertex = new double*[n_vertex];
      for(int i=0; i<n_vertex; i++)
	Vertex[i] = new double[SD];
      
      Edge = new int*[n_edge];
      for(int i=0; i<n_edge; i++)
	Edge[i] = new int[2];
      
      
      Face = new int*[n_face];
      Normal = new double*[n_face];      
      for(int i=0; i<n_face; i++)
	{
	  Face[i] = new int[max_face_vert]; //make sure the list is feasible for all faces
	  Normal[i] = new double[SD];
	}


      face_type = new int[n_face];

      d_cf = new double[n_face];
      d_cv = new double[n_vertex];

      //now we read in the edge pairs
      for(int i=0; i<n_edge; i++)
	{
	  fin>>Edge[i][0]; fin>>Edge[i][1];
	  //cout<<"Edge_"<<i<<" = ("<<Edge[i][0]<<","<<Edge[i][1]<<")"<<endl;
	}
      
      //now we read in the faces...
      //IMPORTANT: THIS IS DIFFERENT THAN P1-P5, THE FIRST Number IS <<<<FACE TYPE>>>>!!!!!!
      for(int i=0; i<n_face; i++)
	{
	  fin>>face_type[i];
	  for(int j=0; j<vert_num_fspecie[face_type[i]]; j++)
	    fin>>Face[i][j];

	  /*
	  cout<<"Face_"<<i<<" = (";
	  for(int j=0; j<vert_num_fspecie[face_type[i]]; j++)
	    cout<<Face[i][j]<<" ";
	  cout<<")"<<endl;
	  */	  
	}

      fin.close();

      //cout<<"Here!"<<endl;
      
    }
  else
    {
      cout<<" Polyhedron "<<temp_name<<" has not been programmed yet!"<<endl;
      exit(1);
    }

  //here just prepare the data type, the values to be read in later
}



void Polyhedron::ShiftVert()
{
  double Centroid[SD] = {0,0,0};

  for(int i=0; i<n_vertex; i++)
    for(int j=0; j<SD; j++)
      Centroid[j] += Vertex[i][j];

  for(int i=0; i<n_vertex; i++)
    for(int j=0; j<SD; j++)
      Vertex[i][j] = Vertex[i][j] - Centroid[j]/(double)n_vertex;
}

void Polyhedron::GetNormal_All()
{
  double v1[SD], v2[SD];
  double temp_norm[SD];

  //cout<<"n_face = "<<n_face<<endl;
  
  //loop over each face ...
  for(int i=0; i<n_face; i++)
    {
      for(int k=0; k<SD; k++) //only involve first three vertices, should work for all shapes
	{
	  v1[k] = Vertex[Face[i][0]][k] - Vertex[Face[i][1]][k];
	  v2[k] = Vertex[Face[i][0]][k] - Vertex[Face[i][2]][k];
	}
      
      Comput.GetCrossProduct(v1, v2, temp_norm);

      //cout<<temp_norm[0]<<"\t"<<temp_norm[1]<<"\t"<<temp_norm[2]<<endl;
      
      Comput.Normalize(temp_norm);

      //Assuming the centroid is at 0, this is the inner direction....
      if(Comput.GetInnerProduct(Vertex[Face[i][0]], temp_norm)<0) 
	{
	  for(int k=0; k<SD; k++)
	    Normal[i][k] = -temp_norm[k];
	}
      else //make sure it is the outward normal...
	{
	  for(int k=0; k<SD; k++)
	    Normal[i][k] = temp_norm[k];
	}
      

    }
}

void Polyhedron::GetNormal_F(int face_index)
{
  double v1[SD], v2[SD];
  double temp_norm[SD];

  //cout<<"n_face = "<<n_face<<endl;
  
  //only for the face specified by face_index

  for(int k=0; k<SD; k++) //only involve first three vertices, should work for all shapes
    {
      v1[k] = Vertex[Face[face_index][0]][k] - Vertex[Face[face_index][1]][k];
      v2[k] = Vertex[Face[face_index][0]][k] - Vertex[Face[face_index][2]][k];
    }
      
  Comput.GetCrossProduct(v1, v2, temp_norm);

  //cout<<temp_norm[0]<<"\t"<<temp_norm[1]<<"\t"<<temp_norm[2]<<endl;
      
  Comput.Normalize(temp_norm);

  //Assuming the centroid is at 0, this is the inner direction....
  if(Comput.GetInnerProduct(Vertex[Face[face_index][0]], temp_norm)<0) 
    {
      for(int k=0; k<SD; k++)
	Normal[face_index][k] = -temp_norm[k];
    }
  else //make sure it is the outward normal...
    {
      for(int k=0; k<SD; k++)
	Normal[face_index][k] = temp_norm[k];
    }
  
  
  
}


void Polyhedron::GetCentDist()
{ 
  
  GetNormal_All(); 

  //now get the center2face distance
  for(int i=0; i<n_face; i++)
    d_cf[i] = fabs(Comput.GetInnerProduct(Vertex[Face[i][0]], Normal[i]));

  //now get the center2vertex distance
  for(int i=0; i<n_vertex; i++)
    d_cv[i] = Comput.GetLength(Vertex[i]);
 
} 

//this only works for the same type of polyhedron, i.e., copy a tetrah from a terah
// one can not copy a tetrah from an octahedron
// thus, only the coordinate data need to be copied, others e.g., n_vertex etc are initialized when poly type is declared
void Polyhedron::CopyPoly(Polyhedron &B)
{
  /*
  n_vertex = B.n_vertex;
  n_edge = B.n_edge;
  n_face = B.n_face;
  n_fs = B.n_fs; //not really need this one, even for complicated shapes...
  max_face_vert = B.max_face_vert;

  for(int i=0; i<n_fs; i++)
    {
      face_num_fspecie[i] = B.face_num_fspecie[i];
      vert_num_fspecie[i] = B.vert_num_fspecie[i];
    }
  */

  for(int i=0; i<n_vertex; i++)
    for(int k=0; k<SD; k++)
      Vertex[i][k] = B.Vertex[i][k];
  
  /*
  for(int i=0; i<n_edge; i++)
    for(int k=0; k<2; k++)
      Edge[i][k] = B.Edge[i][k];
  */
  
  for(int i=0; i<n_face; i++)
    {
      //for(int k=0; k<max_face_vert; k++)
      //Face[i][k] = B.Face[i][k];
      
      for(int k=0; k<SD; k++)
	Normal[i][k] = B.Normal[i][k];
    }
  
  for(int i=0; i<n_face; i++)
    d_cf[i] = B.d_cf[i];
  
  for(int i=0; i<n_vertex; i++)
    d_cv[i] = B.d_cv[i];
        
}


//reload operator =
Polyhedron& Polyhedron::operator = (Polyhedron &B)
{
  CopyPoly(B);

  return *this;
}


//for all complicated shapes, divide the shape into pyrimad with bases being the faces of the polyhedron
double Polyhedron::GetVol()
{
  if(strcmp(poly_name,"tetrahedron")==0 || strcmp(poly_name,"P1")==0)
    {
      double vect1[SD], vect2[SD], vect3[SD], vect4[SD];

      for(int i=0; i<SD; i++)
	{
	  vect1[i] = Vertex[0][i] - Vertex[3][i];
	  vect2[i] = Vertex[1][i] - Vertex[3][i];
	  vect3[i] = Vertex[2][i] - Vertex[3][i];
	}
      
      Comput.GetCrossProduct(vect2, vect3, vect4);
      
      double vol = fabs(Comput.GetInnerProduct(vect1, vect4))/6.0;
      
      return vol;
    }
  else if(strcmp(poly_name,"icosahedron")==0 || strcmp(poly_name,"P2")==0)
    {
      double vect1[SD], vect2[SD], vect3[SD], vect4[SD];
      double vol = 0.0;

      for(int f=0; f<n_face; f++)
	{
	  for(int i=0; i<SD; i++)
	    {
	      vect1[i] = Vertex[Face[f][0]][i] - Vertex[Face[f][1]][i];
	      vect2[i] = Vertex[Face[f][0]][i] - Vertex[Face[f][2]][i];
	      vect3[i] = Vertex[Face[f][0]][i];
	    }


	  Comput.GetCrossProduct(vect1, vect2, vect4);

	  vol += fabs(Comput.GetInnerProduct(vect3, vect4))/6.0;
	}
      
      return vol;
    }
  else if(strcmp(poly_name,"dodecahedron")==0 || strcmp(poly_name,"P3")==0) //this volume computing method is general...
    {
      double vect1[SD], vect2[SD], vect3[SD], vect4[SD];
      double vol = 0.0;

      for(int f=0; f<n_face; f++)
	{
	  double face_cent[SD]={0.0,0.0,0.0};

	  for(int i=0; i<SD; i++)
	    {
	      for(int j=0; j<vert_num_fspecie[0]; j++)
		face_cent[i] += Vertex[Face[f][j]][i];

	      face_cent[i] = face_cent[i]/(double)vert_num_fspecie[0];		 
	    }

	  //divide the pentagon into 5 small triangles and compute volume separately...
	  for(int j=0; j<vert_num_fspecie[0]; j++)
	    {
	      int k=j+1; //the second vertex for an edge
	      if(k>=vert_num_fspecie[0]) k = k - vert_num_fspecie[0];

	      for(int i=0; i<SD; i++)
		{
		  vect1[i] = face_cent[i] - Vertex[Face[f][j]][i];
		  vect2[i] = face_cent[i] - Vertex[Face[f][k]][i];
		}

	      Comput.GetCrossProduct(vect1, vect2, vect4);

	      vol += fabs(Comput.GetInnerProduct(face_cent, vect4))/6.0;

	    }
	 
	}
      
      return vol;
    }
  else if(strcmp(poly_name,"octahedron")==0 || strcmp(poly_name,"P4")==0)
    {
      double vect1[SD], vect2[SD], vect3[SD], vect4[SD];
      double vol = 0.0;

      for(int f=0; f<n_face; f++)
	{
	  double face_cent[SD]={0.0,0.0,0.0};

	  for(int i=0; i<SD; i++)
	    {
	      for(int j=0; j<vert_num_fspecie[0]; j++)
		face_cent[i] += Vertex[Face[f][j]][i];

	      face_cent[i] = face_cent[i]/(double)vert_num_fspecie[0];		 
	    }

	  //divide the polygon into vert_num_fspecie small triangles and compute volume separately...
	  for(int j=0; j<vert_num_fspecie[0]; j++)
	    {
	      int k=j+1; //the second vertex for an edge
	      if(k>=vert_num_fspecie[0]) k = k - vert_num_fspecie[0];

	      for(int i=0; i<SD; i++)
		{
		  vect1[i] = face_cent[i] - Vertex[Face[f][j]][i];
		  vect2[i] = face_cent[i] - Vertex[Face[f][k]][i];
		}

	      Comput.GetCrossProduct(vect1, vect2, vect4);

	      vol += fabs(Comput.GetInnerProduct(face_cent, vect4))/6.0;

	    }

	  //cout<<"f_center ["<<f<<"] = "<<face_cent[0]<<" "<<face_cent[1]<<" "<<face_cent[2]<<endl;
	  //cout<<"vol_f = "<<vol<<endl;
	 
	}
      
      return vol;
    }
  else if(strcmp(poly_name,"cube")==0 || strcmp(poly_name,"P5")==0)
    {
      double vect1[SD], vect2[SD], vect3[SD], vect4[SD];
      double vol = 0.0;

      for(int f=0; f<n_face; f++)
	{
	  double face_cent[SD]={0.0,0.0,0.0};

	  for(int i=0; i<SD; i++)
	    {
	      for(int j=0; j<vert_num_fspecie[0]; j++)
		face_cent[i] += Vertex[Face[f][j]][i];

	      face_cent[i] = face_cent[i]/(double)vert_num_fspecie[0];		 
	    }

	  //divide the polygon into vert_num_fspecie small triangles and compute volume separately...
	  for(int j=0; j<vert_num_fspecie[0]; j++)
	    {
	      int k=j+1; //the second vertex for an edge
	      if(k>=vert_num_fspecie[0]) k = k - vert_num_fspecie[0];

	      for(int i=0; i<SD; i++)
		{
		  vect1[i] = face_cent[i] - Vertex[Face[f][j]][i];
		  vect2[i] = face_cent[i] - Vertex[Face[f][k]][i];
		}

	      Comput.GetCrossProduct(vect1, vect2, vect4);

	      vol += fabs(Comput.GetInnerProduct(face_cent, vect4))/6.0;

	    }
	  
	  //cout<<"f_center ["<<f<<"] = "<<face_cent[0]<<" "<<face_cent[1]<<" "<<face_cent[2]<<endl;
	  //cout<<"vol_f_inner = "<<fabs(Comput.GetInnerProduct(face_cent, vect4))/6.0<<endl;
	  //cout<<"vol_f = "<<vol<<endl;
	 
	}
      
      //cout<<"vol_cube = "<<vol<<endl;
      return vol;
    }
  else if(strcmp(poly_name,"user")==0 || strcmp(poly_name,"U")==0)
    {
      double vect1[SD], vect2[SD], vect3[SD], vect4[SD];
      double vol = 0.0;

      for(int f=0; f<n_face; f++)
	{
	  
	  double face_cent[SD]={0.0,0.0,0.0};

	  for(int i=0; i<SD; i++)
	    {
	      for(int j=0; j<vert_num_fspecie[face_type[f]]; j++)
		face_cent[i] += Vertex[Face[f][j]][i];

	      face_cent[i] = face_cent[i]/(double)vert_num_fspecie[face_type[f]];		 
	    }

	  //divide the polygon into vert_num_fspecie small triangles and compute volume separately...
	  for(int j=0; j<vert_num_fspecie[face_type[f]]; j++)
	    {
	      int k=j+1; //the second vertex for an edge
	      if(k>=vert_num_fspecie[face_type[f]]) k = k - vert_num_fspecie[face_type[f]];

	      for(int i=0; i<SD; i++)
		{
		  vect1[i] = face_cent[i] - Vertex[Face[f][j]][i];
		  vect2[i] = face_cent[i] - Vertex[Face[f][k]][i];
		}

	      Comput.GetCrossProduct(vect1, vect2, vect4);

	      vol += fabs(Comput.GetInnerProduct(face_cent, vect4))/6.0;

	    }
	  
	  //cout<<"f_center ["<<f<<"] = "<<face_cent[0]<<" "<<face_cent[1]<<" "<<face_cent[2]<<endl;
	  //cout<<"vol_f_inner = "<<fabs(Comput.GetInnerProduct(face_cent, vect4))/6.0<<endl;
	  //cout<<"vol_f = "<<vol<<endl;
	 
	}
      
      //cout<<"vol_cube = "<<vol<<endl;
      return vol;
    }
  else
    {
      cout<<"The volume for the shape is not programmed! Return number denstiy only!"<<endl;
      return 1.0;
    }
}



//free all the memory space when the object is destroyed
Polyhedron::~Polyhedron()
{
  for(int i=0; i<n_vertex; i++)
    delete [] Vertex[i];
  delete [] Vertex;

  for(int i=0; i<n_edge; i++)
    delete [] Edge[i];
  delete [] Edge;

  for(int i=0; i<n_face; i++)
    {
      delete [] Face[i];
      delete [] Normal[i];
    }
  delete [] Face;
  delete [] Normal;

  delete [] face_num_fspecie;
  delete [] vert_num_fspecie;

  delete [] d_cf;
  delete [] d_cv;

  //cout<<"Destructor for Polyhedron is called!"<<endl;
}
