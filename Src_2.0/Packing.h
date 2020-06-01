//the data members include the vertices and positions of the packings
//    the lattice vectors, charactersitics of the shape etc.
//the function members include read packing, print packing and overlap check

#include "Cells.h"

class Packing{

 public:
  int N; //the number of particles
  double d0;  
  //for random shapes, this is a charactersitic length
  //   useful when determining the cell numbers
  //   normally choose to be twice the largest d_cv, computed when read in the config.

  double VParticle; //the total volume of all particles, which is a constant for the packing

  double** Lambda; //the lattice vectors...

  double** CenterE; // the centers in Eculdean coordinates...
  double** CenterL; // the centers in Relative coordinates...

  Polyhedron* PackPoly; //the collection of polyhedra, only provide vertice info

  Geometry Comput;

  int*** NNL; //near-neighbor list, use normal 2+1 D array, keep a counter for simplicity, the 3D entry store cell index
  int* NNL_counter; //the NNL counter

  int*** SEP_List; //the separation information list...
  // -1 -1 -1 if m and n are not near neighbor
  // 0  0/1 f_m/f_n; if m, n are neighbors; if separated by a face of m, then 0, 0, f_m; if spearted by a face of n then 0, 1, f_n 
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //the function member of class Packing
  void Rescale(double); //rescale the packing
  void Check_Iconfig(double, double); //check inital configuration rescale the packing if necessary
  Packing(char*); //specify the name, can be generalized to hybrid shape packings
  void GetGlobalPosition();

  double MinDis(double [], double [], double [], int, int, int);
  //the arguments are: (1) relative coordinates of centroid for particle m
  //                   (2) relative coordinates of centroid for particle n
  //                   (3) Euclidean displacement vector pointing from n to m
  //                   (4)-(6) index of box the image of m is in
  double MinDisII(int, int, double [], int &, int &, int &);
  //this is particularly for SEP_List

  double GetMin(double [], int); //passing an array and its size
  double GetMax(double [], int);

  int CheckOverlap(Polyhedron &, Polyhedron &, double [], double [], int indexI, int indexJ, int indexK);
  //the arguments are: (1) object particle m
  //                   (2) object particle n
  //                   (3) relative coordinates of centroid for particle m
  //                   (4) relative coordinates of centroid for particle n
  //                   (5)-(7) index of box the image of m is in

  int GlobalOverlapCheck(int); //this is the brute force check, int is only a flag, print out the overlap pair
  int GlobalOverlapCheck(double); //this is the brute force check, double is only a flag, don't print out the pair
  int GlobalOverlapCheck(Cells &); //check using the cell list

  int CheckOverlap(Polyhedron &, Polyhedron &, double [], double [], int indexI, int indexJ, int indexK, int m, int n);
  //the arguments are: (1) object particle m
  //                   (2) object particle n
  //                   (3) relative coordinates of centroid for particle m
  //                   (4) relative coordinates of centroid for particle n
  //                   (5)-(7) index of box the image of m is in
  //                   (8),(9) index of the particles, this is for use of SEP_List 



  void GetNNL(double); //the parameter is a cut distance...
  void UpdateNNL(double); //the parameters are cut distance ...

  void GetSEP_List(char*);


  void GetVParticle();
  double GetDensity();
  void PrintDensity();
  void PrintDensity(double);

  void PrintOverlapPair(int, int); //specify the index of a pair of polyhedra
  void PrintPacking(); //print the packing 
  //void PlotPacking(); //print an AntiView file for the final packing
};

