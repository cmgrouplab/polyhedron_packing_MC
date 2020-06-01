//use a similar data structure as Packing for coding simplicity...
//use overlap-checking functions to check contacts:
//     (1) obtain the minimal distance (p2f or e2e)
//     (2) if p2f, determine the dist for other points within the same face to see whether possible e2f or f2f contact...
//     (3) if e2e, directly report a e2e contact ...


//#include "Packing.h"
//whether this is needed depends on the order of other directives...


class Contact{

  public:
  int N; //the number of particles
  double d0;  
  //for random shapes, this is a charactersitic length
  //   useful when determining the cell numbers
  //   normally choose to be twice the largest d_cv, computed when read in the config.

  double TOL; //the numerical tolerance determining whether two particles contact or not
  double radial_cutoff; //this is for determing the possible near neighbors, different for different shapes 
  double normal_cutoff; //this is for determing f2f contact based on the inner product of face normals
  double f2f_cutoff; //for determining face related contact based the real separation distance, i.e., f2f_cutoff*edge_len
  double e2f_cutoff;
  double v2f_cutoff;
  double edge_cutoff; //this is for determing e2e contact based on center2center distance between two edges
  //for simplicity, we use a circular disk and its effective radius to approximate the real faces
  //a more sophiscated search based on the exact shape can be established

  double** Lambda; //the lattice vectors...

  double** CenterE; // the centers in Eculdean coordinates...
  double** CenterL; // the centers in Relative coordinates...

  Polyhedron* PackPoly; //the collection of polyhedra, only provide vertice info

  Geometry Comput;

  int** Contact_Array; //a N by N matrix indicating the contact nature between any two particles
                       //if not contact, the entry value is -1
                       //f2f ->1; e2f ->2, p2f ->3, e2e -> 4

  int n_f2f, n_e2f, n_p2f, n_e2e; //the total number of contacts for each contact type...
  int* NP_f2f; int* NP_e2f; int* NP_p2f; int* NP_e2e; //the number of conacts of a particular type for each particle...

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //the function member of class Contact
  Contact(char* temp_name, Packing &); 
                      //the constructor (1) name of polyhedron (2) pass a packing
                      //from the passed packing, initialize all necessary data types
                      //initialize other parameters (TOL) and data

  void GetGlobalPosition();

  double MinDis(double [], double [], double [], int, int, int);
  //the arguments are: (1) relative coordinates of centroid for particle m
  //                   (2) relative coordinates of centroid for particle n
  //                   (3) Euclidean displacement vector pointing from n to m
  //                   (4)-(6) index of box the image of m is in


  double GetMin(double [], int); //passing an array and its size
  double GetMax(double [], int);

  int GetContact(Polyhedron &, Polyhedron &, double [], double [], int indexI, int indexJ, int indexK);
  //the arguments are: (1) object particle m
  //                   (2) object particle n
  //                   (3) relative coordinates of centroid for particle m
  //                   (4) relative coordinates of centroid for particle n
  //                   (5)-(7) index of box the minimal image of m is in

  void GlobalGetContact(); //this is the brute force check, int is only a flag
  //int GlobalGetContact(Cells &); //check using the cell list

  void Print_Stat();
};
