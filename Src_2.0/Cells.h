//for the cell method
//the data members include cell numbers along each direction and 
//    the list of particle index in each cell


class Cells{

 public:
  
  double rho_cell; //the density of last update point

  //the number of cells along each direction, determined automatically later
  //    these are also updated adaptively
  int LXcell; 
  int LYcell;
  int LZcell;
  int Ncell; //the total number of cells
  
  node** CellList; //the cell list containing the particles in the cell...
  //the corner locations are the index...
  //implemented as a one-dimensional array

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //the function members

  Cells();
  int NumMax(int, int); //get the max of two integers
  void DelCellList(); //empty the cell list

  void GetCellList(int, double &, double** &, double**&, int &, double &);
  //the arguments are (1) a flag indicating whether the list need to be destoryed
  //                  (2) a characteristics length of the polyhedron
  //                  (3) the lattice vector
  //                  (4) the relative coordinates of the centroids
  //                  (5) the total number of particles
  //                  (6) the packing density at this update point

  void UpdateCellList(int, double [], double []);
  //the arguments are (1) the index of the moved polyhedron
  //                  (2) the old relative coordinates of the centroid
  //                  (3) the updated relative coordinates of the centroid

  void PrintCellList(); //print the cell list

  ~Cells();

};

