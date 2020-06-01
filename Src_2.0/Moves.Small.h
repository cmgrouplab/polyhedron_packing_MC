//Perform translation and rotation on the particle --- 
// ---- generate a new position or a new particle config, 
//       if the move is accepted, copy the new ones to the old ones
//       if the move is rejected, no operation is needed -- this will speed up the program
//Perform lattice deformation
// ---- since the lattice moves are rare compared to particle moves 
//      we keep the old implementation, otherwise need to pass too many arguments...



class Moves{

 public:
  //the temp position and config.
  double Strain[SD][SD]; //the strain tensor
  double TLambda[SD][SD];
  double TCenterL[SD];
  Polyhedron TParticle; //need to fill in the shape later...

  Geometry Comput;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Moves(char*); //provide the particle shape

  void GetTranslation(double, double* &); //param: (1) TransMod (2) CenterLm
  void RetainTranslation(double* &); // param: (1) CenerLm

  void GetRotation(double, Polyhedron &); //param:(1) RotMod (2) Particle config
  void RetainRotation(Polyhedron &);

  int ParticleMove(double, double, double, int, Packing &);
  //the parameters (1) the probability that the move is a translation
  //               (2) TransMod; (3) RotMod
  //               (4) Index of the moved particle
  //               (5) The Packing configuration including both position and particle config
  //               (6) THe cells config. for the packing

  // int ParticleMoveII(double, double, double, int, Packing &);
   //this uses NNL and SEP_List


  void GetStrainIso(double, double, double); 
  // parameters: (1) p_uphill (2) P_StrainMod (3) S_StrainMod
  void ResetStrain(double, double, int, double, double);
  // parameters: (1) p_uphill (2) Strain_Rescale ratio (3) a rescale power (4) Global_StrainMod (5) Shear_StrainMod
  
  void BoundaryDeform(double** &); //the parameter is the lattice vector of the packing
  void BoundaryRetain(double** &);

  int BoundaryMove(int, double, double, double, double, Packing &);
  //the parameters (1) Maximum number of deformation trials before declare a fail 
  //               (2) p_uphill (3) Strain_Rescale (4) Global_StrainMod (5) Shear_StrainMod
  //               (6) The packing object
  //               (7) The cell object for the packing


};


