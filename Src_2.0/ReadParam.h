//a class handles all the parameters for the MC NPT simulation



class ReadParam{
 
 private:
  
  int Nstage; //on each the box is fixed, at the end of a stage deformation is performed
  int Ncycle; //number of cycles per stage
  int Ndeform; //number of deformation performed at the end of each stage


  //specify the magnitude of the displacements...
  double TransMod; //associated with the relative coordinates...no unity needs to be specified...
  double RotMod;
  double Global_StrainMod; // the default strain mod: >0 is expansion and <0 is compression
  double Shear_StrainMod; //separate shear strain and principle strain

  double Strain_Rescale; //if compression fails, rescale the strain mod by this amount each time
  double Trans_Rescale; //for adpative step sizes, make sure acceptance rate is around 0.5
  double Rot_Rescale;
  //when rescaling the step size, also rescaling the strain...
  //so the specified step size and strain should be comparable, otherwise the real strain would small than expected
  
  double p_uphill; //the generation rate of uphill moves...
  double p_trans; //probabiliy of the current move is a translation 
  
  double relax_ratio; //if the initial configuration contains overlapping, rescale the lattice
  
  double Starting_Density;
  double Terminate_Density;
  
  int Npc; //frequence of statiscs collection, i.e., the number of moves for each particle between two collections
  int flag_comput_pressure;
  int flag_comput_g2;
  int flag_comput_ordermetric;
  int check_option; //indicate how to deal with int config containing overlapping particles...


 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
 public:

  double nnl_cut_dist; //cut off dist for nnl
  double rho_nnl_th; //density threshold for using nnl
  
  //char* poly_name; //the name of the polyhedron...

  ReadParam();
  
  void Rescale_TransMod();
  void Rescale_RotMod();
  void Rescale_Strain();

  int Get_Nstage();
  int Get_Ncycle(); //number of cycles per stage
  int Get_Ndeform(); //number of deformation performed at the end of each stage


  //specify the magnitude of the displacements...
  double Get_TransMod(); //associated with the relative coordinates...no unity needs to be specified...
  double Get_RotMod();
  double Get_Global_StrainMod(); // the default strain mod: >0 is expansion and <0 is compression
  double Get_Shear_StrainMod(); //separate shear strain and principle strain

  double Get_p_uphill(); //the generation rate of uphill moves...
  double Get_p_trans(); //probabiliy of the current move is a translation 
  
  double Get_relax_ratio(); //if the initial configuration contains overlapping, rescale the lattice
  
  double Get_Starting_Density();
  double Get_Terminate_Density();

  int Get_Npc();
  int flag_pressure();
  int flag_g2();
  int flag_ordermetric();
  int Get_check_option();

  double Get_Strain_Rescale();
};




