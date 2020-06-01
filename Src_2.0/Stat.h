//Header file for obtaining statistics of the packing, including:
//     (1) pair correlation functions
//     (2) allignment correlation functions 
//     (3) face-normal correlation functions



//#include "Packing.h" 
//pass the & of a packing to the functions of Stat
//this can be taken care of by listing the header files 

class Stat{

 public:
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //the pair correlation function...
  
  int N_bin; //the number of bins per d0
  double bin; //the bin width...
  double MaxL; //the characteristic length of the simulation box, i.e., the length of the shortest lattice vector
  int NLcounter; //the actual sample length
  int MaxNL; // for the sampling length, the upper limit
  
  //the pair distribution function...
  double*** G; //[MaxNL]
  double* g; //[MaxNL]
  int* g_Num; //[MaxNL], store the number in each bin, for computation of other order metric correlation
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //the local oreintational order parameter... 
  
  double** NormalOP; //[N][N]
  //defined to be fabs(min{<n_i , n_j>}-1)/2; the minimum of the inner product of the face normal
  //the reference system: face-to-face glued di-tetrah -> perfect order, since this local arrangement corresponds to higher density
  //                      perfectly alligned two tetrah -> most disorder, since this corresponds to vertex-contact, leading to low density packing
  
  double* LocalOrder; //[N], NormalOP for contacting neighbours, average to given the value for the central one
  double Ave_LocalOrder;
  double* g_Normal; //[MaxNL] the correlation functions
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //the allignment parameter...
  
  double** AllignOP; //[N][N]
  //defined to be max(<n_i, n_j>+<v(n_j)_i, v(n_j)_j>)/2; the maximum face normal allignement + maximum center2vertex allignment
  
  double* LocalAllign; //[N], averge over all contacting numbers
  double Ave_LocalAllign; //average over all particles...
  double* g_Allign; //[MaxNL] the correlation functions
  

  //*********************************************************************



  Stat();
  Stat(int, Packing &); //the constructor, 
                        //get the sampling parameters either by default (1) or user specified (0)
                        //need the characteristics length of the particle...

  double MinDis(int, int, Packing &, int, int, int); //the minimal distance between two particles, 
                                                     //record the box indices
 
  double Get_NDensity(Packing &); //get the number density of the packing to normalize g2

  void Get_PairCorr(Packing &); //compute both the vector and radial-averaged g2

  void Get_NormalOP(Packing &); //NormalOP is simply the minus of the smallest face normal inner products of any two paritcles
                     //Here, we want to characterize the degree of f2f contact and do not distinguish faces with different types
  
  void Get_LocalNormOrder(Packing &); //compute the face-normal correlation distribution
                    //get the local order face-normal parameter for each particle, contributed from nearest neighbors
                    //then get the average over all particles

  void Get_AllignOP(Packing &); //AllignOP is already the order parameter measuring the allignment of particle
                  //the current way allignOp is computed has given normal allignement more weight
                  //when any two faces of the same kind perfectly allign with each other, the two particles also allign (rigorous)
                  //we simply consider the faces  


  void Get_LocalAllign(Packing &); //compute the allignment correlation distribution
                  //get the local allignment order parameter for each particle, contributed from nearest neighbors
                  //then get the average over all particles

};



