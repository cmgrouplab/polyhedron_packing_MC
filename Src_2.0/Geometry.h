//a class the contains all functions dealing with geometriy 
// and other misc functions such as getmax and get min



class Geometry{

 public:

  void GetCrossProduct(double [], double [], double []);
  double GetInnerProduct(double [], double []);
  void GetPLine(double [], double [], double []);
  //get the common perpendicular line for two lines
  double GetLength(double []); //get the length of a vector
  void Normalize(double []); //normalize a vector
  double Normalize(double [], int); //normalize a vector, return its length, the latter parameter is just a flag
  
};

