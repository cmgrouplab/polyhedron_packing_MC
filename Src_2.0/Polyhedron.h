//A class that defines a generic polyhedron shape
//necessary functions are defined here


class Polyhedron{

 public:

  char* poly_name; //the name of the polyhedron
  int n_vertex; //number of vertices
  int n_edge; //num of edges
  int n_face; //num. of faces
  int n_fs; //num of face species
  int max_face_vert; //the maxium number of vertices for a face if more than one face-type are present
  
  double** Vertex; //vertex coordinates
  int** Edge; //edge configuration
  int** Face; //face configuration
  int* face_type; //which face species Face[i] belongs to, for user specified general shapes
  int* face_num_fspecie; //how many faces for each specie
  int* vert_num_fspecie; //how many vertice for the face of each specie
  
  double* d_cf; //centroid to face distance
  double* d_cv; //centroid to vertex distance

  double** Normal; //outward face normal
  
  Geometry Comput; //for computing the products
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Polyhedron(); //for "new []", which can not call constructor with parameters
  Polyhedron(char*); //initialize polyhedron based on given name
  //NOTE: the order of the vertex is crucial!!!!
  void PolyConstr(char*); //for "new []" which cannot use constructors with parameters
  //exactly the same as Polyhedron(char*)
  void ShiftVert(); //shift the vert to make the centroid at 0;
  void GetNormal_All(); //get the normal for for all faces of polyedron
  void GetNormal_F(int); //only get the normal for face with index int face_index
  //this is a more efficient way for overlap check, rotational move etc...
  void GetCentDist(); //compute d_cf and d_cv
  void CopyPoly(Polyhedron &); //copy the later to the former
  Polyhedron& operator = (Polyhedron &); //reload operator
  double GetVol(); //comput the vol of the particle
  ~Polyhedron(); //desturcturor
};


