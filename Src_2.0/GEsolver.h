//Solves a normal linear algebra system using Guassian Elimination

// global variable declarations
class GESOL{
 public:
int nun; // number of unknowns
int neq; // number of equations
double **sys; // linear algebra system
double *back; // solution vector
int moreeq; // whether the neq>nun

// function prototypes
 GESOL();
void solution();
double evaluate(int);
void result();
int consistent();
void GEsolver(int, int, double [][SD], double []);
~GESOL();
};


