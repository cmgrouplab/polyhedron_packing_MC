//including all the constants shared by several rountines to complete the MC packing algorithm ...

#define SD 3 //the spatial dimensions...
#define MAXX 32000 //a random number mod

#define pi 3.14159265358979 //pi value....


//this is for the structure of cells.
struct node{
  int index_poly; //the index of. in a certain cell...
  node* next;
};


