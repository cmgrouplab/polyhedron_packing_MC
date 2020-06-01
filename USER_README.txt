
Adaptive Shrinking Cell MC Hard-Polyhedron Packing Program:

Author: Yang JIAO, yjiao@princeton.edu

Note: SRC_2.0 version: refined the overlap-check method, but the 
improvement in efficiency not much; tried near-neighbor-list and 
separation list: these almost not helpful at all unless for optimizing 
a perturbed crystal...



This is brief explaination of the C++ codes I developped for simulation hard 
polyhedron systems. The program can deal with ANY Polyhedron shapes including 
regular or irregular, monodisperse or polydisperse, as long as all shapes are 
of the same type (e.g., they are all tetrahedra, but each tetrahedron can be 
different from one another). With slight modifications, the codes can simulate 
systems composed of mixed shapes (e.g., tetrahedra and octahedra, etc.).

Since the codes are designed to simulate equilibrium behavior (quasi NPT ensemble)
and produce disordered jammed packings, cell method is empolyed to improve 
efficiency. However, the cell method only work properly when the number of 
particles is greater than a particular value, i.e., N > 30. Small N is 
crucial for the search of the densest packings of monodisperse shapes, and 
slight modification is needed for the current codes working for small systems. 
Note that small-system codes are developped for each of the Platonic solids separately.

In addition, the codes are designed to work for any general polyhedron. 
There is a balance between generality and efficiency. The key step for the simulation ---
overlap check -- is based on the Separation Axis Theorem (SAT), which requires operations 
on all faces and pair of edges for a pair of particles in the system. For particular 
shapes (e.g., tetrahedron), efficient implementation of SAT has been developped. However, 
such efficiency cannot be achieved for a general polyhedral shape. Therefore, 
currently, we are "doing it correctly, but may not efficiently"...

The following are all the related source code files:

GEsolver: solves linear algebric equations for transformation between relative and global coordinates

Geometry: deal with geometry, i.e., inner produce, cross product, metric, vector length, etc.

Polyhedron: deal with the shape, i.e., vertices, faces, edges, volumes etc.

Packing: for the packing, the central portion for the program; containing lattice vectors, centroid 
         positions, orientation information for each particle, overlap checking functions, density 
         calculator, packing print functions, etc.

Cells: to divide the simulation domain into small cells

Moves: for the MC moves, includes translation, rotation, box deformation/shrinkage/expansion and 
       restoration of the original configuration.

Stat:
	OrderMetric: for orientational correlations, (1) allignment (2) face-normal

	PairCorr: for pair correlation functions (1) along each direction (2) radiually averaged

Pressure: for the pressure of the system

ReadParam: read the dynamical parameter <input.txt> for the simulation

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To compile and produce executables, simply use [make], i.e.,

./make

To run the MC simulation, use the commend:

./MC_NPT.x $poly_name < input.txt

where poly_name = {P1, P2, P3, P4, P5, U} for the Platonic solids and user-specified type.


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#brief instructions for using MC_NPT.x for user-specified shapes:

input files:
	shape_param.txt: the class "Polyhedron" is initialized from this file. It contains alll
		face, edge, vertex configuration information. The actually coordinates of the 
		vertices will be read in from "Iconfig.txt"
		>>>>>> IMPORTANT: the order of vertices in Iconfig.txt must be consistent with that in "shape_param.txt"
	Iconfig.txt: the "classic" initial configuration....
		In the first round, d0 should be the edge length, to make sure it can be also properly 
			processed by the old version of the codes. In the "config.txt", d0 is the circum-diameter
			of the shape...
	input.txt: specifies all the parameters needed by the program...

To run: ../../MC_NPT.x $Poly_Name <input.txt

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#brief instructions for using init.x:

input files: 
	MRJ.txt: configuration of jammed hard spheres, providing centroid information 
		and number of particles...
	Poly_Shape.txt: provide the name of the shape, the number of vertices and 
		 vertex coordinate information...

To run: ./init.x $Poly_Name
	$Poly_Name:= P1, P2, P3, P4, P5, U
	The order of the vertices in Poly_Shape.txt is very important: need to be consistent 
		with the built-in order in class "Polyhedron"
	for the user-specified shape, the objects of class "Polyhedron" are completely
		initialized from input file "shape_param.txt", the vertice in "Poly_Shape.txt"
		need to be consistent with that in "shape_param.txt" 
	The program will ask for a initial density, if the best density it can get is smaller 
		than the specific rho_int, it will ignore that value; otherwise, it will rescale 
		the box to produce the target density...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The input file for dynamical parameters: input.txt

 The number of stages -  Nstage := 
 The number of cycles per stage - Ncycle := 
 The number of deformation per stage - Ndeform := 
 The translation magnitude - TransMod := 
 The rotation magnitude (in unit pi) - RotMod := 
 Normal strain magnitude - Global_StrainMod := 
 Shear strain magnitude - Shear_StrainMod := 
 Strain rescale magnitude - Strain_Rescale := 
 Translation rescale magnitude - Trans_Rescale := 
 Rotation rescale magnitude - Rot_Rescale := 
 Probability of uphill moves - p_uphill := 
 Probability of translation for a trail move - p_trans := 
 Rescale paramter for initial configuration - relax_ratio :=
 Starting density - Starting_Density := 
 Termination density - Terminate_Density := 
 Rescale the packing or Exit if Iconfig.txt containing overlapping particles? Rescale - 0; Exit - 1
 The number of trial moves between each property collection - Npc := 
 Collecting Pressure? 1 - Yes; 0 - No.
 Computing g2? 1 - Yes; 0 - No.
 Computing orientational order metric? 1 - Yes; 0 - No.


$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
The input file for shape parameters: shape_param.txt

poly_name  //name of the shape
n_vertex   //number of vertices
n_edge     //number of edges
n_face     //number of faces
n_fs       //number of different faces species, classified by the number of vertices
max_face_vert  //maximum number of vertices for a face among all face species

//number of faces for each species...
face_num_fspecie[0]
....
face_num_fspecie[n_fs]


//number of vertices for the face in each face species....
vert_num_fspecie[0]
....
vert_num_fspecie[n_fs]


//the two vertices for an edge
Edge[0][0]       Edge[0][1];
....
Edge[n_edge][0]  Edge[n_edge][1];


//the specie of the current face + face vertex configuration
face_type[0]         Face[0][0]      ... Face[0][n_0]
....
face_type[n_face]    Face[n_face][0] ... Face[n_face][n_0]




$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
The input file for initial configuration: Iconfig.txt

% N # the number of particles in the fundamental cell
% d_E  # characteristic length of the particles
% a1_x  a1_y  a1_z  # the three components of lattice vector a_1
% a2_x  a2_y  a2_z  # the three components of lattice vector a_2
% a3_x  a3_y  a3_z  # the three components of lattice vector a_3

% O1_x  O1_y  O1_z  # the relative coordinates of the particle centroid 
% A1_x  A1_y  A1_z  # the coordinates of vertex A of tetrahedron 01
% B1_x  B1_y  B1_z  # the coordinates of vertex B of tetrahedron 01
% C1_x  C1_y  C1_z  # the coordinates of vertex C of tetrahedron 01
% D1_x  D1_y  D1_z  # the coordinates of vertex D of tetrahedron 01

.......

% ON_x  ON_y  ON_z  # the relative coordinates of the particle centroid 
% AN_x  AN_y  AN_z  # the coordinates of vertex A of tetrahedron 72
% BN_x  BN_y  BN_z  # the coordinates of vertex B of tetrahedron 72
% CN_x  CN_y  CN_z  # the coordinates of vertex C of tetrahedron 72
% DN_x  DN_y  DN_z  # the coordinates of vertex D of tetrahedron 72

For other shapes, the format is similar, just need more vertices for each face



$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
The input file for the shape vertex: Poly_Shape.txt

poly_name //name of the polyhedron 
n_vert // number of vertices
x_1      y_1      z_1  // coordinates of each vert...
.....
x_n_vert y_n_vert z_n_vert

This is for generating initial configuration <Iconfig.txt>, i.e., 
input file for [Init_Config.x], which read in the shape information 
from <shape_param.txt> for general shapes (or simply use the built-in 
shape info for the Platonic solids) and read in the vertex of 
prototype shape from <shape_vert.txt>. 

Then [Init.x] will produce an <Iconfig.txt> file with 
user-specified density in a cubic simulation box, either using RSA <RSA.txt> process 
or using a MRJ sphere packing <MRJ.txt> (by inserting a polyhedron with 
arbitary orientation into the sphere and rescale the box properly). 



$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
For ProcessPacking.x:

./ProcessPacking.x P1

which reads in Iconfig.txt, which is the Fconfig.txt of P1.
