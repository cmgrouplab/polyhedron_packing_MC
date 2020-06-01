#makefile for polyhedron codes

MC_NPT: MC_NPT_Polyhedron.o Geometry.o Polyhedron.o Cells.o Packing.o Moves.o ReadParam.o
	g++ -o MC_NPT.x MC_NPT_Polyhedron.o Geometry.o Polyhedron.o Cells.o Packing.o Moves.o ReadParam.o

test: Test.o Geometry.o Polyhedron.o Cells.o Packing.o Moves.o ReadParam.o
	g++ -o test.x Test.o Geometry.o Polyhedron.o Cells.o Packing.o Moves.o ReadParam.o

GEsolver.o: GEsolver.h GEsolver.cpp  const.h
	g++ -c GEsolver.cpp

Geometry.o:  Geometry.h Geometry.cpp  const.h
	g++ -c Geometry.cpp

Polyhedron.o:  Polyhedron.h Polyhedron.cpp Geometry.h const.h
	g++ -c Polyhedron.cpp

Cells.o: Cells.h Cells.cpp const.h
	g++ -c Cells.cpp

Packing.o: Packing.h Packing.cpp Polyhedron.h Geometry.h Cells.h const.h
	g++ -c Packing.cpp

Moves.o: Moves.h Moves.cpp Geometry.h Polyhedron.h Cells.h Packing.h const.h 
	g++ -c Moves.cpp

ReadParam.o: ReadParam.cpp ReadParam.h
	g++ -c ReadParam.cpp

clean:
	rm *.o

#EOF