#makefile for polyhedron codes

MC_NPT: MC_NPT_Polyhedron.o Geometry.o Polyhedron.o Cells.o Packing.o Moves.o ReadParam.o Stat.o
	g++ -o MC_NPT.x MC_NPT_Polyhedron.o Geometry.o Polyhedron.o Cells.o Packing.o Moves.o ReadParam.o Stat.o

MC_NPT.Small: MC_NPT_Polyhedron.Small.o Geometry.o Polyhedron.o Packing.o Moves.Small.o ReadParam.o Stat.o
	g++ -o MC_NPT.Small.x MC_NPT_Polyhedron.Small.o Geometry.o Polyhedron.o Packing.o Moves.Small.o ReadParam.o Stat.o

init: InitConfig.o Polyhedron.o Geometry.o GEsolver.o
	g++ -o init.x InitConfig.o Polyhedron.o Geometry.o GEsolver.o

test: Test.o Geometry.o Polyhedron.o Cells.o Packing.o Moves.o ReadParam.o Contact.o
	g++ -o test.x Test.o Geometry.o Polyhedron.o Cells.o Packing.o Moves.o ReadParam.o Contact.o

ProcessPacking: ProcessPacking.o  Geometry.o Polyhedron.o Packing.o Stat.o Contact.o
	g++ -o ProcessPacking.x ProcessPacking.o Geometry.o Polyhedron.o Packing.o Stat.o Contact.o

view: ViewPacking.o  Geometry.o Polyhedron.o Packing.o Stat.o 
	g++ -o view.x ViewPacking.o Geometry.o Polyhedron.o Packing.o Stat.o 

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

Moves.Small.o: Moves.h Moves.cpp Geometry.h Polyhedron.h Packing.h const.h 
	g++ -c Moves.Small.cpp

ReadParam.o: ReadParam.cpp ReadParam.h
	g++ -c ReadParam.cpp

Stat.o: Stat.h Stat.cpp Packing.h Polyhedron.h Cells.h const.h Geometry.h
	g++ -c Stat.cpp

Contact.o: Contact.h Contact.cpp Packing.h Polyhedron.h const.h Geometry.h
	g++ -c Contact.cpp

clean:
	rm *.o

#EOF