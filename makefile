
export OMP_NUM-THREADS=8
OPT= -fopenmp -O2
exe : main.o cellule.o maillage.o up_wind.o diffuse.o
	g++ $(OPT) -o exe main.o cellule.o maillage.o up_wind.o diffuse.o
cellule.o : cellule.cpp cellule.h
	g++ $(OPT)  -c cellule.cpp -o cellule.o
up_wind.o : cellule.h maillage.h up_wind.h up_wind.cpp
	g++ $(OPT) -c up_wind.cpp -o up_wind.o
diffuse.o : cellule.h maillage.h Diffuse.h diffuse.cpp
	g++ $(OPT) -c diffuse.cpp -o diffuse.o
maillage.o :cellule.h maillage.h maillage.cpp
	g++ $(OPT)  -c maillage.cpp -o maillage.o 
main.o : main.cpp cellule.h
	g++ $(OPT) -c main.cpp -o main.o
clean :
	rm *.o