
exe : main.o cellule.o maillage.o up_wind.o
	g++ -o exe main.o cellule.o maillage.o up_wind.o 
cellule.o : cellule.cpp cellule.h
	g++ -c cellule.cpp -o cellule.o
up_wind.o : cellule.h maillage.h up_wind.h up_wind.cpp
	g++ -c up_wind.cpp -o up_wind.o
maillage.o :cellule.h maillage.h maillage.cpp
	g++ -c maillage.cpp -o maillage.o 
main.o : main.cpp cellule.h
	g++ -c main.cpp -o main.o
clean :
	rm *.o