#include <iostream>
#include <cstdio>
#include <ctime>
#include <omp.h>

#include "cellule.h"
#include "maillage.h"
#include "up_wind.h"
#include "Diffuse.h"

int main()
{
  //Maillage test(0,1,0,1,10);
  double a_x,a_y , b_x, b_y,T_f,montre1,montre2,CFL,Beta;
  int nb_cel, init ,methode,CL;
  char ok;
  std::cout << "lire le temps final"<<std::endl;
  std::cin >> T_f;
  std::cout<<" lire borne inf de x et borne sup de x"<<std::endl;
  std::cin>>a_x >>b_x;
  std::cout<<" lire borne inf de y et borne sup de y"<<std::endl;
  std::cin>>a_y >>b_y;
  std::cout<<" Entrer 0 pour droite 1 pour cercle en translation 2 pour cercle en rotation 3 pour KOTHE "<<std::endl;
  std::cin>>init;
  std::cout<<" Entrer 0 pour sharp 1 pour diffuse "<<std::endl;
  std::cin>>methode;
  std::cout<<" le nombre de cellule en X "<<std::endl;
  std::cin>>nb_cel;
  std::cout<<" Entrer le nombre CFL  "<<std::endl;
  std::cin>>CFL;
  std::cout<<"Entrer 0 our condition periodique 1 pour Neumann" <<std::endl;
  std::cin>>CL;
  std::cout<< "Entrer Beta superieur à 0"<<std::endl;
  std::cin>>Beta;
  if(methode==0){
    Up_wind test2(T_f,a_x,b_x,a_y,b_y,nb_cel,init,CFL,CL);
    montre1 = omp_get_wtime();
    test2.solve_sharp();
    montre2 =omp_get_wtime();
    std::cout<<montre2-montre1<<std::endl;
    test2.solution();
    test2.saveMaillage();
  }
  else{
    Diffuse test2(T_f,a_x,b_x,a_y,b_y,nb_cel,init,CFL,CL,Beta);
    montre1 = omp_get_wtime();
    test2.solve_diffuse();
    montre2 =omp_get_wtime();
    std::cout<<montre2-montre1<<std::endl;
    test2.solution();
    ok =a_x;
    test2.saveMaillage();
    std::cout<<"ok"<<std::endl;
  }
   
  return 0;

}
