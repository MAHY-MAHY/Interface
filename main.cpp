#include <iostream>
#include <cstdio>
#include <ctime>
#include <omp.h>

#include "cellule.h"
#include "maillage.h"
#include "up_wind.h"

int main()
{
  //Maillage test(0,1,0,1,10);
  double a_x,a_y , b_x, b_y,T_f,montre1,montre2,CFL;
  int nb_cel, init ,methode;
  std::cout << "lire le temps final"<<std::endl;
  std::cin >> T_f;
  std::cout<<" lire borne inf de x et borne sup de x"<<std::endl;
  std::cin>>a_x >>b_x;
  std::cout<<" lire borne inf de y et borne sup de y"<<std::endl;
  std::cin>>a_y >>b_y;
  std::cout<<" Entrer 0 pour droite 1 pour cercle en rotation 2 pour cercle en translation "<<std::endl;
  std::cin>>init;
  std::cout<<" Entrer 0 pour sharp 1 pour diffuse "<<std::endl;
  std::cin>>methode;
  std::cout<<" le nombre de cellule en X "<<std::endl;
  std::cin>>nb_cel;
  std::cout<<" Entrer le nombre CFL  "<<std::endl;
  std::cin>>CFL;

  Up_wind test2(T_f,a_x,b_x,a_y,b_y,nb_cel,init,CFL);//les paramètres sont(Tfin, Xdebut,Xfin, nombre de point initial du mailage,taile minimal d'une cellule, choix de la fonction initial,saut maximum ou raffiner)
  //std::cout << test <<std::endl;
  montre1 = omp_get_wtime();

  std::cout<< "ok"<< std::endl;
  test2.solve_sharp();
  montre2 =omp_get_wtime();
  std::cout<<montre2-montre1<<std::endl;
   test2.solution();
   test2.saveMaillage();
   
  return 0;

}
