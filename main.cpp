#include <iostream>

#include "cellule.h"
#include "maillage.h"
#include "up_wind.h"

int main()
{
  //Maillage test(0,1,0,1,10);
  Up_wind test2(2,0,2,0,2,200,0);//les paramètres sont(Tfin, Xdebut,Xfin, nombre de point initial du mailage,taile minimal d'une cellule, choix de la fonction initial,saut maximum ou raffiner)
  //std::cout << test <<std::endl;
   test2.solve();
   test2.solution();
   test2.saveMaillage();
   
  return 0;

}
