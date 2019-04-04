#ifndef UP_WIND_H
#define UP_WIND_H

#include <string>
#include <fstream>
#include <iostream>

#include "maillage.h"
#include "cellule.h"

class Up_wind
{
private:
    Maillage *m_M;
    Maillage *m_tmp;
    double m_Tfin;
    int m_nFonc;
    double a;


public:
    Up_wind(double Tfin,double inf_x, double sup_x,double in_y,double sup_y,int nb_cell,int n_fonc);   
      // void Raffinement();
      void solve();
      void solution();
      void saveMaillage();

private:
    void init_maill();
    double fonc(double x);
    //   double get_dt();
    // double get_Dx();
    // double get_Dy();
    //
    double delta_zmx(int i);
    double delta_zpx(int i);
    double phi(double a , double b);
    double delta_zmy(int i);
    double delta_zpy(int i);
    double sjx(int i);
    double sjy(int i);
    double z_plusx(int i);
    double z_moinsx(int i);
    double psix(int i,double u);
    double z_plusy(int i);
    double z_moinsy(int i);
    double psiy(int i, double u);
};

#endif // UP_WIND_H

