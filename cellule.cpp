#include "cellule.h"

//Constructeurs de cellule
Cellule::Cellule():a_x(-1.),b_x(1.),a_y(-1),b_y(1), uin(0.){
}

Cellule::  Cellule(double inf_x, double sup_x,double inf_y,double sup_y ,double valeur){
    a_x=inf_x;
    b_x=sup_x;
    a_y=inf_y;
    b_y=sup_y;
    uin = valeur;
}

Cellule::Cellule(const Cellule &c){
    a_x = c.a_x;
    b_x = c.b_x;
    a_y = c.a_y;
    b_y = c.b_y;
    uin = c.uin;
}

//retourne la valeur de la cellule
double Cellule::getValue() const
{
    return uin;
}

//modifie la valeur de la cellule
void Cellule::setValue(double v)
{
    uin = v;
}

//retourne la borne inf de la cellule
double Cellule::getBinf_x() const
{
    return a_x;
}

//modifie la borne inf de la cellule
void Cellule::setBinf_x(double inf)
{
    a_x = inf;
}

//Retourne la borne sup de la cellule

double Cellule::getBsup_x() const
{
    return b_x;
}

//modifie la borne sup de la cellule
void Cellule::setBsup_x(double sup)
{
    b_x=sup;
}

//retourne la borne inf de la cellule
double Cellule::getBinf_y() const
{
    return a_y;
}

//modifie la borne inf de la cellule
void Cellule::setBinf_y(double inf)
{
    a_y = inf;
}

//Retourne la borne sup de la cellule

double Cellule::getBsup_y() const
{
    return b_y;
}

//modifie la borne sup de la cellule
void Cellule::setBsup_y(double sup)
{
    b_y=sup;
}

//Surcharge del'operateur d'affichage
std::ostream &operator<<(std::ostream &os, const Cellule &C)
{
  os << "[" << C.getBinf_x() << "," << C.getBinf_y() <<"] \n" << C.getValue() << "\n [" << C.getBsup_x() << "," << C.getBsup_y() << "]";
    return os;
}
