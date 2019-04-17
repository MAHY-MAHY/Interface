#ifndef MAILLAGE_H
#define MAILLAGE_H

#include <vector>
#include <iostream>

#include "cellule.h"

class Maillage
{
public:
    Maillage();
    Maillage(double binf_x, double bsup_x,double binf_y,double bsup_y, int nb_cellx,int nb_celly,int CL);
    Maillage(const Maillage &m); //par recopie


    // MÉTHODES DE MODIFICATION
    // Ajout a une position
    void add_position(double valeur, int position);
    // Retrait a une position
    void remove_position(int position);
    // METHDOES DE CONSULTATION
    int consult_position(int position, int & valeur) const;
    // SURCHARGE D’OPERATEUR
    double operator[] (int position);

    int getNbCell() const;
    int getNbCellx() const;
    int getNbCelly() const;
    double getBinfCell_x(int i) const;
    double getBsupCell_x(int i) const;
    double getBinfCell_y(int i) const;
    double getBsupCell_y(int i) const;
    double getDxCell(int i) const;
    double getDyCell(int i) const;
    double getValueCell(int i) const;
    void init_maill(int m_nFonc);
    void setValueCell(int i, double val);
    double getValxy(int i,int j) const;

friend std::ostream & operator<<(std::ostream &os, const Maillage &M);

protected:
    int nb_cellules;
    int nb_cellx;
    int nb_celly;
    int CL;
    // Nombre de cellules du maillage.
    std::vector<Cellule*> m_V;
};


#endif // MAILLAGE_H
