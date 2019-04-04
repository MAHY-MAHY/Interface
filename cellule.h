#include <iostream>
#ifndef cellule_h
#define cellule_h

class Cellule{
  friend class Maillage;
public:
  Cellule();
  Cellule(const Cellule &c);
  Cellule(double inf_x, double sup_x,double inf_y,double sup_y ,double valeur);

  double getValue() const;
  void setValue(double v);

  double getBinf_x() const;
  void setBinf_x(double inf);

  double getBsup_x() const;
  void setBsup_x(double sup);

  double getBinf_y() const;
  void setBinf_y(double inf);

  double getBsup_y() const;
  void setBsup_y(double sup);


private:
  double a_x,b_x;
  double a_y,b_y;
  double uin;


};

std::ostream &operator<<(std::ostream &os, const Cellule &C);

#endif
