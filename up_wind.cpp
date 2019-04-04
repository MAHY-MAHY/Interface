#include "up_wind.h"
#include <cmath>

Up_wind::Up_wind(double Tfin,double inf_x, double sup_x,double inf_y,double sup_y,int nb_cell, int n_fonc)
{
    m_Tfin = Tfin;
    m_nFonc = n_fonc;

    a =2.;// coefficient pour la fonction u

    m_tmp = new Maillage(inf_x,sup_x,inf_y,sup_y,nb_cell,nb_cell);
    m_M   = new Maillage(inf_x,sup_x,inf_y,sup_y,nb_cell,nb_cell);
    init_maill();
}

void Up_wind::solve()
{
  double t=0;
  int compt=0;
  double dt;
  double Dx;
  double Val;
  std::vector<double> f_tmp;
  f_tmp.reserve(m_M->getNbCell());
  f_tmp.resize(m_M->getNbCell());
  dt=m_M->getDxCell(0)*0.25*0.5;
  for(int i=0;i<m_M->getNbCell();i++){
    f_tmp.at(i)=0;
  }
    while(t<m_Tfin){
      for(int i=0;i<m_M->getNbCell();i++){
	f_tmp.at(i)=m_M->getValueCell(i)-dt*(psix(i,2))/m_M->getDxCell(i);
	f_tmp.at(i)=f_tmp.at(i)-dt*psiy(i,1)/m_M->getDyCell(i);
      } 
      for(int i=0;i<m_M->getNbCell();i++){
	Val=f_tmp.at(i);
	m_M->setValueCell(i,Val);
      }
      t=t+dt;
      if(compt%25==0){
	std::cout<<t<<std::endl;
      }
      compt=compt+1;
    }   
}

void Up_wind::solution()
{
    for(int i=0; i<m_M->getNbCell();i++){
      double mil_x,mil_y;
      mil_x=m_tmp->getBinfCell_x(i)+m_tmp->getDxCell(i)*0.5-2*m_Tfin;
      mil_y=m_tmp->getBinfCell_y(i)+m_tmp->getDyCell(i)*0.5-m_Tfin;
      if(mil_y<=0.5*mil_x){
	m_tmp->setValueCell(i,1);
      }
      else{
	m_tmp->setValueCell(i,0);
      }
    }
  
}
      

//this->Raffinement();
//std::swap(m_tmp, m_M);//Echange les deux maillages




//Initialisation du maillage et choix de la fonction initiale
void Up_wind::init_maill()
{
  if(m_nFonc==0){
    for(int i=0; i<m_M->getNbCell();i++){
      double mil_x,mil_y;
      mil_x=m_M->getBinfCell_x(i)+m_M->getDxCell(i)*0.5;
      mil_y=m_M->getBinfCell_y(i)+m_M->getDyCell(i)*0.5;
      if(mil_y<=0.5*mil_x){
	m_M->setValueCell(i,1);
      }
      else{
	m_M->setValueCell(i,0);
      }
    }
  }
  else{

    for(int i=0; i<m_M->getNbCell();i++){
      double mil;
      mil=m_M->getBinfCell_x(i)+m_M->getDxCell(i)*0.5;
      if(mil<-0.5 || mil>0.5){
	m_M->setValueCell(i,0);
      }
      else{
	m_M->setValueCell(i,1);
      }
    }
  }
}


//****Calcul des delta Z****//
double Up_wind:: delta_zmx(int i){
  double delta;
  if(((i)%m_M->getNbCellx())==0)
    {
      delta=0;
    }
  else{
    delta =m_M->getValueCell(i)-m_M->getValueCell(i-1);
      }
  return delta ;   
  }
double Up_wind:: delta_zpx(int i){
  double delta;
  if(((i+1)%m_M->getNbCellx())==0)
    {
      if(m_M->getBinfCell_y(i) + m_M->getDyCell(i)*0.5<\
	 0.5*(m_M->getBsupCell_x(i) + m_M->getDxCell(i)*1.5)\ 
	 && m_M->getBinfCell_y(i) + m_M->getDyCell(i)*0.5>\
	 0.5*( m_M->getBinfCell_x(i) + m_M->getDxCell(i)*0.5))
	{
	  delta=1;
	}
      else
	{
	  delta=0;
	}
    }
  else{
    delta =m_M->getValueCell(i+1)-m_M->getValueCell(i);
      }
  return delta ; 
}

double Up_wind :: delta_zmy(int i){
  double delta;
  if(i<m_M->getNbCellx())
    {
      delta=m_M->getValueCell(i)-1;
    }
  else{
    delta=m_M->getValueCell(i)-m_M->getValueCell(i-m_M->getNbCellx());
      }
  return delta ;   
}
double Up_wind :: delta_zpy(int i){
  double delta;
  if(i+m_M->getNbCellx()>m_M->getNbCell()-1)
    {
      delta=0;
    }
  else{
    delta =m_M->getValueCell(i+m_M->getNbCellx())-m_M->getValueCell(i);
  }
  return delta ; 
}

double Up_wind:: phi(double a, double b){
  double phi;
  double tmp1=std::abs(a);
  double tmp2=std::abs(b);
  /*if(a*b>0){
    phi=a/(tmp1)*(std::max(std::min(tmp1,2*tmp2),std::min(2*tmp1,tmp2)));
  }
  else{
  phi=0;*/
  if(a*b>0){
    phi=2*a/(tmp1)*std::min(tmp1,tmp2);
  }
  else{
  phi=0;
    
  }
  return phi;
}

double Up_wind :: sjx(int i){
  double s;
  double delta_m, delta_p;
  delta_m=delta_zmx(i);
  delta_p=delta_zpx(i);
  s=phi(delta_m,delta_p)/(m_M->getDxCell(i));
  return s;
  
}

double Up_wind :: sjy(int i){
  double s;
  double delta_m, delta_p;
  delta_m=delta_zmy(i);
  delta_p=delta_zpy(i);
  s=phi(delta_m,delta_p)/(m_M->getDyCell(i));
  return s;  
}

double Up_wind :: z_plusx(int i){
  double z;
  z=m_M->getValueCell(i+1)-0.5*m_M->getDxCell(i+1)*sjx(i+1);
  return z;    
}

double Up_wind :: z_moinsx(int i){
  double z;
  z=m_M->getValueCell(i)+0.5*m_M->getDxCell(i)*sjx(i); ///// A verif le  i+1 ?!!!
  return z;    
}

double Up_wind ::psix(int i, double u){
  double psi;
  double psi_m=0;
  double psi_p=0;
  if(i%m_M->getNbCellx()==0){
    psi_m=0;//u*m_M->getValueCell(i);
    psi_p=u*z_moinsx(i);
  }
  else if((i+1)%m_M->getNbCellx()==0){
    if(m_M->getBinfCell_y(i) + m_M->getDyCell(i)*0.5<\
       0.5*(m_M->getBsupCell_x(i) + m_M->getDxCell(i)*1.5)\ 
       && m_M->getBinfCell_y(i) + m_M->getDyCell(i)*0.5>\
       0.5*( m_M->getBinfCell_x(i) + m_M->getDxCell(i)*0.5))
      {
	psi_p=u;
      }
    else
      {
	psi_p=u*m_M->getValueCell(i);
      }
      psi_m=u*z_moinsx(i);
  }  
  else{
    psi_p=u*z_moinsx(i);
    psi_m=u*z_moinsx(i-1);
  }
  psi=psi_p-psi_m;
  return psi;
} 

double Up_wind :: z_plusy(int i){
  double z;
  z=m_M->getValueCell(i+1)-0.5*m_M->getDyCell(i+1)*sjy(i+1);
  return z;    
}

double Up_wind :: z_moinsy(int i){
  double z;
  z=m_M->getValueCell(i)+0.5*m_M->getDyCell(i)*sjy(i);
  return z;    
}

double Up_wind ::psiy(int i,double u){
  double psi;
  double psi_m=0;
  double psi_p=0;
  if(i<m_M->getNbCellx()){
      if(m_M->getBinfCell_y(i) - m_M->getDyCell(i)*0.5<\
	 0.5*(m_M->getBinfCell_x(i) + m_M->getDxCell(i)*0.5)\ 
	 && m_M->getBinfCell_y(i) + m_M->getDyCell(i)*0.5>\
	 0.5*( m_M->getBinfCell_x(i) + m_M->getDxCell(i)*0.5)){
	psi_m=u;//*m_M->getValueCell(i);
      }
      else{
	psi_m=u*m_M->getValueCell(i);
      }
      psi_p=u*z_moinsy(i);
  }
  else if(i+m_M->getNbCellx()>m_M->getNbCell()-1){
    psi_p=u*m_M->getValueCell(i);
    psi_m=u*z_moinsy(i);
  }  
  else{
    psi_p=u*z_moinsy(i);
    psi_m=u*z_moinsy(i-m_M->getNbCellx());
  }
  psi=psi_p-psi_m;
  return psi;
} 

void Up_wind::saveMaillage()
{
    std::ofstream myfile;
    std::ofstream myfile2;
    std::ofstream myfile3;
    myfile.open("outY.txt",std::ios::trunc);
    myfile2.open("outX.txt",std::ios::trunc);
    myfile3.open("outYex.txt",std::ios::trunc);
 

    for(int i=0; i<m_M->getNbCell(); ++i)
      {
	myfile <<m_M->getValueCell(i) << "\n";
	myfile2 << m_M->getBinfCell_x(i) + m_M->getDxCell(i)*0.5 << "   " \
		<< m_M->getBinfCell_y(i) + m_M->getDyCell(i)*0.5 << "\n" ;
	myfile3 <<m_tmp->getValueCell(i) << "\n";
	if(i<m_M->getNbCell() - 1)
	  {
	    myfile3 << " ";
	    myfile2 << " ";
	    myfile << " ";
	  }
      }
    myfile << std::endl;
    myfile2 << std::endl;
    myfile3 << std :: endl;
    myfile.close();
    myfile2.close();
    myfile3.close();
}



//////////********////////

/* void Up_wind:: Raffinement(){
    for(int i=0; i<m_M->getNbCell()-1 ;i++){
        if(std::abs(m_M->getValueCell(i)-m_M->getValueCell(i+1))>m_Eps && m_M->getDxCell(i)>m_dx_crit){// on raffine le maillage sous certaines conditions
            double Val_i=m_M->getValueCell(i);
            m_M->add_position(Val_i,i);
            m_tmp->add_position(Val_i,i);

            i++;
        }
    }

}
*/
/*
double Up_wind::get_dt()
{
    double dx_min=m_M->getDxCell(0);
    double u_max=m_M->getValueCell(0);
    for (int i=1;i<m_M->getNbCell(); i++){
        if(u_max<m_M->getValueCell(i)){
            if(m_nFonc==0){
                u_max=m_M->getBsupCell_x(i);
		  }
            else{
                u_max=a;
            }
        }
        if(dx_min>m_M->getDxCell(i)){
            dx_min=m_M->getDxCell(i);
        }
    }
    return dx_min/u_max;

}
double Up_wind::get_Dx(){
    double dx_min=m_M->getDxCell(0);
    double u_max=m_M->getValueCell(0);
    for (int i=1;i<m_M->getNbCell(); i++){
            if(dx_min>m_M->getDxCell(i)){
            dx_min=m_M->getDxCell(i);
        }
    }
    return dx_min;

}

double Up_wind::get_Dy(){
    double dy_min=m_M->getDyCell(0);
    double u_max=m_M->getValueCell(0);
    for (int i=1;i<m_M->getNbCell(); i++){
            if(dy_min>m_M->getDyCell(i)){
            dy_min=m_M->getDyCell(i);
        }
    }
    return dy_min;
}

*/
