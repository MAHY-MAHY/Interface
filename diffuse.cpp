#include "Diffuse.h"
#include <omp.h>


Diffuse:: Diffuse(double Tfin,double inf_x, double sup_x,double inf_y,double sup_y,int nb_cell, int n_fonc,double CFL,int Cond_lim,double Beta)
{
  m_Tfin = Tfin;
  m_nFonc = n_fonc;
  a =2.;// coefficient pour la fonction u
  twelth=1./12.;
  third=1./3.;
  m_CFL=CFL;
  m_CL=Cond_lim;
  m_nbCellx=nb_cell;
  m_Beta=Beta;
  g_pi=3.141592653589793;
  std::cout<<twelth<<third<<std::endl;
  m_tmp = new Maillage(inf_x,sup_x,inf_y,sup_y,nb_cell,nb_cell,m_CL);
  m_M   = new Maillage(inf_x,sup_x,inf_y,sup_y,nb_cell,nb_cell,m_CL);
  m_M->init_maill(m_nFonc);
  m_tmp->init_maill(m_nFonc);
  Dx=m_M->getDxCell(0);
  h_x=1./Dx;
}

void Diffuse::solve_diffuse()
{
  double t=0;
  int compt=0;
  double dt;
  //double Dx;
  double Val;
  double u_x,u_y;
  int x,y;
  std::vector<double> f_tmp;
  std::vector<double> f_tmp2;
  f_tmp.reserve(m_M->getNbCell());
  f_tmp.resize(m_M->getNbCell()); 
  f_tmp2.reserve(m_M->getNbCell());
  f_tmp2.resize(m_M->getNbCell());  
  dt=Dx*m_CFL*0.5;

# pragma omp parallel for shared(f_tmp)
  for(int i=0;i<m_M->getNbCell();i++){
    f_tmp.at(i)=0;
  }
  while(t<m_Tfin){
# pragma omp parallel for shared(f_tmp2)
    for(int i=0;i<m_M->getNbCell();i++){
      f_tmp2.at(i)= m_M->getValueCell(i);
    }
# pragma omp parallel for shared(f_tmp) private(u_x,u_y,x,y)
    for(int i=0;i<m_M->getNbCell();i++){ 
      switch(m_nFonc){ 
      case 0:
	u_x=2;
	u_y=1;
	break;
      case 1:
	u_x=1;
	u_y=1;
	break;
      case 2:
	u_y=m_M->getBinfCell_x(i)+0.5*Dx;
	u_x=-(m_M->getBinfCell_y(i)+0.5*Dx);
	break;
      case 3:
	double pos_x,pos_y;
	pos_x=g_pi*m_M->getBinfCell_x(i)+0.5*Dx;
	pos_y=g_pi*m_M->getBinfCell_y(i)+0.5*Dx;
	u_x=-sin(pos_y)*cos(pos_y)*sin(pos_x)*sin(pos_x)*cos(g_pi*t/m_Tfin);
	u_y=sin(pos_x)*cos(pos_x)*sin(pos_y)*sin(pos_y)*cos(g_pi*t/m_Tfin);
	break;
      default:
	std::cout<<"mauvaise fonction"<<std::endl;
	break;
      }
      x=i%m_nbCellx;
      y=i/m_nbCellx;
      f_tmp.at(i)=m_M->getValueCell(i)-dt*0.5*psix(x,y,u_x)*h_x
	-dt*0.5*psiy(x,y,u_y)*h_x;
    } 
# pragma omp parallel for shared(f_tmp) private(Val)
    for(int i=0;i<m_M->getNbCell();i++){
      Val=f_tmp.at(i);
      m_M->setValueCell(i,Val);
    }
# pragma omp parallel for shared(f_tmp) private(u_x,u_y,x,y)
    for(int i=0;i<m_M->getNbCell();i++){ 
      switch(m_nFonc){ 
      case 0:
	u_x=2;
	u_y=1;
	break;
      case 1:
	u_x=1;
	u_y=1;
	break;
      case 2:
	u_y=m_M->getBinfCell_x(i)+0.5*Dx;
	u_x=-(m_M->getBinfCell_y(i)+0.5*Dx);
	break;
      case 3:
	double pos_x,pos_y;
	pos_x=g_pi*m_M->getBinfCell_x(i)+0.5*Dx;
	pos_y=g_pi*m_M->getBinfCell_y(i)+0.5*Dx;
	u_x=-sin(pos_y)*cos(pos_y)*sin(pos_x)*sin(pos_x)*cos(g_pi*t/m_Tfin);
	u_y=sin(pos_x)*cos(pos_x)*sin(pos_y)*sin(pos_y)*cos(g_pi*t/m_Tfin);
	break;	
      default:
	std::cout<<"mauvaise fonction"<<std::endl;
	break;
      }
      x=i%m_nbCellx;
      y=i/m_nbCellx;
      f_tmp.at(i)=f_tmp2.at(i)-dt*psix(x,y,u_x)*h_x
	-dt*psiy(x,y,u_y)*h_x;
    } 
# pragma omp parallel for shared(f_tmp) private(Val)
    for(int i=0;i<m_M->getNbCell();i++){
      Val=f_tmp.at(i);
      m_M->setValueCell(i,Val);
    }
    t=t+dt;
    if(compt%100==0){
      std::cout<<t <<std::endl;     
    }
    compt=compt+1; 
    /*if(compt==7000){
      saveMaillage();
      std::cout<<"ok"<<t<<std::endl;
      }  */
  }

}
//*** Ok solution complete ***//
void Diffuse::solution()
{
  switch(m_nFonc){
  case(0):
# pragma omp parallel for
    for(int i=0; i<m_M->getNbCell();i++){
      double mil_x,mil_y;
      mil_x=m_tmp->getBinfCell_x(i)+Dx*0.5-2*m_Tfin;
      mil_y=m_tmp->getBinfCell_y(i)+Dx*0.5-m_Tfin;
      if(mil_y<=0.5*mil_x){
	m_tmp->setValueCell(i,1);
      }
      else{
	m_tmp->setValueCell(i,0);
      }
    }
    break;
  case(1):
# pragma omp parallel for
    for(int i=0; i<m_M->getNbCell();i++){
      double mil_x,mil_y;
      mil_x=m_tmp->getBinfCell_x(i)+Dx*0.5-m_Tfin;
      mil_y=m_tmp->getBinfCell_y(i)+Dx*0.5-m_Tfin;
      if(mil_y*mil_y+mil_x*mil_x<=0.2){
	m_tmp->setValueCell(i,1);
      }
      else{
	m_tmp->setValueCell(i,0);
      }
    }
    break;
  case(2):
# pragma omp parallel for
    for(int i=0; i<m_M->getNbCell();i++){ 
      double mil_x,mil_y;
      double tmp_x,tmp_y;
      mil_x=m_tmp->getBinfCell_x(i)+Dx*0.5;
      mil_y=m_tmp->getBinfCell_y(i)+Dx*0.5;
      tmp_x=mil_y*std::sin(m_Tfin)+mil_x*std::cos(m_Tfin);
      tmp_y=mil_y*std::cos(m_Tfin)-mil_x*std::sin(m_Tfin);
      if(std::pow(tmp_x-0.5,2)+std::pow(tmp_y,2)<=0.15){
	m_tmp->setValueCell(i,1);
      }
      else{
	m_tmp->setValueCell(i,0);
      }
    }
    break;
  default:
    std::cout<<"mauvais choux de fonction"<<std::endl;
    break;

  }
}      
 
/////**** fin de solution****////

//****Calcul des delta Z****//

double Diffuse:: phi_demipp(int i,int j,double a , double b){
  double phi;
  double tmp;
  double grad_tmp;
  double zsup,zsup_tmp;
  tmp=0;
  grad_tmp=a+b;
  if(grad_tmp>0)
    {
      zsup=std::max(m_M->getValxy(i,j),m_M->getValxy(i,j+1));
      zsup_tmp=std::max(m_M->getValxy(i+1,j),m_M->getValxy(i+1,j+1));
      zsup=std::max(zsup,zsup_tmp);
      tmp=(zsup-m_M->getValxy(i,j))/(grad_tmp);
    }
  else if(grad_tmp<0)
    {
      zsup=std::min(m_M->getValxy(i,j),m_M->getValxy(i,j+1));
      zsup_tmp=std::min(m_M->getValxy(i+1,j),m_M->getValxy(i+1,j+1));
      zsup=std::min(zsup,zsup_tmp);
      tmp=(zsup-m_M->getValxy(i,j))/(grad_tmp);
    }
  phi=std::min(m_Beta,tmp);
  return phi;
}
double Diffuse:: phi_demipm(int i,int j,double a, double b){
  double phi;
  double tmp;
  double grad_tmp;
  double zsup,zsup_tmp;
  tmp=0;
  grad_tmp=a-b;
  if(grad_tmp>0)
    {
      zsup=std::max(m_M->getValxy(i,j),m_M->getValxy(i,j-1));
      zsup_tmp=std::max(m_M->getValxy(i+1,j-1),m_M->getValxy(i+1,j));
      zsup=std::max(zsup,zsup_tmp);
      tmp=(zsup-m_M->getValxy(i,j))/(grad_tmp);
    }
  else if(grad_tmp<0)
    {
      zsup=std::min(m_M->getValxy(i,j),m_M->getValxy(i,j-1));
      zsup_tmp=std::min(m_M->getValxy(i+1,j),m_M->getValxy(i+1,j-1));
      zsup=std::min(zsup,zsup_tmp);
      tmp=(zsup-m_M->getValxy(i,j))/(grad_tmp);
    }
  phi=std::min(m_Beta,tmp);
  return phi;
}
double Diffuse:: phi_demimm(int i,int j,double a , double b){
  double phi;
  double tmp;
  double grad_tmp;
  double zsup,zsup_tmp;
  tmp=0;
  grad_tmp=-a-b;
  if(grad_tmp>0)
    {
      zsup=std::max(m_M->getValxy(i,j),m_M->getValxy(i,j-1));
      zsup_tmp=std::max(m_M->getValxy(i-1,j),m_M->getValxy(i-1,j-1));
      zsup=std::max(zsup,zsup_tmp);
      tmp=(zsup-m_M->getValxy(i,j))/(grad_tmp);
    }
  else if(grad_tmp<0)
    {
      zsup=std::min(m_M->getValxy(i,j),m_M->getValxy(i,j-1));
      zsup_tmp=std::min(m_M->getValxy(i-1,j),m_M->getValxy(i-1,j-1));
      zsup=std::min(zsup,zsup_tmp);
      tmp=(zsup-m_M->getValxy(i,j))/(grad_tmp);
    }
  phi=std::min(m_Beta,tmp);
  return phi;
}
double Diffuse:: phi_demimp(int i,int j, double a , double b){
  double phi;
  double tmp;
  double grad_tmp;
  double zsup,zsup_tmp;
  tmp=0;
  grad_tmp=-a+b;
  if(grad_tmp>0)
    {
      zsup=std::max(m_M->getValxy(i,j),m_M->getValxy(i,j+1));
      zsup_tmp=std::max(m_M->getValxy(i-1,j),m_M->getValxy(i-1,j+1));
      zsup=std::max(zsup,zsup_tmp);
      tmp=(zsup-m_M->getValxy(i,j))/(grad_tmp);
    }
  else if(grad_tmp<0)
    {
      zsup=std::min(m_M->getValxy(i,j),m_M->getValxy(i,j+1));
      zsup_tmp=std::min(m_M->getValxy(i-1,j),m_M->getValxy(i-1,j+1));
      zsup=std::min(zsup,zsup_tmp);
      tmp=(zsup-m_M->getValxy(i,j))/(grad_tmp);
    }
  phi=std::min(m_Beta,tmp);
  return phi;
}

double Diffuse:: phi(int i , int j,double grad_x,double grad_y){
  double phi;
  double tmp;
  // tmp=std::max(phi_demipp(i,j),phi_demimp(i,j));
  //phi=std::max(phi_demimm(i,j),phi_demipm(i,j));
  //phi=std::max(phi,tmp);;
  //phi=(phi_demipp(i,j)+phi_demimp(i,j)+phi_demimm(i,j)+phi_demipm(i,j))*0.25;
  tmp=std::min(phi_demipp(i,j,grad_x,grad_y),phi_demimp(i,j,grad_x,grad_y));
  phi=std::min(phi_demimm(i,j,grad_x,grad_y),phi_demipm(i,j,grad_x,grad_y));
  phi=std::min(phi,tmp);
  return phi;

}

/////*****calcul du gradient*****/////


double Diffuse :: gradient_x(int i,int j){// changer la fonction grad  faire grad_x grad_y
  double grad;
  grad=twelth*(m_M->getValxy(i+1,j+1)-m_M->getValxy(i-1,j+1))+
    third*(m_M->getValxy(i+1,j)-m_M->getValxy(i-1,j))+
    twelth*(m_M->getValxy(i+1,j-1)-m_M->getValxy(i-1,j-1));
  return grad*h_x;
}
double Diffuse :: gradient_y(int i,int j){// changer la fonction grad  faire grad_x grad_y
  double grad;
  grad=twelth*(m_M->getValxy(i+1,j+1)-m_M->getValxy(i+1,j-1))+
    third*(m_M->getValxy(i,j+1)-m_M->getValxy(i,j-1))+
    twelth*(m_M->getValxy(i-1,j+1)-m_M->getValxy(i-1,j-1));
  return grad*h_x;
}
////*****fin de calcul *****////

double Diffuse::z_plus(int i , int j,int coor){ /// rassembler zplus et zmoins pour 1 seul calcul de gradient 
  double z;
  double tmp_gradx, tmp_grady;
  tmp_gradx=gradient_x(i,j)*Dx*0.5;
  tmp_grady=gradient_y(i,j)*Dx*0.5;
  if(coor==0){
    z=m_M->getValxy(i,j)-phi(i,j,tmp_gradx,tmp_grady)*tmp_gradx;
  }
  else
    {
      z=m_M->getValxy(i,j)-phi(i,j,tmp_gradx,tmp_grady)*tmp_grady;
    }
  return z;  
}

double Diffuse::z_moins(int i, int j,int coor){
  double z;
  double tmp_gradx, tmp_grady;
  tmp_gradx=gradient_x(i,j)*Dx*0.5;
  tmp_grady=gradient_y(i,j)*Dx*0.5;
  if(coor==0){
    z=m_M->getValxy(i,j)+phi(i,j,tmp_gradx,tmp_grady)*tmp_gradx;
  }
  else{
    z=m_M->getValxy(i,j)+phi(i,j,tmp_gradx,tmp_grady)*tmp_grady;
  }
  return z;  
}

double Diffuse:: psix(int i, int j,double u){
  double psi;
  double psi_m;
  double psi_p;
  if(u>0){
    psi_m=u*z_moins(i-1,j,0);
    psi_p=u*z_moins(i,j,0);
  }
  else{
    psi_m=u*z_plus(i,j,0);
    psi_p=u*z_plus(i+1,j,0);
  }
  psi=psi_p-psi_m;
  return psi;
}

double Diffuse:: psiy(int i, int j,double u){
  double psi;
  double psi_m;
  double psi_p;
  if(u>0){
    psi_m=u*z_moins(i,j-1,1);
    psi_p=u*z_moins(i,j,1);
  }
  else{
    psi_m=u*z_plus(i,j,1);
    psi_p=u*z_plus(i,j+1,1);
  }
  psi=psi_p-psi_m;
  return psi;
}

void Diffuse::saveMaillage()
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
