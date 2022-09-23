// Simular el movimiento de N moleculas en un gas 2D
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "../../vector.h"
#include "../../Random64.h"
using namespace std;

//---- declarar constantes ---
const double Lx=160, Ly=60;
const int Ns=80;

//--- declarar clases -----
class Cuerpo;

//---- interface e implementacion de clases ----
//---- clase cuerpo ---
class Cuerpo{
private:
  vector3D r,V,F; double m,R; double theta,omega,tau,I;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,
	      double theta0,double omega0,double m0,double R0);
  void BorreFuerza(){F.load(0,0,0); tau=0;};
  void AdicioneFuerza(vector3D F0,double tau0){F+=F0; tau+=tau0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(bool p,int lt);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  double GetR(void){return R;}; //inline
  double Getheta(void){return theta;}; //inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,
		    double theta0,double omega0,double m0,double R0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
  theta=theta0; omega=omega0; I=2.0/5*m*R*R;
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt); theta+=omega*(Coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m); omega+=tau*(Coeficiente*dt/I);
}
void Cuerpo::Dibujese(bool p, int lt){
  if(p==false){
    cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
	<<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
  }
  else{
    cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) lc "<<lt<<" , "
	<<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t lc "<<lt;
  }
}

//-----------  Programa Principal --------------
int main(void){
  ifstream fin("f_piladearena(seed2).txt");
  int N;fin>> N;
  std::vector<Cuerpo> syst(N);

  //Inicializar los granitos
  for(int i=0;i<syst.size();i++){
    double x,y,R,theta;
    fin>> x;fin>> y;fin>> theta;fin>> R;
    syst[i].Inicie(x,y,0,0,theta,0,0,R);
  }

  int Nx=40; //particiones del eje x
  double dx=Lx/Nx;

  for(int i=0;i<Nx;i++){
    double xmax=0.0;
    double ymax=0.0;
    for(int j=0;j<syst.size();j++){
      if(syst[j].Getx()>=dx*i && syst[j].Getx()<dx*(i+1)){
	if(syst[j].Gety()>ymax){
	  ymax=syst[j].Gety();
	  xmax=syst[j].Getx();
	  syst.erase(syst.begin()+j);
	  j--;
	}
      }
    }
    cout<<xmax<<"\t"<<ymax<<"\n";
  }
  
  return 0;
}
