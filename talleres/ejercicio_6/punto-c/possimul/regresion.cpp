// Simular el movimiento de N moleculas en un gas 2D
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "../../vector.h"
#include "../../Random64.h"
using namespace std;

/*---ctes globales---*/
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

//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  cout<<"set terminal pdf"<<endl; 
  cout<<"set output 'regresion(seed2).pdf'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-1:"<<Lx+1<<"]"<<endl;
  cout<<"set yrange[-1:"<<Ly+1<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;
  //cout<<"set terminal dumb"<<endl; 
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
  cout<<endl;
}

//-----------  Programa Principal --------------
int main(){
  ifstream fin("f_piladearena(seed2).txt");
  int N;fin>> N;
  std::vector<Cuerpo> syst(N);
  std::vector<Cuerpo> piso(Ns);

  //Inicializar el suelo de granitos
  double Rs=Lx/(2.0*Ns);
  for(int i=0;i<piso.size();i++){
    piso[i].Inicie((1+2*i)*Rs,0.0,0.0,0.0,0.0,0.0,0,Rs);
  }

  //Inicializar los granitos
  for(int i=0;i<syst.size();i++){
    double x,y,R,theta;
    fin>> x;fin>> y;fin>> theta;fin>> R;
    syst[i].Inicie(x,y,0,0,theta,0,0,R);
  }


  
  InicieAnimacion();
  InicieCuadro();
  for(int i=0;i<syst.size();i++){
    syst[i].Dibujese(false,1);
  }
  for(int i=0;i<piso.size();i++){
    piso[i].Dibujese(true,10);
  }
  TermineCuadro();
  
  return 0;
}
