#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//--------- DECLARAR CONSTANTES GLOBALES ------------//
const double Lx=60, Ly=60;
const int Nx=5, Ny=5, N=Nx*Ny;
const double re=10, epsi=1.0; //valores dentro de la fuerza de Lennard-Jones//
const double K=10000; //cte de hook de fuerza contra la pared circular

//------------------ CONSTANTES DEL PERFL ------------//
const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//----------------- DECLARAR CLASES ----------------//
class Cuerpo;
class Colisionador;

//---- INTERFASE E IMPLEMENTACION DE LAS CLASES ----//
//---------------- CLASE CUERPO --------------------//
class Cuerpo{
private:
  vector3D r,V,F; double m,R,intensidad; 
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void BorreFuerza(){F.load(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void BorreIntensidad(double I0){intensidad=I0;}; //Valor de la intensidad de las colisiones de una particula contra una pared//
  void AdicioneIntensidad(double P){intensidad+=P;}; //Suma de la intensidad para cada partícula//
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  double GetVx(void){return V.x();}; //inline
  double GetIntensidad(void){return intensidad;};
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt);  
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)"; //Solo dibujar las paredes y quitar los granos que las representan// 
}

//---------------- CLASE COLISIONADOR ----------------//
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Grano);
  void CalculeFuerzaPared(Cuerpo & Grano);
  void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2);
};

void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2){
  vector3D r21, n, F; double d;
  r21 = Grano2.r-Grano1.r; d=r21.norm(); n= r21/d;
  F = n*(12*epsi/d)*(pow((re/d),12)-pow((re/d),6)); //Aquí se define la fuerza deLennard Jones que experimentará la partícula //
  Grano2.AdicioneFuerza(F); Grano1.AdicioneFuerza(F*(-1));    
}

void Colisionador::CalculeFuerzaPared(Cuerpo & Grano){ //Método para calcular la interacción con las paredes//
  vector3D Fp; //fuerza de la pared
  double R=50; //radio de la pared circular
  double r_g=Grano.r.norm(); //posición radial del grano

  double s=r_g+Grano.R-R; //interpenetración

  if(s>0){
    Fp=(-K*s/Grano.r.norm())*Grano.r;
    Grano.AdicioneFuerza(Fp);
  }
}

void Colisionador::CalculeFuerzas(Cuerpo * Grano){
  int i,j; vector3D Fg;
  //---------- BORRA TODAS LAS FUERZAS --------//
  for(i=0;i<N;i++){
    Grano[i].BorreFuerza();
  }
  //--- CALCULAR FUERZAS ENTRE PARES DE GRANOS ---//
  for(i=0;i<N;i++){
    CalculeFuerzaPared(Grano[i]);
    for(j=i+1;j<N;j++){
      CalculeFuerzaEntre(Grano[i], Grano[j]);
    }
  }
}

//-------------- FUNCIONES DE ANIMACION Y DIBUJO ----------//
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'd-lennard-liquido.gif'"<<endl; //Nombre d-lennard-(solido/liquido/gas)
  cout<<"unset key"<<endl;
  cout<<"set xrange[-60:60]"<<endl;
  cout<<"set yrange[-60:60]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , 50*cos(t),50*sin(t)";
}
void TermineCuadro(void){
  //cout<<"\n pause 5\n";
  cout<<endl;
}

//------------------- PROGRAMA PRINCIPAL --------------------//
int main(void){
  Cuerpo Grano[N];
  Colisionador Lennard;
  Crandom ran64(1);//----kT=0.05 solido, 0.5 liquido, 10.0 gas
  double m0=1.0, R0=3.0, kT=0.5, V0=sqrt(2*kT/m0);
  int i,ix,iy;
  double t,tdibujo,dt=5e-4;
  double tmax=100;
  double tcuadro=tmax/2000; //cambiar a 2000 para enviar
  double dl=10.0;
  double Theta, Vx;

  InicieAnimacion(); //Se prepara la terminal para realizar una animación en .gif//
  
  //------------------- INICIALIZAR LAS MOLÉCULAS -------------------//
  for(ix=0;ix<Nx;ix++){
    for(iy=0;iy<Ny;iy++){
      Theta=2*M_PI*ran64.r();
      //--------------------(       x0,       y0,          Vx0,          Vy0, m0,R0)
      Grano[Nx*iy+ix].Inicie(-20+iy*dl,-20+ix*dl,V0*cos(Theta),V0*sin(Theta), m0,R0);
    }
  }
  
  //-------------------- INICIA EL CICLO PRINCIPAL ---------------//
  for(t=0,tdibujo=0 ; t<tmax ; t+=dt,tdibujo+=dt){
    
    //------------------- DIBUJO DE CADA CUADRO --------------------//
    if(tdibujo>tcuadro){
      
      InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      
      tdibujo=0;
    }    
    
    //---------------------- MUEVASE POR PERFL -----------------//
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);
    Lennard.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Lennard.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chiepsilon);
    Lennard.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Lennard.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);
  }
  
  return 0;
}
