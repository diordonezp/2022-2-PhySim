#include<iostream>
#include<cmath>
#include"vector.h"

//CONSTANTES GLOBALES
const double g=9.8;
const double GM=1;
//constantes de pefrl
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Lambda);

//DECLARACIÓN DE CLASES
class Cuerpo; // se declara la interfaz, luego se declara el interior
//--------------Clase Cuerpo--------------------------
class Cuerpo{ 
private: /*datos que no son accesibles, solo son importantes las ordenes*/
  vector3D r,V,F;
  double m,R; 
  
public:  //datos que son accesibles 
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0);
  void CalculeFuerza(void); //para indicar que está vacío
  void Arranque(double dt);
  void Muevar(double dt,double theta);
  void MuevaV(double dt,double theta);
  double Getx(void){return r.x();}; //Función Inline
  double Gety(void){return r.y();};
  double Getz(void){return r.z();}; 
};
//aquí acabaría un eventual .h

//----------------------------------------------------
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0){
  /*tengo acceso a todas las variables privada
    porque estoy con Cuerpo::*/
  
  r.load(x0,y0,z0); V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
}

void Cuerpo::CalculeFuerza(void){
  F=-GM*m/pow(r.norm(),3)*r;
}

void Cuerpo::Muevar(double dt,double theta){
  r+=V*(dt*theta);
}

void Cuerpo::MuevaV(double dt,double theta){
  V+=F*(dt*theta/m);
}

//--------------Funciones Globales--------------------
int main(){
  Cuerpo Planeta; //Ejemplares de la clase Cuerpo (Instance)
  double t, dt=1;
  double omega, r0, V0, T;
  r0=10; omega=sqrt(GM*pow(r0,-3)); V0=omega*r0; T=2*M_PI/omega;
  
  //----------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0)
  Planeta.Inicie( r0, 0,0, 0, V0/2,0, 0.45, 0.15);
  
  for(t=0;t<2.1*T;t+=dt){
    std::cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<"\n";
    //mover usando forest_ruth
    Planeta.Muevar(dt,Zeta);
    Planeta.CalculeFuerza();
    Planeta.MuevaV(dt,Coeficiente1);
    Planeta.Muevar(dt,Chi);
    Planeta.CalculeFuerza();
    Planeta.MuevaV(dt,Lambda);
    Planeta.Muevar(dt,Coeficiente2);
    Planeta.CalculeFuerza();
    Planeta.MuevaV(dt,Lambda);
    Planeta.Muevar(dt,Chi);
    Planeta.CalculeFuerza();
    Planeta.MuevaV(dt,Coeficiente1);
    Planeta.Muevar(dt,Zeta);
  }
  
  return 0;
}
