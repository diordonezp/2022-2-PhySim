#include<iostream>
#include<cmath>
#include"vector.h"

using namespace std;
//CONSTANTES GLOBALES
const double g=9.8;
const double GM=1;
//DECLARACIÓN DE CLASES
class Cuerpo; // se declara la interfaz, luego se declara el interior
//--------------Clase Cuerpo--------------------------
class Cuerpo{ 
private: /*datos que no son accesibles, solo son importantes las ordenes*/
  vector3D r,rold,V,F;
  double m,R; 
  
public:  //datos que son accesibles 
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0);
  void CalculeFuerza(void); //para indicar que está vacío
  void Arranque(double dt);
  void Muevase(double dt);
  double Getx(void){return r.x();}; //Función Inline
  double Gety(void){return r.y();};
  double Getz(void){return r.z();}; 
};
//aquí acabaría un eventual .h

//----------------------------------------------------
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0){/*tengo acceso a todas las variables privada
													  porque estoy con Cuerpo::*/
  r.load(x0,y0,z0); V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
}

void Cuerpo::CalculeFuerza(void){
  F=-GM*m/pow(r.norm(),3)*r;
}

void Cuerpo::Arranque(double dt){
  rold=r-V*dt+F*(dt*dt)/(2*m);
}

void Cuerpo::Muevase(double dt){
  vector3D rnew;
  rnew=2*r-rold+F*(dt*dt/(2*m));
  V=(r-rold)/(2*dt);
  rold=r;
  r=rnew;
}
//--------------Funciones Globales--------------------
int main(){
  Cuerpo Planeta; //Ejemplares de la clase Cuerpo (Instance)
  double t, dt=0.01;
  double omega, r0, V0, T;
  r0=10; omega=sqrt(GM*pow(r0,-3)); V0=omega*r0; T=2*M_PI/omega;

  //----------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0)
  Planeta.Inicie( r0, 0,0, 0, V0/2,0, 0.45, 0.15);
  Planeta.CalculeFuerza();
  Planeta.Arranque(dt);
  for(t=0;t<2.1*T;t+=dt){
    cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<"\n";
    Planeta.CalculeFuerza();
    Planeta.Muevase(dt);
  }
  
  return 0;
}
