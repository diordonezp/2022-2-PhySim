#include <iostream>
#include <cmath>

using namespace std;
//CONSTANTES GLOBALES
const double g=9.8;
const double GM=1;
//DECLARACIÓN DE CLASES
class Cuerpo; // se declara la interfaz, luego se declara el interior
//--------------Clase Cuerpo--------------------------
class Cuerpo{ 
private: //datos que no son accesibles, solo son importantes las ordenes  
  double x,y,Vx,Vy,Fx,Fy,m,R; 
public:  //datos que son accesibles 
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void CalculeFuerza(void); //para indicar que está vacío
  void Muevase(double dt);
  double Getx(void){return x;}; //Función Inline
  double Gety(void){return y;}; 
};
//aquí acabaría un eventual .h

//----------------------------------------------------
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){// tengo acceso a todas las variables privadas porque estoy con Cuerpo::
  x=x0; y=y0; Vx=Vx0; Vy=Vy0; m=m0; R=R0;
}

void Cuerpo::CalculeFuerza(void){
  double r;
  r=sqrt(x*x+y*y);
  Fx=-GM*m*x/pow(r,3); Fy=-GM*m*y/pow(r,3);
}

void Cuerpo::Muevase(double dt){
   x+=Vx*dt;    y+=Vy*dt;
  Vx+=Fx/m*dt; Vy+=Fy/m*dt;
}
//--------------Funciones Globales--------------------
int main(){
  Cuerpo Planeta; //Ejemplares de la clase Cuerpo (Instance)
  double t, dt=0.0001;
  double omega, r0, V0, T;
  r0=10; omega=sqrt(GM*pow(r0,-3)); V0=omega*r0; T=2*M_PI/omega;
  //----------(x0,y0,Vx0,Vy0,m0,R0)
  Planeta.Inicie( r0, 0, 0, V0/2, 0.45, 0.15);
  for(t=0;t<1.1*T;t+=dt){
    cout<<Planeta.Getx()<<" "<<Planeta.Gety()<<endl;
    Planeta.CalculeFuerza();
    Planeta.Muevase(dt);
  }
  return 0;
}
