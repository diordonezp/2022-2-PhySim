#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

//------DECLARACIÓN DE CONSTANTES GLOBALES------//
const int N = 2; //Número de partículas//
const int epsi = 1.0; //Constante de la fuerza de Lenard-Jones//
const double re = 10.0;//Distancia de equilibrio de la fuerza de Lenard-Jones//

//---------CONSTANTES PERFL------------//
const double Zeta = 0.1786178958448091;
const double Lambda = -0.2123418310626054;
const double Chi = -0.6626458266981849e-1;

const double Lambda1 = 0.5 - Lambda;
const double Zetachi = 1 - 2*(Chi+Zeta);

//------------DECLARACIÓN DE CLASES----------//
class Cuerpo;
class Colisionador;

//-------------CLASE CUERPO------------------//
class Cuerpo{
private:
  vector3D r,V,F;
  double m, R;

public:
  void Inicie(vector3D r0, vector3D V0, double m0, double R0);
  void SumeFuerza(vector3D F0);
  void BorreFuerza(void);
  void MuevaR(double dt, double Coeficiente);
  void MuevaV(double dt, double Coeficiente);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  friend class Colisionador;
};

void Cuerpo::Inicie(vector3D r0, vector3D V0, double m0, double R0){
  r=r0; V=V0; m=m0; R=R0;
}

void Cuerpo::BorreFuerza(void){
  F.load(0,0,0);
}

void Cuerpo::SumeFuerza(vector3D F0){
  F+=F0;
}

void Cuerpo::MuevaR(double dt, double a){
  r+= V*a*dt;
}

void Cuerpo::MuevaV(double dt, double a){
  V+= F*(a*dt/m);
}


//-----------------CLASE COLISIONADOR----------------//
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Grano);
  void CalculeFuerzasEntre(Cuerpo & Grano1, Cuerpo & Grano2);
};

void Colisionador::CalculeFuerzasEntre(Cuerpo & Grano1, Cuerpo & Grano2){
  vector3D r21, n, F; double d;
  r21 = Grano2.r-Grano1.r; d=r21.norm(); n= r21/d;
  F = n*(12*epsi/d)*(pow((re/d),12)-pow((re/d),6)); //Aquí se define la fuerza deLennard Jones que experimentará la partícula //
  Grano2.SumeFuerza(F); Grano1.SumeFuerza(F*(-1));    
}

void Colisionador::CalculeFuerzas(Cuerpo * Grano){
  int i,j; vector3D F;
  for(i=0;i<N;i++){
    Grano[i].BorreFuerza();
  }
  for(i=0;i<N;i++){
    for(j=i+1;j<N;j++){
      CalculeFuerzasEntre(Grano[i],Grano[j]); //La fuerza se da entre dos partículas //
    }
  }
}

//------------------Función Principal-----------------------//
int main(void){
  Cuerpo Grano[N];
  Colisionador Lennard;
  double m=1.0, R=3.0, x=10, y=0, kT=0.5, V0=sqrt(2*kT/m), Vx=V0, Vy=0; //condiciones iniciales//
  double t, tmax=100, dt=0.1;
  vector3D r0, r, Vo, V;

  r.load(x,y,0); r0.load(0,0,0); //se fija una partícula en el origen que no experimentará ninguna fuerza //
  V.load(Vx,Vy,0);  Vo.load(0,0,0); 
  
  Grano[0].Inicie(r0,Vo,0,1);//se establecen las condiciones iniciales para la partícula que si se mueve y la que no//
  Grano[1].Inicie(r,V,m,R);
  Lennard.CalculeFuerzas(Grano);
  
  for(t=0;t<tmax;t+=dt){ //En la aplicación del PERFL solo se considera la partícula en movimiento y no la que está en el origen  //
    cout<<t<<" "<<Grano[1].Getx()<<endl; // Se extrae la posición en x en función del tiempo para la partícula en movimiento//
    Grano[1].MuevaR(dt,Zeta); //1 paso
    Lennard.CalculeFuerzas(Grano);
    Grano[1].MuevaV(dt,Lambda1); //2 Paso
    Grano[1].MuevaR(dt,Chi); //3 Paso
    Lennard.CalculeFuerzas(Grano);
    Grano[1].MuevaV(dt,Lambda); //4 Paso
    Grano[1].MuevaR(dt,Zetachi); //5 Paso
    Lennard.CalculeFuerzas(Grano);
    Grano[1].MuevaV(dt,Lambda); //6 Paso
    Grano[1].MuevaR(dt,Chi); //7 Paso
    Lennard.CalculeFuerzas(Grano);
    Grano[1].MuevaV(dt,Lambda1); //8 Paso
    Grano[1].MuevaR(dt,Zeta);// 9 Paso
  }
  return 0;
}
