#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

//Constantes globales

const int N=5;
const double g=980;//g en cm/s2
const double k=1e10; //cte de fuerza de herz
const double alpha=3/2; //exp de fuerza de herz

//constantes de PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//Declaración de las clases
class Cuerpo;
class Colisionador;

//---------- Clase Cuerpo --------------
class Cuerpo{
private:
  //vector3D r,V,F;
  double theta,omega,tau;
  double m,R,l,x0,I;
public:
  void Inicie(double theta0,double omega0,double l0,
	      double x00,double m0,double R0);
  void BorreTorques(void){tau=0;};
  void SumeTorques(double tau1){tau+=tau1;};
  void Mueva_theta(double dt,double coeficiente);
  void Mueva_omega(double dt,double coeficiente);
  double Getx(void){return x0+l*sin(theta);}; //Inline
  double Gety(void){return -l*cos(theta);}; //Inline
  void Dibujese(void);
  friend class Colisionador;
};
void Cuerpo::Inicie(double theta0,double omega0,double l0,double x00,double m0,double R0){
  theta=theta0; omega=omega0;
  m=m0; R=R0; l=l0; x0=x00; I=m*l*l;
}
void Cuerpo::Mueva_theta(double dt,double coeficiente){
  theta+=omega*(dt*coeficiente);
}
void Cuerpo::Mueva_omega(double dt,double coeficiente){
  omega+=tau*(dt*coeficiente/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+"<<R<<"*sin(t)";
  cout<<" , "<<(l-R)*sin(theta)<<"*t/7+"<<x0<<","<<-(l-R)*cos(theta)<<"*t/7";
}
//---------- Clase Colisionador --------------
class Colisionador{
private:
public:
  void CalculeTorques(Cuerpo * Pendulo);
  void CalculeTorquesEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2);    
};
void Colisionador::CalculeTorques(Cuerpo * Pendulo){
  int i,j;
  //Borrar torques y sume de gravedad
  for(i=0;i<N;i++){
    Pendulo[i].BorreTorques();
    Pendulo[i].SumeTorques(-Pendulo[i].m*g*Pendulo[i].l*sin(Pendulo[i].theta));
  }
  //Calcular las fuerzas entre todas las parejas de planetas
  for(i=1;i<N;i++){
    CalculeTorquesEntre(Pendulo[i],Pendulo[i-1]);
  }
}
void Colisionador::CalculeTorquesEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2){
  double r21=sqrt(pow((Pendulo2.Getx()-Pendulo1.Getx()),2)
		 +pow((Pendulo2.Gety()-Pendulo1.Gety()),2));
  double s=Pendulo2.R+Pendulo1.R-r21;
  if(s>0){
    double tau=k*pow(s,alpha)/r21*((Pendulo2.x0-Pendulo1.x0)*Pendulo2.l*cos(Pendulo2.theta)+
				   Pendulo2.l*Pendulo1.l*sin(Pendulo2.theta-Pendulo1.theta));
    Pendulo2.SumeTorques(tau);
    tau=k*pow(s,alpha)/r21*(-(Pendulo2.x0-Pendulo1.x0)*Pendulo1.l*cos(Pendulo1.theta)-
			    Pendulo2.l*Pendulo1.l*sin(Pendulo2.theta-Pendulo1.theta));
    Pendulo1.SumeTorques(tau);
  }
}

//----------- Funciones Globales -----------

void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl; 
  //cout<<"set output 'Balon.gif'"<<endl;
  cout<<"unset key"<<endl; //para que no aparezca title
  cout<<"set xrange[-15:31]"<<endl; //rangos
  cout<<"set yrange[-25:5]"<<endl;
  cout<<"set size ratio -1"<<endl; //relacion en pantalla 1:1
  cout<<"set parametric"<<endl; //que sea una figura parametrizca x(t), y(t)
  cout<<"set trange [0:7]"<<endl; //rango de t
  cout<<"set isosamples 12"<<endl; //rango de t dividido en 12 puntos
  cout<<"set terminal dumb"<<endl;
}

void InicieCuadro(void){
    cout<<"plot 0,0 ";
}

void TermineCuadro(void){
  cout<<"\n"<<"pause 0.01";
  cout<<endl;
}

int main(){
  Cuerpo Pendulo[N];
  Colisionador Newton;
  double m0=50, l=12, R=2;
  double theta=-M_PI/9;
  double omega0=0, T=2*M_PI*sqrt(l/g);
  double t,tmax=3*T,dt=1e-5;
  double tdibujo,tcuadro=T/1000;
  int i;
  
  //---------------(theta,omega,l,x0,m0,R0)
  Pendulo[0].Inicie(theta,omega0,l,0,m0,R);
  Pendulo[1].Inicie(theta,omega0,l,4,m0,R);
  for(i=2;i<N;i++){
    Pendulo[i].Inicie(0,0,l,i*2*R,m0,R);
  }

  InicieAnimacion();
 
  
  for(t=0,tdibujo=0; t<tmax; t+=dt,tdibujo+=dt){
    //Dibujar
    if(tdibujo>tcuadro){
      InicieCuadro();
      for(i=0;i<N;i++) Pendulo[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    
    //cout<<Pendulo[1].Getx()<<" "<<Pendulo[1].Gety()<<endl;
    // Mover por PEFRL
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Zeta);
    Newton.CalculeTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Chi);
    Newton.CalculeTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Lambda);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Coeficiente2);
    Newton.CalculeTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Lambda);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Chi);
    Newton.CalculeTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Zeta);   
  }
  
  return 0;
}
