// Simular el movimiento de N moleculas en un gas 2D
#include <iostream>
#include <fstream>
#include <cmath>
#include "../../vector.h"
#include "../../Random64.h"
using namespace std;

//---- declarar constantes ---
const double K=1.0e4;
const double Lx=160, Ly=60;
const int N=200, Ns=80, Ntot=N+Ns+3;

const double g=9.8, Gamma=150, Kcundall=500, mu=0.4;

const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//--- declarar clases -----
class Cuerpo;
class Colisionador;

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

//--- clase Colisionador ----
class Colisionador{
private:
  double xCundall[Ntot][Ntot],sold[Ntot][Ntot];
public:
  void Inicie(void);
  void CalculeFuerzas(Cuerpo * Grano,double dt,int Nlive);
  void CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2
			  ,double & x_Cundall,double & s_old,double dt);
};
void Colisionador::Inicie(void){
  int i,j; //j>i
  for(i=0;i<Ntot;i++)
    for(j=0;j<Ntot;j++)
      xCundall[i][j]=sold[i][j]=0;
}
void Colisionador::CalculeFuerzas(Cuerpo * Grano,double dt,int Nlive){
  int i,j; 
  //--- Borrar todas las fuerzas ---
  for(i=0;i<Ntot;i++)
    Grano[i].BorreFuerza();
  //--- Añadir fuerza de gravedad ---
  vector3D F0;
  for(i=0;i<Nlive;i++){
    F0.load(0,-Grano[i].m*g,0);
    Grano[i].AdicioneFuerza(F0,0);
  }
  //--- Calcular Fuerzas entre pares de planetas ---
  for(i=0;i<Nlive;i++){
    for(j=i+1;j<Nlive;j++){
      CalculeFuerzaEntre(Grano[i],Grano[j],xCundall[i][j],sold[i][j],dt);
    }
    for(j=N;j<Ntot;j++){
      CalculeFuerzaEntre(Grano[i],Grano[j],xCundall[i][j],sold[i][j],dt);
    }
  }
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2
				      ,double & x_Cundall,double & s_old,double dt){
  
  //Cantidades generales para saber si hay contacto
  vector3D r21=Grano2.r-Grano1.r; double R1=Grano1.R, R2=Grano2.R;
  double d=r21.norm(), s=R1+R2-d;
  
  if(s>0){//si hay contacto
    
    //Variables a calcular
    vector3D F2,F1,tau2,tau1;
    
    //Vectores unitarios
    vector3D n=r21*(1.0/d),t,k; t.load(n.y(),-n.x(),0); k.load(0,0,1);

    //Velocidades relativas
    vector3D V21=Grano2.V-Grano1.V;
    vector3D Rw; Rw.load(0,0,R1*Grano1.omega+R2*Grano2.omega);
    vector3D Vc=V21-(Rw^n); double Vn=Vc*n, Vt=Vc*t;

    //Fn (Fuerza de Hertz-Kuwabara-Kono) 
    double m12=(Grano1.m*Grano2.m)/(Grano1.m+Grano2.m); 
    double Fn=K*pow(s,1.5)-Gamma*m12*sqrt(s)*Vn; if(Fn<0) Fn=0;
    
    //Ft (Fuerza de Cundall)
    x_Cundall+=Vt*dt; double Ft=-Kcundall*x_Cundall; double Ftmax=mu*fabs(Fn);
    if(fabs(Ft)>Ftmax) Ft=Ft/fabs(Ft)*Ftmax;
    
    //Calcula y Cargue las fuerzas
    F2=n*Fn+t*Ft; tau2=((n*(-R2))^F2); F1=F2*(-1); tau1=((n*R1)^F1);
    Grano2.AdicioneFuerza(F2,tau2*k);   Grano1.AdicioneFuerza(F1,tau1*k);
  }

  if(s_old>=0 && s<0) x_Cundall=0;
    
  s_old=s;
}

//-----------  Programa Principal --------------
int main(void){
  Cuerpo Grano[Ntot];
  Colisionador Hertz; Hertz.Inicie();
  Crandom Ran64(1);
  Crandom ran64(1);
  double m0=1;
  int i,Nlive;
  double cuadros=5,t,tdibujo,dt=1e-3,tmax=cuadros*sqrt(Ly/g),tcuadro=tmax/(4*cuadros);
  double Omega0,OmegaMax=8.0;
  double R0,RMin=1.6,RMax=2.4;
  
  //Inicializar las paredes
  double Rpared=100*Lx, Mpared=100*m0, Rs=Lx/(2.0*Ns);
  //------------------(  x0,       y0,Vx0,Vy0,theta0,omega0,m0,R0)
  Grano[N+Ns+0].Inicie(Lx/2,Ly+Rpared,  0,  0,     0,     0,Mpared,Rpared); //Pared de arriba
  Grano[N+Ns+1].Inicie(Lx+Rpared,Ly/2,  0,  0,     0,     0,Mpared,Rpared); //Pared derecha
  Grano[N+Ns+2].Inicie(  -Rpared,Ly/2,  0,  0,     0,     0,Mpared,Rpared); //Pared izquierda
  
  //Inicializar el suelo de granitos
  for(i=N;i<N+Ns;i++){
    Grano[i].Inicie((1+2*(i-N))*Rs,0.0,0.0,0.0,0.0,0.0,m0,Rs);
  }
  
  //Inicializar las moléculas
  for(i=0;i<N;i++){
    Omega0=2*OmegaMax*Ran64.r()-OmegaMax;
    R0=(RMax-RMin)*ran64.r()+RMin;
    //---------------(x0,     y0,Vx0,Vy0,theta0,omega0,m0,R0)
    Grano[i].Inicie(Lx/2,Ly-2*R0,  0,  0,     0,Omega0,m0,R0);//OJO
  }
  
  for(Nlive=1;Nlive<N;Nlive++){
    for(t=0; t<tmax; t+=dt){
      //--- Muevase por PEFRL ---
      for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,epsilon);
      Hertz.CalculeFuerzas(Grano,dt,Nlive);
      for(i=0;i<Nlive;i++)Grano[i].Mueva_V(dt,lambda2);
      for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,chi);
      Hertz.CalculeFuerzas(Grano,dt,Nlive);
      for(i=0;i<Nlive;i++)Grano[i].Mueva_V(dt,lambda);
      for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,chiepsilon);
      Hertz.CalculeFuerzas(Grano,dt,Nlive);
      for(i=0;i<Nlive;i++)Grano[i].Mueva_V(dt,lambda);
      for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,chi);
      Hertz.CalculeFuerzas(Grano,dt,Nlive);
      for(i=0;i<Nlive;i++)Grano[i].Mueva_V(dt,lambda2);
      for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,epsilon); 
    }   
  }
  Nlive=N;
  
  for(t=0; t<2*tmax; t+=dt){   
    //--- Muevase por PEFRL ---
    for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,epsilon);
    Hertz.CalculeFuerzas(Grano,dt,Nlive);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano,dt,Nlive);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,chiepsilon);
    Hertz.CalculeFuerzas(Grano,dt,Nlive);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano,dt,Nlive);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,epsilon); 
  }

  ofstream fout("../possimul/f_piladearena(seed2).txt");
  fout<<N<<"\n";
  for(i=0;i<N;i++){
    fout<<Grano[i].Getx()<<"\t"<<Grano[i].Gety()<<"\t"<<Grano[i].Getheta()<<"\t"<<Grano[i].GetR()<<"\n";
  }
  fout.close();
  return 0;
}
