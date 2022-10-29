#include<iostream>
#include<cmath>
#include<vector>
#include"vector.h"

//Constantes globables de simulación
const double G=1; //Cte gravitacional
//constantes globales de integración
const double Xi=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e00;
const double Chi=-0.6626458266981849e-1;
const double coef1=(1-2*Lambda)/2;
const double coef2=1-2*(Chi+Xi);

//-----------------------------------------------------------------------------------

/*Clase Cuerpo: objeto que guarda las coordenadas r,v,F,m y R (radio) de cada cuerpo
  del sistema
  rutinas: Init=>inicia las variables del cuerpo
           print=>printea en pantalla la posición del cuerpo
	   printv=>printea en pantalla la velocidad del cuerpo
   	   printF=>printea en pantalla la Fuerza del cuerpo
*/

class Cuerpo{
private:
  vector3D r,v,F;
  double m,R;
  
public:
  //inicia las variables del cuerpo
  void Init(double x0,double y0,double z0,
	    double vx0,double vy0,double vz0,
	    double m0,double R0);
  //printea las variables del cuerpo
  void printr();
  void printv();
  void printF();
  void dibujese();

  friend class Colisionador;
};

//Implementaciones de Cuerpo

void Cuerpo::Init(double x0,double y0,double z0,
		  double vx0,double vy0,double vz0,
		  double m0,double R0){
  r.load(x0,y0,z0); v.load(vx0,vy0,vz0);
  m=m0; R=R0;
}

void Cuerpo::printr(){
  std::cout<<r.x()<<"\t"<<r.y()<<"\t"<<r.z()<<"\t";
}

void Cuerpo::printv(){
  std::cout<<v.x()<<"\t"<<v.y()<<"\t"<<v.z()<<"\t";
}

void Cuerpo::printF(){
  std::cout<<F.x()<<"\t"<<F.y()<<"\t"<<F.z()<<"\t";
}

void Cuerpo::dibujese(){
  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(2*pi*t),"<<r.y()<<"+"<<R<<"*sin(2*pi*t)";
}

//-----------------------------------------------------------------------------------

/*Clase Colisionador: accede a las variables privadas de los cuerpos de un sistema
  (vector de cuerpos) y calcula las fuerzas de cada cuerpo. También evoluciona un
  paso en el teimpo usando PEFRL.
 */

class Colisionador{
public:
  //calcula y asiga a cada cuerpo del sistema la fuerza que debe sentir
  void Fuerza_syst(std::vector<Cuerpo> &syst,double omega);
  //Calcula solo las fuerzas ficticias (de inercia) que debe sentir cada cuerpo
  void Ficti_syst(std::vector<Cuerpo> &syst,double omega);
  //calcula solo la fuerza entre dos cuerpos del sistema
  void F_entre(Cuerpo &c1,Cuerpo &c2);
  //Ejecuta el paso en el tiempo con una integración PEFRL
  void Paso_syst(std::vector<Cuerpo> &syst,double dt,double omega);
};

//Implementaciones de Colisionador

void Colisionador::Fuerza_syst(std::vector<Cuerpo> &syst,double omega){
  for(int i=0;i<syst.size();i++){
    syst[i].F.load(0,0,0);
  }
  for(int i=0;i<syst.size();i++){
    for(int j=i+1;j<syst.size();j++){
      F_entre(syst[i],syst[j]);
    }
  }
  Ficti_syst(syst,omega);
}

void Colisionador::Ficti_syst(std::vector<Cuerpo> &syst,double omega){
  vector3D v_o;
  v_o.load(0,0,omega);
  for(int i=0;i<syst.size();i++){
    syst[i].F+=-2*syst[i].m*(v_o^syst[i].v)-syst[i].m*(v_o^(v_o^syst[i].r)); 
  }
}

void Colisionador::F_entre(Cuerpo &c1,Cuerpo &c2){
  vector3D r21=c2.r-c1.r;
  double F=G*c1.m*c2.m*std::pow(r21.norm(),-3);
  c1.F+=F*r21; c2.F+=-F*r21;
}

void Colisionador::Paso_syst(std::vector<Cuerpo> &syst,double dt,double omega){
  for(int i=0;i<syst.size();i++){
    syst[i].r+=Xi*dt*syst[i].v;
  }
  Fuerza_syst(syst,omega);
  for(int i=0;i<syst.size();i++){
    syst[i].v+=coef1*dt*syst[i].F/syst[i].m;
  }
  for(int i=0;i<syst.size();i++){
    syst[i].r+=Chi*dt*syst[i].v;
  }
  Fuerza_syst(syst,omega);
  for(int i=0;i<syst.size();i++){
    syst[i].v+=Lambda*dt*syst[i].F/syst[i].m;
  }
  for(int i=0;i<syst.size();i++){
    syst[i].r+=coef2*dt*syst[i].v;
  }
  Fuerza_syst(syst,omega);
  for(int i=0;i<syst.size();i++){
    syst[i].v+=Lambda*dt*syst[i].F/syst[i].m;
  }
  for(int i=0;i<syst.size();i++){
    syst[i].r+=Chi*dt*syst[i].v;
  }
  Fuerza_syst(syst,omega);
  for(int i=0;i<syst.size();i++){
    syst[i].v+=coef1*dt*syst[i].F/syst[i].m;
  }
  for(int i=0;i<syst.size();i++){
    syst[i].r+=Xi*dt*syst[i].v;
  }
}
//----------------------------------------------------------------------------------
//funciones de graficación
void set_graph(){
  std::cout<<"set terminal pdf\n";
  std::cout<<"set output 'troyanos_rot.pdf'\n";
  std::cout<<"unset key\n";
  std::cout<<"set size ratio -1\n";
  std::cout<<"set xrange[-100:1100]\n";
  std::cout<<"set yrange[-100:100]\n";
  std::cout<<"set parametric\n";
  std::cout<<"set trange[0:1]\n";
  std::cout<<"set isosamples 12\n";
}

void init_graph(){
  std::cout<<"plot 0,0";
}

void end_graph(){
  std::cout<<"\n";
}


//-----------------------------------------------------------------------------------

int main(){
  int N=2; //número de cuerpos
  double dt=0.1; //paso de tiempo
  
  std::vector<Cuerpo> syst(N); //creación de sistema=>un vector de cuerpos
  Colisionador Gravity; //el colisionador=>la gravedad

  //condiciones iniciales
  double m0=1047,m1=1,r0=1000,omega,V0,V1,T;
  double M=m0+m1;
  double x0=-m1*r0/M;
  double x1=m0*r0/M;
  omega=std::sqrt(G*M*pow(r0,-3));
  T=2*M_PI/omega;

  //condiciones de graficación
  double tgraph=T;
  double taux;
    
  
  //----------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0)
  syst[0].Init( x0, 0.0, 0.0, 0.0, 0.0, 0.0, m0, 50);
  syst[1].Init( x1, 0.0, 0.0, 0.0, 0.0, 0.0, m1, 25);

  
  set_graph();
  init_graph();

  for(double t=0,taux=0;t<20*T;t+=dt,taux+=dt){
    if(taux>tgraph){
      for(int i=0;i<syst.size();i++){
	syst[i].dibujese();
      }
      taux=0;
    }
    Gravity.Paso_syst(syst,dt,omega);
  }
  
  end_graph();
  
  return 0;
}
