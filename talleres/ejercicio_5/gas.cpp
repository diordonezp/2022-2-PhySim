#include<iostream>
#include<cmath>
#include<vector>
#include"vector.h"

//Constantes globables de simulación
const double K=1e6; //Cte de fuerza de hertz
const double Gamma=100;//Cte gamma de fuerza plástica
const double g=9.8; //aceleración gravitacional

//constantes globales de integración
const double Xi=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e00;
const double Chi=-0.6626458266981849e-1;
const double coef1=(1-2*Lambda)/2;
const double coef2=1-2*(Chi+Xi);

//-----------------------------------------------------------------------------------

/*Clase Cuerpo: objeto que guarda las coordenadas r,v,F,m y R (radio) de cada cuerpo
  del sistema, theta el ángulo del cuerpo, omega, y tau varaibles angulares, e I su
  momento de inercia.
  rutinas: Init=>inicia las variables del cuerpo
           dibujese=>dibuja el cuerpo como un círculo con un radio
*/

class Cuerpo{
private:
  vector3D r,v,F;
  double theta,omega,tau;
  double m,R,I;
  
public:
  //inicia las variables del cuerpo
  void Init(double x0,double y0,double z0,
	    double vx0,double vy0,double vz0,
	    double theta0,double omega0,
	    double m0,double R0);
  //pinta el cuerpo
  void dibujese();

  friend class Colisionador;
};

//Implementaciones de Cuerpo

void Cuerpo::Init(double x0,double y0,double z0,
		  double vx0,double vy0,double vz0,
		  double theta0,double omega0,
		  double m0,double R0){
  r.load(x0,y0,z0); v.load(vx0,vy0,vz0);
  theta=theta0; omega=omega0;
  m=m0; R=R0;
  I=m*R*R/2; //momento de inercia de un disco
}

void Cuerpo::dibujese(){
  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(2*pi*t),"<<r.y()<<"+"<<R<<"*sin(2*pi*t)";
  std::cout<<" , "<<r.x()<<"+"<<R<<"*t*cos("<<theta<<"),"<<r.y()<<"+"<<R<<"*t*sin("<<theta<<")";
}

//-----------------------------------------------------------------------------------

/*Clase Colisionador: accede a las variables privadas de los cuerpos de un sistema
  (vector de cuerpos) y calcula las fuerzas de cada cuerpo. También evoluciona un
  paso en el teimpo usando PEFRL.
 */

class Colisionador{
public:
  //calcula y asiga a cada cuerpo del sistema la fuerza que debe sentir
  void Fuerza_syst(std::vector<Cuerpo> &syst);
  //calcula solo la fuerza entre dos cuerpos del sistema
  void F_entre(Cuerpo &c1,Cuerpo &c2);
  //Ejecuta el paso en el tiempo con una integración PEFRL
  void Paso_syst(std::vector<Cuerpo> &syst,double dt);
};

//Implementaciones de Colisionador

void Colisionador::Fuerza_syst(std::vector<Cuerpo> &syst){
  for(int i=0;i<syst.size();i++){
    syst[i].F.load(0,0,0);
    syst[i].tau=0;
    
  }
  
  //fuerza entre cuerpos
  for(int i=0;i<syst.size();i++){
    for(int j=i+1;j<syst.size();j++){
      F_entre(syst[i],syst[j]);
    }
  }

  //fuerza de gravedad
  for(int i=0;i<syst.size();i++){
    vector3D Fg; Fg.load(0,-syst[i].m*g,0);
    syst[i].F+=Fg;
  }
}

void Colisionador::F_entre(Cuerpo &c1,Cuerpo &c2){
  vector3D r21=c2.r-c1.r;
  double s=c1.R+c2.R-r21.norm(); //interpenetración entre cuerpos
  
  if(s>0){ 
    /*---fuerza de hertz---*/
    double F=K*std::pow(s,1.5)/r21.norm();
    c2.F+=F*r21; c1.F+=-F*r21;
    /*---fuerza plástica---*/
    vector3D v21=c2.v-c1.v;
    double m21=c2.m*c1.m/(c2.m+c1.m);
    vector3D F2=-Gamma*std::sqrt(s)*m21*v21;
    c2.F+=F2; c1.F-=F2;



    /*---EN CASO DE QUE EL PROFE SE EQUIVOQUE---*/
    /*---fuerza de hertz---*/
    //double F=K*std::pow(s,1.5)/r21.norm();
    /*---fuerza plástica---*/
    //vector3D v21=c2.v-c1.v;
    //double m21=c2.m*c1.m/(c2.m+c1.m);
    //F+=-Gamma*std::sqrt(s)*m21*(v21*r21/r21.norm2());
    
    //c2.F+=F*r21; c1.F+=-F*r21;
  }
}

void Colisionador::Paso_syst(std::vector<Cuerpo> &syst,double dt){
  int N=syst.size()-4; //numero de cuerpos sin las paredes
  
  for(int i=0;i<N;i++){
    syst[i].r+=Xi*dt*syst[i].v;
    syst[i].theta+=Xi*dt*syst[i].omega;
  }
  Fuerza_syst(syst);
  for(int i=0;i<N;i++){
    syst[i].v+=coef1*dt*syst[i].F/syst[i].m;
    syst[i].omega+=coef1*dt*syst[i].tau/syst[i].I;
  }
  for(int i=0;i<N;i++){
    syst[i].r+=Chi*dt*syst[i].v;
    syst[i].theta+=Chi*dt*syst[i].omega;
  }
  Fuerza_syst(syst);
  for(int i=0;i<N;i++){
    syst[i].v+=Lambda*dt*syst[i].F/syst[i].m;
    syst[i].omega+=Lambda*dt*syst[i].tau/syst[i].I;
  }
  for(int i=0;i<N;i++){
    syst[i].r+=coef2*dt*syst[i].v;
    syst[i].theta+=coef2*dt*syst[i].omega;
  }
  Fuerza_syst(syst);
  for(int i=0;i<N;i++){
    syst[i].v+=Lambda*dt*syst[i].F/syst[i].m;
    syst[i].omega+=Lambda*dt*syst[i].tau/syst[i].I;
  }
  for(int i=0;i<N;i++){
    syst[i].r+=Chi*dt*syst[i].v;
    syst[i].theta+=Chi*dt*syst[i].omega;
  }
  Fuerza_syst(syst);
  for(int i=0;i<N;i++){
    syst[i].v+=coef1*dt*syst[i].F/syst[i].m;
    syst[i].theta+=coef1*dt*syst[i].tau/syst[i].I;
  }
  for(int i=0;i<N;i++){
    syst[i].r+=Xi*dt*syst[i].v;
    syst[i].theta+=Xi*dt*syst[i].omega;
  }
}

//-----------------------------------------------------------------------------------

void Init_animation(double lx,double ly){
  std::cout<<"set terminal gif animate\n";
  std::cout<<"set output 'balones_fric.gif'\n";
  std::cout<<"set nokey\n";
  std::cout<<"set size ratio -1\n";
  std::cout<<"set xrange[-1:"<<lx<<"+1]\n";
  std::cout<<"set yrange[-1:"<<ly<<"+1]\n";
  std::cout<<"set parametric\n";
  std::cout<<"set trange[0:1]\n";
  std::cout<<"set isosamples 20\n";
  //std::cout<<"set terminal dumb 74 37\n";
}

void Init_frame(double lx,double ly){
  std::cout<<"plot 0,0 , ";
  std::cout<<lx<<"*t,0 , 0,"<<ly<<"*t , "<<lx<<","<<ly<<"*t , "<<lx<<"*t,"<<ly;
}

void End_frame(){
  std::cout<<"\n";
}

int main(){
  /*---variables de simulación---*/
  int Nx=3;               //número de cuerpos en la primera fila x
  int Ny=3;               //número de cuerpos en la primera fila y
  int N=Nx*Ny;            //número de cuerpos totales
  double t=0;             //tiempo inicial
  double tmax=20;         //tiempo final
  double dt=1e-5;         //paso de tiempo
  double m=5;             //masa de los cuerpos
  double R=2;             //radio de los cuerpos
  double lx=(2*Nx+1)*2*R; //lado x de la caja
  double ly=(2*Ny+1)*2*R; //lado y de la caja
  
  /*---variables de animación---*/
  double t1=0;        //timepo auxiliar
  double t_frame=0.1; //tiempo entres frames

  /*---variables de inicialización de la simulación---*/
  std::vector<Cuerpo> syst(N+4);
  Colisionador Newton;

  /*---inicialización de las paredes---*/
  //Init=>(x,y,z,vx,vy,vz,theta,omega,m,R)
  syst[N+0].Init(lx+1000*ly,ly/2.0,0,0,0,0,0,0,100*m,1000*ly); //pared derecha
  syst[N+1].Init(lx/2.0,ly+1000*lx,0,0,0,0,0,0,100*m,1000*lx); //pared arriba
  syst[N+2].Init(-1000*ly,ly/2.0,0,0,0,0,0,0,100*m,1000*ly); //pared izquierda
  syst[N+3].Init(lx/2.0,-1000*lx,0,0,0,0,0,0,100*m,1000*lx); //pared abajo

  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      syst[Ny*i+j].Init((3+4*i)*R,(3+4*j)*R,0,0,0,0,0,1,m,R);
    }
  }
  
  
  Init_animation(lx,ly);
  
  for(t,t1;t<tmax;t+=dt,t1+=dt){
    if(t1>t_frame){
      Init_frame(lx,ly);
      for(int i=0;i<N;i++){
	syst[i].dibujese();
      }
      End_frame();
      t1=0;
    }
    Newton.Paso_syst(syst,dt);
  }
  
  return 0;
}
