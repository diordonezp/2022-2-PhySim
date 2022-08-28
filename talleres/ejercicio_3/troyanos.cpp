#include<iostream>
#include<cmath>
#include<vector>
#include"vector.h"

//Constantes globables
double G=1; //Cte gravitacional

//-----------------------------------------------------------------------------------

/*Clase Cuerpo: objeto que guarda las coordenadas r,v,F,m y R (radio) de cada cuerpo
  del sistema
  rutinas: Init=>inicia las variables del cuerpo
           print=>printea en pantalla todas las variables de cuerpo
           Paso_euler=>evoluciona las varibles del cuerpo integrando por euler
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
  //printea las variables r del cuerpo
  void print();
  //evolución temporal del cuerpo
  void Paso_euler(double dt);

  friend class Colisionador;
};

//Implementaciones de Cuerpo
void Cuerpo::Init(double x0,double y0,double z0,
		  double vx0,double vy0,double vz0,
		  double m0,double R0){
  r.load(x0,y0,z0); v.load(vx0,vy0,vz0);
  m=m0; R=R0;
}

void Cuerpo::print(){
  std::cout<<r.x()<<"\t"<<r.y()<<"\t"<<r.z()<<"\t"<<F.x()<<"\t"<<F.y();
}

void Cuerpo::Paso_euler(double dt){
  r+=v*dt; v+=F/m*dt;
}

//-----------------------------------------------------------------------------------

/*Clase Colisionador: accede a las variables privadas de los cuerpos de un sistema
  (vector de cuerpos)y calcula las fuerzas de cada cuerpo.
 */
class Colisionador{
public:
  //calcula y asiga a cada cuerpo del sistema la fuerza que debe sentir
  void Fuerza_syst(std::vector<Cuerpo> &syst);
  //calcula solo la fuerza entre dos cuerpos del sistema
  void F_entre(Cuerpo c1,Cuerpo c2);
};

//Implementaciones de Colisionador
void Colisionador::Fuerza_syst(std::vector<Cuerpo> &syst){
  for(int i=0;i<syst.size();i++){
    syst[i].F=(-G*syst[i].m/std::pow(syst[i].r.norm(),3))*syst[i].r;
  }
}

//-----------------------------------------------------------------------------------



int main(){
  int N=1; //número de cuerpos
  double dt=0.001; //paso de tiempo
  //Condiciones iniciales
  double R=1;
  double m=5;
  double x0=10;
  double y0=0;
  double z0=0;
  double vx0=0;
  double vy0=std::sqrt(G*std::pow(x0,-3))*x0/2;
  double vz0=0;
  std::vector<Cuerpo> syst(N); //creación de sistema=>un vector de cuerpos
  Colisionador G; //el colisionador=>la gravedad
  
  syst[0].Init(x0,y0,z0,vx0,vy0,vz0,m,R);
  for(double t=0;t<4*M_PI*x0/vy0;t+=dt){
    syst[0].print();std::cout<<"\n";
    G.Fuerza_syst(syst);
    syst[0].Paso_euler(dt);
  }
  
  return 0;
}
