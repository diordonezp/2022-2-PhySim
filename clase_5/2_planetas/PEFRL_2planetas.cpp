#include<iostream>
#include<cmath>
#include"vector.h"

//CONSTANTES GLOBALES
const double G=1.0;
int N=2;//numero de planetas

//constantes de pefrl
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Lambda);


//DECLARACIÓN DE CLASES
class Cuerpo; // se declara la interfaz, luego se declara el interior
class Colisionador;
//--------------Clase Cuerpo--------------------------
class Cuerpo{ 
private: /*datos que no son accesibles, solo son importantes las ordenes*/
  vector3D r,V,F;
  double m,R; 
  
public:  //datos que son accesibles 
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0);
  void BorreFuerza(void); //para indicar que está vacío
  void SumeFuerza(vector3D F0);
  void Arranque(double dt);
  void Muevar(double dt,double theta);
  void MuevaV(double dt,double theta);
  double Getx(void){return r.x();}; //Función Inline
  double Gety(void){return r.y();};
  double Getz(void){return r.z();};

  /*Esto quiere decir que cuerpo va a leer todo dato de colisionador como amigo, y entonces las rutinas de colisionador
    pueden acceder libremente a las variables de los objetos de cuerpo*/
  friend class Colisionador;
};
//aquí acabaría un eventual .h

//----------------------------------------------------
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0){
  /*tengo acceso a todas las variables privada
    porque estoy con Cuerpo::*/
  
  r.load(x0,y0,z0); V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
}

void Cuerpo::BorreFuerza(void){
  F.load(0,0,0);
}

void Cuerpo::SumeFuerza(vector3D F0){
  F+=F0;
}

void Cuerpo::Muevar(double dt,double theta){
  r+=V*(dt*theta);
}

void Cuerpo::MuevaV(double dt,double theta){
  V+=F*(dt*theta/m);
}

/*El colisionador (interaccionador) se ecnarga de ser un ente externo a los cuerpos que calcula sus fuerzas.
  De esta forma solo hay que calcular una vez la fuerza entre los planetas así entrega esa fuerza cuando se
  le pide*/
class Colisionador{
private:
public:
  void Calculefuerzas(Cuerpo *Planeta);
  void Calcularfuerzaentre(Cuerpo &Planeta1, Cuerpo &Planeta2);
};

void Colisionador::Calculefuerzas(Cuerpo *Planeta){
  //borra fuerzas
  for(int i=0;i<N;i++){
    Planeta[i].BorreFuerza();
  }
  //calcula las fuerzas entre todos los pares de planetas
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      Calcularfuerzaentre(Planeta[i],Planeta[j]);
    }
  }
}

void Colisionador::Calcularfuerzaentre(Cuerpo &Planeta1, Cuerpo &Planeta2){
  vector3D r21,F1,n;
  double d21, F;
  r21=Planeta2.r-Planeta1.r; d21=r21.norm(); n=r21/d21;
  F=G*Planeta1.m*Planeta2.m*std::pow(d21,-2.0);
  F1=F*n;
  Planeta1.SumeFuerza(F1); Planeta2.SumeFuerza((-1)*F1);
}


//--------------Funciones Globales--------------------
int main(){
  Cuerpo Planeta[N]; //Ejemplares de la clase Cuerpo (Instance)
  Colisionador newton;
  double t, dt=0.01;
  double m0=1,m1=1,omega, r0=1,V0,V1, T;
  double M=m0+m1;
  double x0=-m1*r0/M;
  double x1=m0*r0/M;
  omega=std::sqrt(G*M*pow(r0,-3)); V0=omega*x0; V1=omega*x1; T=2*M_PI/omega;
  
  //----------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0)
  Planeta[0].Inicie( x0, 0.0, 0.0, 0.0, V0, 0.0, m0, 0.15);
  Planeta[1].Inicie( x1, 0.0, 0.0, 0.0, V1, 0.0, m1, 0.15);
  
  for(t=0;t<2.1*T;t+=dt){
    std::cout<<Planeta[1].Getx()<<" "<<Planeta[1].Gety()<<"\n";
    //mover usando forest_ruth
    for(int i=0;i<N;i++){Planeta[i].Muevar(dt,Zeta);}
    newton.Calculefuerzas(Planeta);
    for(int i=0;i<N;i++){Planeta[i].MuevaV(dt,Coeficiente1);}
    for(int i=0;i<N;i++){Planeta[i].Muevar(dt,Chi);}
    newton.Calculefuerzas(Planeta);
    for(int i=0;i<N;i++){Planeta[i].MuevaV(dt,Lambda);}
    for(int i=0;i<N;i++){Planeta[i].Muevar(dt,Coeficiente2);}
    newton.Calculefuerzas(Planeta);
    for(int i=0;i<N;i++){Planeta[i].MuevaV(dt,Lambda);}
    for(int i=0;i<N;i++){Planeta[i].Muevar(dt,Chi);}
    newton.Calculefuerzas(Planeta);
    for(int i=0;i<N;i++){Planeta[i].MuevaV(dt,Coeficiente1);}
    for(int i=0;i<N;i++){Planeta[i].Muevar(dt,Zeta);}
  }
  
return 0;
}
