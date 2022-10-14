#include<iostream>
#include<cmath>
#include"Random64.h"

const int Lx=1024;
const double p=0.5; //probabilidad de no cmabiar dirección
const int Q=2; //cantidad de direcciones

/*---clase laticce---*/
class LatticeGas{
private:
  int n[Lx][Q],nnew[Lx][Q];
  int V[Q]; //Q=0->derecha
  
public:
  LatticeGas(void);
  void Borrese(void);
  void Inicie(int N,double mu,double sigma,Crandom &ran64);
  void Show(void);
  void Colisione(Crandom &ran64);
  void Adveccione(void);
  double rho(int ix){return n[ix][0]+n[ix][1];};
  double Varianza();
};

/*---implementaciones---*/
LatticeGas::LatticeGas(void){
  V[0]=1; V[1]=-1;
}

void LatticeGas::Borrese(void){
  for(int ix=0;ix<Lx;ix++){
    for(int i=0;i<Q;i++){
      n[ix][i]=nnew[ix][i]=0;
    }
  }
}

void LatticeGas::Inicie(int N,double mu,double sigma,Crandom &ran64){
  int ix,i;
  while(N>0){
    ix=(int) ran64.gauss(mu,sigma); //Escoger un celda al azar con distr normal
    if(ix<0) ix=0;   //condiciones de bordes
    if(ix>Lx) ix=Lx; //para evitar segmentation fault
    i=(int) Q*ran64.r(); //dirección al azar
    if(n[ix][i]==0){ //si esta vacía
      n[ix][i]=1;
    }
    N--;
  }
}

void LatticeGas::Show(){
  for(int i=0;i<Q;i++){
    for(int ix=0;ix<Lx;ix++){
      std::cout<<n[ix][i];
    }
    std::cout<<endl;
  }
  std::cout<<endl;
}

void LatticeGas::Colisione(Crandom &ran64){
  for(int ix=0;ix<Lx;ix++){//para cada celda
    if(ran64.r()<p){
      nnew[ix][0]=n[ix][1];
      nnew[ix][1]=n[ix][0];
    }
    else{
      nnew[ix][0]=n[ix][0];
      nnew[ix][1]=n[ix][1];
    }
  }
}
void LatticeGas::Adveccione(void){
  for(int ix=0;ix<Lx;ix++){
    for(int i=0;i<Q;i++){
      n[(ix+V[i]+Lx)%Lx][i]=nnew[ix][i]; //condiciones de frontera periódicas
    }
  }
}

double LatticeGas::Varianza(){
  int ix;double N,Xprom,sigma2;
  //calcular N
  for(N=0,ix=0;ix<Lx;ix++){
    N+=rho(ix);
  }
  //calcular Xprom
  for(Xprom=0,ix=0;ix<Lx;ix++){
    Xprom+=ix*rho(ix);
  }
  Xprom/=N;
  //calcular Sigma2
  for(sigma2=0,ix=0;ix<Lx;ix++){
    sigma2+=std::pow(ix-Xprom,2.0)*rho(ix);
  }
  sigma2/=(N-1);

  return sigma2;
}

/*---Funcion Principal---*/
int main(){
  LatticeGas Difusion;
  Crandom ran64(1);

  int N=6; double mu=Lx/2,sigma=Lx/8;
  int t, tmax=400;
  
  Difusion.Borrese();
  Difusion.Inicie(N,mu,sigma,ran64);

  for(t=0;t<tmax;t++){
    std::cout<<t<<" "<<Difusion.Varianza()<<endl;
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
  }
  
  return 0;
}
