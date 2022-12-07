#include <iostream>
#include <fstream>
#include <cmath>
#include "Random64.h"
using namespace std;

const int Lx=64;
const int Ly=64;
const int Lz=64; //Se agrega la dimensión z
const int Vol=Lx*Ly*Lz;

const int Q=7;
const double W0=1.0/4; //Se cambia el peso a 1/4 en lugar de 1/3
const double D=0.3;    //Coeficiente de absorción

const double C=0.5; // C<0.707 cells/click
const double C2=C*C;
const double AUX0=1-4*C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

//-----Clase LatticeBoltzmann----//
class LatticeBoltzmann{
private:
  double w[Q]; //pesos
  int Vx[Q],Vy[Q],Vz[Q]; //vectores de vel, se agrega en dirección z
  double *f,*fnew; //func. de distribución
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix,int iy,int iz,int i){return ((Ly*ix+iy)*Lz+iz)*Q+i;}; //Se modifica la función indice, así no está mal
  //Se agrega el índice iz para la declaración de todas las funciones 
  double rho(int ix,int iy,int iz,bool UseNew);
  double Jx(int ix,int iy,int iz,bool UseNew);
  double Jy(int ix,int iy,int iz,bool UseNew);
  double Jz(int ix,int iy,int iz,bool UseNew); //Se agrega el campo macroscópico en dirección z
  double feq(double rho0,double Jx0, double Jy0,double Jz0,int i); //Se agrega Jz0
  void Start(double rho0, double Jx0, double Jy0, double Jz0);//Se agrega Jz0
  void Collision(void);
  void ImposeFieldsClassic(int t);
  void ImposeFieldsRand(int N,double *Source,int t);
  void Advection(void);
  void MaxMinRho(double maxrho[Vol], double minrho[Vol]);
  double T_presion(double *maxrho,double *minrho);
  void Print(const char * NameFile);
  void PrintR(const char * NameFile);
  void CalcPoR(const char * NameFile, double maxrho[Vol],double minrho[Vol]);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  //Set the weights
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=w[5]=w[6]=W0/2; //Se agregan los pesos correspondientes a los nuevos vectores del D3Q7 
					       //y se cambian para recuperar ondas acústicas
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=-1; Vx[3]=0; Vx[4]=0;  Vx[5]=0; Vx[6]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=0;  Vy[3]=1; Vy[4]=-1; Vy[5]=0; Vy[6]=0;
  Vz[0]=0;  Vz[1]=0;  Vz[2]=0;  Vz[3]=0; Vz[4]=0;  Vz[5]=1; Vz[6]=-1;
  //Create the dynamic arrays
  int ArraySize=Lx*Ly*Lz*Q; //Se agrega la dimension z
  f=new double [ArraySize];  fnew=new double [ArraySize];
}
LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f;  delete[] fnew;
}

double LatticeBoltzmann::rho(int ix,int iy,int iz,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,iz,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}

double LatticeBoltzmann::Jx(int ix, int iy, int iz, bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,iz,i);
    if(UseNew) sum+=Vx[i]*fnew[n0]; else sum+=Vx[i]*f[n0];
  }
  return sum;
}

double LatticeBoltzmann::Jy(int ix, int iy, int iz, bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,iz,i);
    if(UseNew) sum+=Vy[i]*fnew[n0]; else sum+=Vy[i]*f[n0];
  }
  return sum;
}

double LatticeBoltzmann::Jz(int ix, int iy, int iz, bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,iz,i);
    if(UseNew) sum+=Vz[i]*fnew[n0]; else sum+=Vz[i]*f[n0];
  }
  return sum;
}

double LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,double Jz0, int i){
  if(i>0)
    return 4*w[i]*(C2*rho0+Vx[i]*Jx0+Vy[i]*Jy0+Vz[i]*Jz0);
  else 
    return rho0*AUX0; 
}

void LatticeBoltzmann::Start(double rho0, double Jx0, double Jy0,double Jz0){
  int ix,iy,iz,i,n0;
  for(ix=0;ix<Lx;ix++)//for each cell
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
	for(i=0;i<Q;i++){ //on each direction 
	    n0=n(ix,iy,iz,i);
	  f[n0]=feq(rho0,Jx0,Jy0,Jz0,i);
	}
      }
}

void LatticeBoltzmann::Collision(void){
  int ix,iy,iz,i,n0,n0f; double rho0,Jx0,Jy0,Jz0;
  for(ix=1;ix<Lx-1;ix++){
    for(iy=1;iy<Ly-1;iy++){
      for(iz=1;iz<Lz-1;iz++){
	      rho0=rho(ix,iy,iz,false); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
	      for(i=0;i<Q;i++){
	        n0=n(ix,iy,iz,i);
	        fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,Jz0,i);
	      }
      }
    }
  }

  //condiciones de Bounce back
  //pared izquierda (ix=0)
  for(ix=0,iy=0;iy<Ly;iy++){
    for(iz=0;iz<Lz;iz++){
      for(i=1;i<Q;i++){
	      n0=n(ix,iy,iz,i);
	      n0f=n(ix,iy,iz,(i-2*((int)floor((i+1)/2)-1))%2+2*((int)floor((i+1)/2)-1)+1);
	      fnew[n0f]=D*f[n0];
      }
    }
  }

  //pared derecha (ix=Lx-1)
  for(ix=Lx-1,iy=0;iy<Ly;iy++){
    for(iz=0;iz<Lz;iz++){
      for(i=1;i<Q;i++){
	      n0=n(ix,iy,iz,i);
	      n0f=n(ix,iy,iz,(i-2*((int)floor((i+1)/2)-1))%2+2*((int)floor((i+1)/2)-1)+1);
	      fnew[n0f]=D*f[n0];
      }
    }
  }

  //pared frente (iy=0)
  for(iy=0,ix=0;ix<Lx;ix++){
    for(iz=0;iz<Lz;iz++){
      for(i=1;i<Q;i++){
	      n0=n(ix,iy,iz,i);
	      n0f=n(ix,iy,iz,(i-2*((int)floor((i+1)/2)-1))%2+2*((int)floor((i+1)/2)-1)+1);
	      fnew[n0f]=D*f[n0];
      }
    }
  }

  //pared atrás (iy=Ly-1)
  for(iy=Ly-1,ix=0;ix<Lx;ix++){
    for(iz=0;iz<Lz;iz++){
      for(i=1;i<Q;i++){
	      n0=n(ix,iy,iz,i);
	      n0f=n(ix,iy,iz,(i-2*((int)floor((i+1)/2)-1))%2+2*((int)floor((i+1)/2)-1)+1);
	      fnew[n0f]=D*f[n0];
      }
    }
  }

  //piso (iz=0)
  for(iz=0,ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      for(i=1;i<Q;i++){
	      n0=n(ix,iy,iz,i);
	      n0f=n(ix,iy,iz,(i-2*((int)floor((i+1)/2)-1))%2+2*((int)floor((i+1)/2)-1)+1);
	      fnew[n0f]=D*f[n0];
      }
    }
  }

  //techo (iz=Lz-1)
  for(iz=Lz-1,ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      for(i=1;i<Q;i++){
	      n0=n(ix,iy,iz,i);
	      n0f=n(ix,iy,iz,(i-2*((int)floor((i+1)/2)-1))%2+2*((int)floor((i+1)/2)-1)+1);
	      fnew[n0f]=D*f[n0];
      }
    }
  }
}
void LatticeBoltzmann::ImposeFieldsClassic(int t){
  int i, ix, iy, iz, n0, nf0;
  double lambda,omega,rho0,Jx0,Jy0,Jz0,Aux; lambda=10; omega=2*M_PI/lambda*C;
  ix=Lx/2; iy=Ly/2; iz=Lz/2;
  rho0=10*sin(omega*t); Jx0=Jx(ix,iy,iz,false); Jy0=Jx(ix,iy,iz,false); Jz0=Jx(ix,iy,iz,false);
  for(i=0;i<Q;i++){
    n0=n(ix,iy,iz,i);
    fnew[n0]=feq(rho0,Jx0,Jy0,Jz0,i);
  }
}
void LatticeBoltzmann::ImposeFieldsRand(int N,double *Source,int t){
  int ix,iy,iz,i,n0; double lambda,omega,rho0,Jx0,Jy0,Jz0,fase;
  lambda=10; omega=2*M_PI/lambda*C;

  for(int p=0;p<N;p++){ //Para esta función, se toma la información de las fuentes del array Source, que contiene las
    ix=Source[4*p];     //posiciones y fases de las fuentes y ejecuta el imposefield sobre cada una de ellas.
    iy=Source[4*p+1];
    iz=Source[4*p+2];
    fase=Source[4*p+3];
    rho0=10*sin(omega*t+fase); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
    for(i=0;i<Q;i++){
      n0=n(ix,iy,iz,i);
      fnew[n0]=feq(rho0,Jx0,Jy0,Jz0,i);
    }
  }
}
void LatticeBoltzmann::Advection(void){
  int ix,iy,iz,i,ixnext,iynext,iznext,n0,n0next;
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      for(iz=0;iz<Lz;iz++){
	for(i=0;i<Q;i++){
	  ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly; iznext=(iz+Vz[i]+Lz)%Lz;
	  n0=n(ix,iy,iz,i); n0next=n(ixnext,iynext,iznext,i);
	  f[n0next]=fnew[n0];
	}
      }
    }
  }
}
void LatticeBoltzmann::MaxMinRho(double maxrho[Vol], double minrho[Vol]){
  double rho0; int ix, iy, iz;
  for(ix=0;ix<Lx-1;ix++){
    for(iy=0;iy<Ly-1;iy++){
      for(iz=0;iz<Lz-1;iz++){
        //if(ix!=Lz/2 && iy!=Ly/2 && iz!=Lz/2){
          rho0=rho(ix,iy,iz,true);
          if(rho0>=maxrho[((Ly*ix+iy)*Lz+iz)]){maxrho[((Ly*ix+iy)*Lz+iz)]=rho0;}
          else if(rho0<=minrho[((Ly*ix+iy)*Lz+iz)]){minrho[((Ly*ix+iy)*Lz+iz)]=rho0;}
        //}
      }
    }
  }
}

double LatticeBoltzmann::T_presion(double *maxrho,double *minrho){
  int ix,iy,iz;
  double A,P,suma;
  
  for(suma=0,ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      for(iz=0;iz<Lz;iz++){
	A=(maxrho[((Ly*ix+iy)*Lz+iz)]-minrho[((Ly*ix+iy)*Lz+iz)])/2.0;
	P=A*A;
	suma+=P;
      }
    }
  }

  return suma;
}

void LatticeBoltzmann::Print(const char * NameFile){
  ofstream MyFile(NameFile); double rho0; int ix,iy,iz;
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      for(iz=0;iz<Lz;iz++){
	      rho0=rho(ix,iy,iz,true);
	      MyFile<<ix<<" "<<iy<<" "<<iz<<" "<<rho0<<endl;
	    }
      MyFile<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}
void LatticeBoltzmann::CalcPoR(const char * NameFile, double maxrho[Vol],double minrho[Vol]){
  int ix; double A,P; int iy=Ly/2, iz=Lz/2;
  ofstream MyFile(NameFile);
  //for(int i=0;i<Vol;i++){cout<<maxrho[i]<<""<<minrho[i]<<endl;}
  for(ix=0;ix<Lx;ix++){
    A=(maxrho[((Ly*ix+iy)*Lz+iz)]-minrho[((Ly*ix+iy)*Lz+iz)])/2.0;
    P=A*A;
    //cout<<A<<endl;
    MyFile<<ix<<" "<<P<<endl;
  }
}

/*---Funciones extras---*/
/*función de fuente: Dado un numero de fuentes y 4 características: ix,iy,iz y fase,
  se rellena el array con valores random.*/
void RanSource(int N, double *Source){
  Crandom Ran64(1);
  double Aux;

  for(int n=0;n<N;n++){
    Aux=Ran64.r()*Lx;
    Source[4*n]=(int) Aux;
    Aux=Ran64.r()*Ly;
    Source[4*n+1]=(int) Aux;
    Aux=Ran64.r()*Lz;
    Source[4*n+2]=(int) Aux;
    Aux=Ran64.r()*2*M_PI;
    Source[4*n+3]=Aux;
  }
}

//------Programa Principal------//
int main(void){
  LatticeBoltzmann Ondas;
  int t,tp,tsaturado=650,tmax=1100,N=6;
  double rho0=0, Jx0=0, Jy0=0, Jz0=0;
  double maxrho[Vol], minrho[Vol], Source[4*N];

  RanSource(N,Source);
  
  Ondas.Start(rho0,Jx0,Jy0,Jz0);
  for(int i=0;i<Vol;i++){
    maxrho[i]=0.0;
    minrho[i]=0.0;
  }
  
  for(t=0,tp=0;t<tmax;t++,tp++){
    Ondas.Collision();
    if(t <= tsaturado){Ondas.ImposeFieldsRand(N,Source,t);}
    Ondas.Advection();
    Ondas.MaxMinRho(maxrho,minrho);
    if(tp>25){
      cout<<t<<" "<<Ondas.T_presion(maxrho,minrho)<<"\n";
      tp=0;
      for(int i=0;i<Vol;i++){
	maxrho[i]=0.0;
	minrho[i]=0.0;
      }
    }
  }
  
  //Ondas.CalcPoR("PotenciaX.dat",maxrho,minrho);
  Ondas.Print("Ondas.dat");
  //Ondas.PrintR("Perfil.dat");

  return 0;
}
