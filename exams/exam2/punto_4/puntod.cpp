#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;

const int Lx=256;
const int Ly=64;

const int Q=9; //número de direcciones por celda

const double tau=0.55; //tiempo de relajación
const double Utau=1.0/tau;  //variables auxiliares
const double UmUtau=1-Utau;

const double dp=0.001;
const double eta=10;

//-----Clase LatticeBoltzmann----//
class LatticeBoltzmann{
private:
  double W[Q]; //pesos
  int Vx[Q],Vy[Q]; //vectores de vel
  double *f, *fnew; //func. de distribución (punteros)
public:
  LatticeBoltzmann(void); //constructor
  ~LatticeBoltzmann(void); //destructor
  int n(int ix,int iy,int i){return (ix*Ly+iy)*Q+i;}; //funcion de índice
  double rho(int ix,int iy,bool UseNew); //cálculo de campos macroscópicos
  double Jx(int ix,int iy,bool Usenew);
  double Jy(int ix,int iy,bool Usenew);
  double f_eq(double rho0,double Ux0,double Uy0,int i); //función de equilibrio
  void Start(double mu,double sigma,double Ux0,double Uy0); //funcion de inicio del lattice
  void Collision(void);
  void Advection(void);
  double Varianza(void);
  double Detector(int ix);
  void Print(const string filename);
};

/*---Implementaciones de la clase---*/
LatticeBoltzmann::LatticeBoltzmann(void){
  //Set the weights
  W[0]=4.0/9; W[1]=W[2]=W[3]=W[4]=1.0/9; W[5]=W[6]=W[7]=W[8]=1.0/36;
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;
  
  Vx[5]=1;  Vx[6]=-1; Vx[7]=-1; Vx[8]=1;
  Vy[5]=1;  Vy[6]=1;  Vy[7]=-1; Vy[8]=-1;
  
  //Create the dynamic arrays
  int ArraySize=Lx*Ly*Q;
  f=new double [ArraySize];  fnew=new double [ArraySize];
}
LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f;  delete[] fnew;
}

/*función de densidad: le entra la casilla y retorna rho: la suma de las f's para
  las distintas direcciones de la casilla*/
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double sum=0; int i,n0;
  
  for(i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}

/*funciones de corriente: le entra la casilla y retorna corriente: la suma de
  las f's por las velocidades para las distintas direcciones de la casilla.
  J=sum_i V_i*f_i.*/
double LatticeBoltzmann::Jx(int ix,int iy,bool Usenew){
  return 0.5*rho(ix,iy,Usenew)*dp/eta*iy*(Ly-1-iy);
}

double LatticeBoltzmann::Jy(int ix,int iy,bool Usenew){
  return 0;
}

//función de equilibrio de LB para fluidos.
double LatticeBoltzmann::f_eq(double rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*Vx[i]+Uy0*Vy[i];
  double U2=Ux0*Ux0+Uy0*Uy0;
  
  return W[i]*rho0*(1+3*UdotVi+9.0/2*UdotVi*UdotVi-3.0/2*U2);
}

/*inicialización de la grilla*/
void LatticeBoltzmann::Start(double mu,double sigma,double Ux0,double Uy0){
  int ix,iy,i,n0; double rho0;
  
  for(ix=0;ix<Lx;ix++){  //para cada celda ix,iy
    for(iy=0;iy<Ly;iy++){
      rho0=1/(sigma*sqrt(2*M_PI))*exp(-1/2.0*pow((ix-mu)/sigma,2)); //función inicial de densidad
      for(i=0;i<Q;i++){  //para cada dirección
	n0=n(ix,iy,i);
	f[n0]=f_eq(rho0,Ux0,Uy0,i); // se rellena inicialmente cada
      }                             // celda con los valores de f_eq
    }                               // según los valores iniciales de
  }                                 // rho y J.
}

/*paso de colision del LB: se evalúan los valores de rho y J a partir
  de f y se escriben sobre fnew.*/
void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0; double rho0,Ux0,Uy0;
  
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      Ux0=Jx(ix,iy,false)/rho0;
      Uy0=Jy(ix,iy,false)/rho0;
      for(i=0;i<Q;i++){
	n0=n(ix,iy,i);
	fnew[n0]=UmUtau*f[n0]+Utau*f_eq(rho0,Ux0,Uy0,i);
      }
    }
  }
}

/*función de advección: mueve los valores de las flechas a sus respectivas casillas en f.
  Tiene condiciones periódicas de frontera.*/
void LatticeBoltzmann::Advection(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      for(i=0;i<Q;i++){
	n0=n(ix,iy,i);
	ixnext=(ix+Vx[i]+Lx)%Lx;
	iynext=(iy+Vy[i]+Ly)%Ly;
	n0next=n(ixnext,iynext,i);
	f[n0next]=fnew[n0];
      }
    }
  }
}

double LatticeBoltzmann::Varianza(void){
  int ix,iy,i,n0;
  double N,var,xprom;

  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      N+=rho(ix,iy,true);
      xprom+=rho(ix,iy,true)*ix;
    }
  }
  xprom/=N;

  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      var+=rho(ix,iy,true)*(ix-xprom)*(ix-xprom);
    }
  }
  var/=N;

  return var;
}

double LatticeBoltzmann::Detector(int ix){
  int iy;
  double rho0;

  for(rho0=0,iy=0;iy<Ly;iy++){
    rho0+=rho(ix,iy,true);
  }

  return rho0;
}

void LatticeBoltzmann::Print(const string filename){
  ofstream Myfile(filename);
  int ix,iy;
  double rho0;

  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true);
      Myfile<<ix<<"\t"<<iy<<"\t"<<rho0<<"\n";
    }
    Myfile<<"\n";
  }
  Myfile.close();
}

//------Programa Principal------//
int main(void){
  LatticeBoltzmann fluid;
  int t, tmax=5000;
  double mu=Lx/8.0,sigma=Ly/9.0;

  /*---inicialización---*/
  fluid.Start(mu,sigma,0,0);

  /*---simulación/evolución---*/
  for(t=0;t<=tmax;t++){
    fluid.Collision();
    fluid.Advection();
    cout<<t<<"\t"<<fluid.Detector(150)<<"\n";
  }
  
  fluid.Print("d_Densidad5000.dat");
  
  return 0;
}
