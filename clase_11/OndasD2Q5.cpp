#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=128;
const int Ly=128;

const int Q=5; //número de direcciones por celda
const double W0=1.0/3; //peso de la velocidad 0

const double C=0.5;  // C<0.707 cells/click velocidad de las ondas, c no puede
const double C2=C*C; /* ser mayor que 1 ya que permitiría información de
			segundos vecinos e inestabilidades numéricas.*/
const double AUX0=1-3*C2*(1-W0); //parte de la ecuacion de equilibrio

const double tau=0.5; //tiempo de relajación
const double Utau=1.0/tau;  //variables auxiliares
const double UmUtau=1-Utau;

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
  double f_eq(double rho0,double Jx0,double Jy0,int i); //función de equilibrio
  void Start(double rho0,double Jx0,double Jy0); //funcion de inicio del lattice
  void Collision(void);
  void ImposeField(int t);
  void Advection(void);
  void Print(const char *filename);
};

/*---Implementaciones de la clase---*/
LatticeBoltzmann::LatticeBoltzmann(void){
  //Set the weights
  W[0]=W0; W[1]=W[2]=W[3]=W[4]=(1.0-W0)/4;
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;
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
  double sum; int i,n0;
  
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(Usenew){
      sum+=Vx[i]*fnew[n0];
    }
    else{
      sum+=Vx[i]*f[n0];
    }
  }

  return sum;
}

double LatticeBoltzmann::Jy(int ix,int iy,bool Usenew){
  double sum; int i,n0;
  
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(Usenew){
      sum+=Vy[i]*fnew[n0];
    }
    else{
      sum+=Vy[i]*f[n0];
    }
  }

  return sum;
}

//función de equilibrio de LB para ondas.
double LatticeBoltzmann::f_eq(double rho0,double Jx0,double Jy0,int i){
  if(i>0){
    return 3*C2*rho0*W[i]+3*W[i]*(Vx[i]*Jx0+Vy[i]*Jy0);
  }
  else{
    return rho0*AUX0; 
  }
}

/*inicialización de la grilla*/
void LatticeBoltzmann::Start(double rho0,double Jx0,double Jy0){
  int ix,iy,i,n0;
  
  for(ix=0;ix<Lx;ix++){  //para cada celda ix,iy
    for(iy=0;iy<Ly;iy++){
      for(i=0;i<Q;i++){  //para cada dirección
	n0=n(ix,iy,i);
	f[n0]=f_eq(rho0,Jx0,Jy0,i); // se rellena inicialmente cada
      }                             // celda con los valores de f_eq
    }                               // según los valores iniciales de
  }                                 // rho y J.
}
/*paso de colision del LB: se evalúan los valores de rho y J a partir
  de f y se escriben sobre fnew.*/
void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0; double rho0,Jx0,Jy0;
  
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      Jx0=Jx(ix,iy,false);
      Jy0=Jy(ix,iy,false);
      for(i=0;i<Q;i++){
	n0=n(ix,iy,i);
	fnew[n0]=UmUtau*f[n0]+Utau*f_eq(rho0,Jx0,Jy0,i);
      }
    }
  }
}

/*Con esta función se imponen los campos microscópicos en los puntos deseados de acuerdo
  con los campos macro deseados. Por ejemplo en una condición de forzamiento o de frontera.
  Si la condición solo es sobre rho o solo sobre J, la otra variable se calcula usualmente
  como la suma de las fi. Para el caso de esta onda, se impone un forzamiento sinusoidal
  solo sobre la densidad rho y en el centro de la malla.*/
void LatticeBoltzmann::ImposeField(int t){
  int ix,iy,i,n0; double rho0,Jx0,Jy0;
  double lambda=10; //longitud de onda asociada al periodo del forzamiento
  double omega=2*M_PI/lambda*C; //vel angular del forzamiento
  
  //fuente oscilante en la mitad
  ix=Lx/2; iy=Ly/2;
  rho0=10*sin(omega*t); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
  for(i=0;i<Q;i++){
    n0=n(ix,iy,i);
    fnew[n0]=f_eq(rho0,Jx0,Jy0,i);
  }
}

/*función de advección: mueve los valores de las flechas a sus respectivas casillas en f.
  Tiene condiciones periódicas de frontera.*/
void LatticeBoltzmann::Advection(){
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

void LatticeBoltzmann::Print(const char *filename){
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
  LatticeBoltzmann Ondas;
  int t, tmax=100;
  double rho0=0, Jx0=0, Jy0=0;

  /*---inicialización---*/
  Ondas.Start(rho0,Jx0,Jy0);
  /*---simulación/evolución---*/
  for(t=0;t<tmax;t++){
    Ondas.Collision();
    Ondas.ImposeField(t);
    Ondas.Advection();
  }
  /*---print---*/
  Ondas.Print("Ondas.txt");
  
  return 0;
}
