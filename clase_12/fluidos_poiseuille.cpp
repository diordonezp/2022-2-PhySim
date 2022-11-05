#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=1;
const int Ly=64;

const int Q=9; //número de direcciones por celda

const double tau=1.2; //tiempo de relajación
const double Utau=1.0/tau;  //variables auxiliares
const double UmUtau=1-Utau;
const double ThreeUmU2tau=3*(1-1/(2*tau));

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
  double Jx(int ix,int iy,bool Usenew,double Fx);
  double Jy(int ix,int iy,bool Usenew,double Fy);
  double Fi(double Ux0,double Uy0,double Fx,double Fy,int i);
  double f_eq(double rho0,double Ux0,double Uy0,int i); //función de equilibrio
  void Start(double rho0,double Ux0,double Uy0); //funcion de inicio del lattice
  void Collision(double gx,double gy);
  void ImposeField(void);
  void Advection(void);
  void Print(const char *filename,double gx,double gy);
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
double LatticeBoltzmann::Jx(int ix,int iy,bool Usenew,double Fx){
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

  return sum+0.5*Fx;
}

double LatticeBoltzmann::Jy(int ix,int iy,bool Usenew,double Fy){
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

  return sum+0.5*Fy;
}

double LatticeBoltzmann::Fi(double Ux0,double Uy0,double Fx,double Fy,int i){
  double VidotF=Vx[i]*Fx+Vy[i]*Fy;
  double UdotF=Ux0*Fx+Uy0*Fy;
  double VidotU=Vx[i]*Ux0+Vy[i]*Uy0;
    
  return W[i]*ThreeUmU2tau*(VidotF-UdotF+3*VidotU*VidotF);
}

//función de equilibrio de LB para ondas.
double LatticeBoltzmann::f_eq(double rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*Vx[i]+Uy0*Vy[i];
  double U2=Ux0*Ux0+Uy0*Uy0;
  
  return W[i]*rho0*(1+3*UdotVi+9.0/2*UdotVi*UdotVi-3.0/2*U2);
}

/*inicialización de la grilla*/
void LatticeBoltzmann::Start(double rho0,double Ux0,double Uy0){
  int ix,iy,i,n0;
  
  for(ix=0;ix<Lx;ix++){  //para cada celda ix,iy
    for(iy=0;iy<Ly;iy++){
      for(i=0;i<Q;i++){  //para cada dirección
	n0=n(ix,iy,i);
	f[n0]=f_eq(rho0,Ux0,Uy0,i); // se rellena inicialmente cada
      }                             // celda con los valores de f_eq
    }                               // según los valores iniciales de
  }                                 // rho y J.
}
/*paso de colision del LB: se evalúan los valores de rho y J a partir
  de f y se escriben sobre fnew.*/
void LatticeBoltzmann::Collision(double gx,double gy){
  int ix,iy,i,n0; double rho0,Ux0,Uy0;
  
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      double Fx=rho0*gx;
      double Fy=rho0*gy;
      Ux0=Jx(ix,iy,false,Fx)/rho0;
      Uy0=Jy(ix,iy,false,Fy)/rho0;   
      
      for(i=0;i<Q;i++){
	n0=n(ix,iy,i);
	fnew[n0]=UmUtau*f[n0]+Utau*f_eq(rho0,Ux0,Uy0,i)+Fi(Ux0,Uy0,Fx,Fy,i);
      }
    }
  }
}

/*Con esta función se imponen los campos microscópicos en los puntos deseados de acuerdo
  con los campos macro deseados. Por ejemplo en una condición de forzamiento o de frontera.
  Si la condición solo es sobre rho o solo sobre J, la otra variable se calcula usualmente
  como la suma de las fi. Para el caso del fluido, se impone que las velocidades de la parte
  superior e inferior de la lámina sean cero.*/
void LatticeBoltzmann::ImposeField(void){
  int ix,iy,i,n0,ixc,iyc,R; double rho0;

  //celda de abajo
  iy=0;
  for(ix=0;ix<Lx;ix++){
    rho0=rho(ix,iy,false);
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=f_eq(rho0,0,0,i);
    }
  }

  //celda de arriba
  iy=Ly-1;
  for(ix=0;ix<Lx;ix++){
    rho0=rho(ix,iy,false);
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=f_eq(rho0,0,0,i);
    }
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

void LatticeBoltzmann::Print(const char *filename,double gx,double gy){
  ofstream Myfile(filename);
  int ix,iy;
  double rho0,Ux0,Uy0;
  
  ix=0;
  for(iy=0;iy<Ly;iy++){
    rho0=rho(ix,iy,true); double Fx=rho0*gx, Fy=rho0*gy;
    Ux0=Jx(ix,iy,true,Fx)/rho0;
    Uy0=Jy(ix,iy,true,Fy)/rho0;
    Myfile<<iy<<"\t"<<Ux0<<"\n";
  }
  
  Myfile.close();
}

//------Programa Principal------//
int main(void){
  LatticeBoltzmann Air;
  int t, tmax=100000;
  double rho0=1.0, g=0.01;

  /*---inicialización---*/
  Air.Start(rho0,0,0);
  /*---simulación/evolución---*/
  for(t=0;t<tmax;t++){
    Air.Collision(g,0);
    Air.ImposeField();
    Air.Advection();
  }
  /*---print---*/
  Air.Print("Poiseuille.txt",g,0);
  
  return 0;
}
