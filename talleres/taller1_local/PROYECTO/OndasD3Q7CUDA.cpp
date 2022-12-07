//LB D2Q5 para Ondas en CUDA
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define Lx 64
#define Ly 64
#define Lz 64
#define N 32 //Threads per Block
const int NCells=Lx*Ly*Lz;
const int M=(NCells+N-1)/N; //Blocks per Grid
#define Q 7
const int ArraySize=Lx*Ly*Lz*Q;

const float W0=1.0/4;

const float C=0.5; // C<0.707 cells/click
const float C2=C*C;
const float AUX0=1-4*C2*(1-W0);

const float tau=0.5;
const float Utau=1.0/tau;
const float UmUtau=1-Utau;

//------------------------------------------------------
//------------ PROGRAMMING ON THE DEVICE ----------------
//---------------Constants (Symbols)----------------
__constant__ float d_w[7];
__constant__ int d_Vx[7];
__constant__ int d_Vy[7];
__constant__ int d_Vz[7];
__constant__ float d_C[3];   // d_C[0]=C,  d_C[1]=C2,  d_C[2]=AUX0, 
__constant__ float d_tau[3]; // d_tau[0]=tau,  d_tau[1]=Utau,  d_tau[2]=UmUtau, 

//----------Functions called by the device itself (__device__)
//data index
__device__ int d_n(int ix,int iy,int iz,int i){
  return ((ix*Ly+iy)*Lz+iz)*Q+i;
}
//macroscopic fields
__device__ float d_rho(int ix,int iy,int iz,float *d_f){
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=d_n(ix,iy,iz,i); sum+=d_f[n0];
  }
  return sum;
}
__device__ float d_Jx(int ix,int iy,int iz,float *d_f){
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=d_n(ix,iy,iz,i); sum+=d_Vx[i]*d_f[n0];
  }
  return sum;
}  
__device__ float d_Jy(int ix,int iy,int iz,float *d_f){
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=d_n(ix,iy,iz,i); sum+=d_Vy[i]*d_f[n0];
  }
  return sum;
}
__device__ float d_Jz(int ix,int iy,int iz,float *d_f){
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=d_n(ix,iy,iz,i); sum+=d_Vz[i]*d_f[n0];
  }
  return sum;
}
//equilibrium functions
__device__ float d_feq(float rho0,float Jx0,float Jy0,float Jz0,int i){
  return 4*d_w[i]*(d_C[1]*rho0+d_Vx[i]*Jx0+d_Vy[i]*Jy0+d_Vz[i]*Jz0);
}  
__device__ float d_f0eq(float rho0){
  return rho0*d_C[2];
}  
//---------------------KERNELS----------------------------
__global__ void d_Collision(float *d_f,float *d_fnew){
  //Define internal registers
  int icell,ix,iy,iz,i,n0;  float rho0,Jx0,Jy0,Jz0;
  //Find which thread an which cell should I work
  icell=blockIdx.x*blockDim.x+threadIdx.x;
  ix=icell/(Ly*Lz); iy=(icell/Lz)%Ly; iz=icell%Lz;
  //Compute the macroscopic fields
  rho0=d_rho(ix,iy,iz,d_f); //rho
  Jx0=d_Jx(ix,iy,iz,d_f);   //Jx0
  Jy0=d_Jy(ix,iy,iz,d_f);   //Jy0
  Jz0=d_Jz(ix,iy,iz,d_f);   //Jz0
  //Collide and compute fnew
  n0=d_n(ix,iy,iz,0); d_fnew[n0]=d_tau[2]*d_f[n0]+d_tau[1]*d_f0eq(rho0);
  for(i=1;i<Q;i++){ //on each direction
    n0=d_n(ix,iy,iz,i); d_fnew[n0]=d_tau[2]*d_f[n0]+d_tau[1]*d_feq(rho0,Jx0,Jy0,Jz0,i);
  }
}
__global__ void d_ImposeFields(float *d_f,float *d_fnew,float RhoSource){
  //Define internal registers
  int ix,iy,iz,i,n0;  float rho0,Jx0,Jy0,Jz0;
  //There is only one thread and for one cell in the center
  ix=Lx/2; iy=Ly/2; iz=Lz/2; 
  //Compute the macroscopic fields
  rho0=RhoSource; //rho
  Jx0=d_Jx(ix,iy,iz,d_f);   //Jx0
  Jy0=d_Jy(ix,iy,iz,d_f);   //Jy0
  Jz0=d_Jz(ix,iy,iz,d_f);   //Jz0
  //Collide and compute fnew
  n0=d_n(ix,iy,iz,0); d_fnew[n0]=d_f0eq(rho0);
  for(i=1;i<Q;i++){ //on each direction
    n0=d_n(ix,iy,iz,i); d_fnew[n0]=d_feq(rho0,Jx0,Jy0,Jz0,i);
  }
}
__global__ void d_Advection(float *d_f,float *d_fnew){
  //Define internal registers
  int icell,ix,iy,iz,i,ixnext,iynext,iznext,n0,n0next;
  //Find which thread an which cell should I work
  icell=blockIdx.x*blockDim.x+threadIdx.x;
  ix=icell/(Ly*Lz); iy=(icell/Lz)%Ly; iz=icell%Lz;
  //Move the contents to the neighboring cells
  for(i=0;i<Q;i++){ //on each direction
    ixnext=(ix+d_Vx[i]+Lx)%Lx; iynext=(iy+d_Vy[i]+Ly)%Ly; iznext=(iz+d_Vz[i]+Lz)%Lz;//periodic boundaries
    n0=d_n(ix,iy,iz,i); n0next=d_n(ixnext,iynext,iznext,i);
    d_f[n0next]=d_fnew[n0]; 
  }
}
//------------------------------------------------------
//------------ PROGRAMMING ON THE HOST ----------------
//--------------------- Clase LatticeBoltzmann ------------
class LatticeBoltzmann{
private:
  float h_C[3];   // h_C[0]=C,  h_C[1]=C2,  h_C[2]=AUX, 
  float h_tau[3]; // h_tau[0]=tau,  h_tau[1]=Utau,  h_tau[2]=UmUtau, 
  float h_w[Q];      //Weights 
  int h_Vx[Q],h_Vy[Q],h_Vz[Q];  //Velocity vectors
  float *h_f, *h_fnew;  float *d_f,*d_fnew; //Distribution Functions
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int h_n(int ix,int iy,int iz,int i){return ((ix*Ly+iy)*Lz+iz)*Q+i;};
  float h_rho(int ix,int iy,int iz);
  float h_Jx(int ix,int iy,int iz);
  float h_Jy(int ix,int iy,int iz);
  float h_Jz(int ix,int iy,int iz);
  float h_feq(float rho0,float Jx0,float Jy0,float Jz0,int i);
  void Start(float rho0,float Jx0,float Jy0,float Jz0);
  void Collision(void);
  void ImposeFields(int t);
  void Advection(void);
  void Print(const char * NameFile);

};
LatticeBoltzmann::LatticeBoltzmann(void){
  //CONSTANTS(d_Symbols)
  //---Charge constantes on the Host-----------------
  //running constants
  h_C[0]=C;  h_C[1]=C2;  h_C[2]=AUX0;
  h_tau[0]=tau;  h_tau[1]=Utau;  h_tau[2]=UmUtau;
  //Set the weights
  h_w[0]=W0; h_w[1]=h_w[2]=h_w[3]=h_w[4]=h_w[5]=h_w[6]=W0/2;
  //Set the velocity vectors
  h_Vx[0]=0;  h_Vx[1]=1;  h_Vx[2]=-1; h_Vx[3]=0; h_Vx[4]=0;  h_Vx[5]=0; h_Vx[6]=0;
  h_Vy[0]=0;  h_Vy[1]=0;  h_Vy[2]=0;  h_Vy[3]=1; h_Vy[4]=-1; h_Vy[5]=0; h_Vy[6]=0;
  h_Vz[0]=0;  h_Vz[1]=0;  h_Vz[2]=0;  h_Vz[3]=0; h_Vz[4]=0;  h_Vz[5]=1; h_Vz[6]=-1;
  //------Send to the Device-----------------
  cudaMemcpyToSymbol(d_w,h_w,Q*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vx,h_Vx,Q*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vy,h_Vy,Q*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vz,h_Vz,Q*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_C,h_C,3*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_tau,h_tau,3*sizeof(float),0,cudaMemcpyHostToDevice);
  
  //DISTRIBUTION FUNCTIONS
  //Build the dynamic matrices on the host
  h_f=new float [ArraySize];  h_fnew=new float [ArraySize];
  //Build the dynamic matrices on the device
  cudaMalloc((void**) &d_f,ArraySize*sizeof(float));
  cudaMalloc((void**) &d_fnew,ArraySize*sizeof(float));
}
LatticeBoltzmann::~LatticeBoltzmann(void){
  delete[] h_f;  delete[] h_fnew;
  cudaFree(d_f);  cudaFree(d_fnew);
}
float LatticeBoltzmann::h_rho(int ix,int iy,int iz){
  //Note: Please import data from device before running
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=h_n(ix,iy,iz,i); sum+=h_fnew[n0];
  }
  return sum;
}  
float LatticeBoltzmann::h_Jx(int ix,int iy,int iz){
  //Note: Please import data from device before running
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=h_n(ix,iy,iz,i); sum+=h_Vx[i]*h_fnew[n0];
  }
  return sum;
}  
float LatticeBoltzmann::h_Jy(int ix,int iy,int iz){
  //Note: Please import data from device before running
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=h_n(ix,iy,iz,i); sum+=h_Vy[i]*h_fnew[n0];
  }
  return sum;
}
float LatticeBoltzmann::h_Jz(int ix,int iy,int iz){
  //Note: Please import data from device before running
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=h_n(ix,iy,iz,i); sum+=h_Vz[i]*h_fnew[n0];
  }
  return sum;
}
float  LatticeBoltzmann::h_feq(float rho0,float Jx0,float Jy0,float Jz0,int i){
  if(i>0)
    return 4*h_w[i]*(h_C[1]*rho0+h_Vx[i]*Jx0+h_Vy[i]*Jy0+h_Vz[i]*Jz0);
  else
    return rho0*AUX0;
}  
void LatticeBoltzmann::Start(float rho0,float Jx0,float Jy0,float Jz0){
  //Charge on the Host
  int ix,iy,iz,i,n0;
  for(ix=0;ix<Lx;ix++){ //for each cell
    for(iy=0;iy<Ly;iy++){
      for(iz=0;iz<Lz;iz++){
        for(i=0;i<Q;i++){ //on each direction
          n0=h_n(ix,iy,iz,i); h_f[n0]=h_feq(rho0,Jx0,Jy0,Jz0,i);
        }
      }
    }
  }  
  //Send to the Device
  cudaMemcpy(d_f,h_f,ArraySize*sizeof(float),cudaMemcpyHostToDevice);
}  
void LatticeBoltzmann::Collision(void){
  //Do everything on the Device
  dim3 ThreadsPerBlock(N,1,1);
  dim3 BlocksPerGrid(M,1,1);
  d_Collision<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f,d_fnew);
}
void LatticeBoltzmann::ImposeFields(int t){
  float lambda=10, omega=2*M_PI/lambda*C;
  float RhoSource=10*sin(omega*t);
  dim3 ThreadsPerBlock(1,1,1); //A single thread (in this case)
  dim3 BlocksPerGrid(1,1,1);
  d_ImposeFields<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f,d_fnew,RhoSource);
}
void LatticeBoltzmann::Advection(void){
  //Do everything on the Device
  dim3 ThreadsPerBlock(N,1,1);
  dim3 BlocksPerGrid(M,1,1);
  d_Advection<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f,d_fnew);
}
void LatticeBoltzmann::Print(const char * NameFile){
  ofstream MyFile(NameFile); float rho0; int ix,iy,iz;
  //Bring back the data from Device to Host
  cudaMemcpy(h_fnew,d_fnew,ArraySize*sizeof(float),cudaMemcpyDeviceToHost);
  //Print for gnuplot splot  
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      for(iz=0;iz<Lz;iz++){
        rho0=h_rho(ix,iy,iz);
        MyFile<<ix<<" "<<iy<<" "<<iz<<" "<<rho0<<endl;
      }
      MyFile<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}

//------------------- Funciones Globales ------------

int main(void){
  LatticeBoltzmann Ondas;
  int t,tmax=60;
  float rho0=0,Jx0=0,Jy0=0,Jz0=0;

  //Start
  Ondas.Start(rho0,Jx0,Jy0,Jz0);
  //Run
  for(t=0;t<tmax;t++){
    Ondas.Collision();
    Ondas.ImposeFields(t);
    Ondas.Advection();
  }
  //Print
  Ondas.Print("Ondas.dat");

  return 0;
}