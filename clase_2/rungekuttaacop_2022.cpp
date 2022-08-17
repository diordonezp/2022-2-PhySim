#include<iostream>
#include<cmath>

const double omega=3;

double f1(double t,double x1,double x2){
  return -omega*omega*x2;
}

double f2(double t,double x1, double x2){
  return x1;
}

void Unpasoderkuta(double &t,double &x1,double &x2,double &dt){
  double dx11=dt*f1(t,x1,x2);                          double dx12=dt*f2(t,x1,x2);
  double dx21=dt*f1(t+dt/2,x1+dx11/2,x2+dx12/2);       double dx22=dt*f2(t+dt/2,x1+dx11/2,x2+dx12/2);
  double dx31=dt*f1(t+dt/2,x1+dx21/2,x2+dx22/2);       double dx32=dt*f2(t+dt/2,x1+dx12/2,x2+dx22/2);
  double dx41=dt*f1(t+dt,x1+dx31,x2+dx32);             double dx42=dt*f2(t+dt,x1+dx31,x2+dx32);  
  
  x1+=(dx11+2*dx21+2*dx31+dx41)/6;                     x2+=(dx12+2*dx22+2*dx32+dx42)/6;
  t+=dt;
}

int main(){
  double t=0,x1=1, x2=0, dt=0.01;
  
  for(t;t<2+dt/2;){
    std::cout<<t<<"\t"<<x1<<"\t"<<x2<<"\n";
    Unpasoderkuta(t,x1,x2,dt);
  }
   
  return 0;
}
