#include<iostream>
#include<cmath>

double f(double t,double x){
  return x;
}

void Unpasoderkuta(double &t,double &x,double &dt){
  double dx1=dt*f(t,x);
  double dx2=dt*f(t+dt/2,x+dx1/2);
  double dx3=dt*f(t+dt/2,x+dx2/2);
  double dx4=dt*f(t+dt,x+dx3);
  
  x+=(dx1+2*dx2+2*dx3+dx4)/6;
  t+=dt;
}

int main(){
  double t=0,x=1, dt=0.01;
  
  for(t,x;t<2+dt/2;){
    std::cout<<t<<"\t"<<x<<"\t"<<std::exp(t)<<"\n";
    Unpasoderkuta(t,x,dt);
  }
   
  return 0;
}
