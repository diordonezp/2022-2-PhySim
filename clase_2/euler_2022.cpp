#include<iostream>
#include<cmath>

double f(double t,double x){
  return x;
}

void Unpasodeeuler(double &t,double &x,double &dt){
  double dx=dt*f(t,x);
  x+=dx;
  t+=dt;
}

int main(){
  double t=0,x=1, dt=0.1;
  
  for(t,x;t<2+dt/2;){
    std::cout<<t<<"\t"<<x<<"\t"<<std::exp(t)<<"\n";
    Unpasodeeuler(t,x,dt);
  }
   
  return 0;
}
