#include <iostream>
#include <cmath>

using namespace std;

double f(double x,double t,int n){
  return cos(n*t-x*sin(t));
}

double bessel(int n,double x,int N){
  N*=2;
  double h=M_PI/N;
  double sum=0;

  for(int i=0;i<=N;i++){
    double t=i*h;
    if(i==0 || i==N){
      sum+=f(x,t,n); //si es el primero o el últmo súmelo una vez
    }
    else if(i%2==0){
      sum+=2*f(x,t,n); //si es par súmelo 2 veces
    }
    else{
      sum+=4*f(x,t,n); //si es impar súmelo 4 veces
    }
  }

  return 1/M_PI*h/3*sum;
}

int main(){
  int N=50;
  int n=0;

  for(double x=0.0;x<=10.0;x+=0.01){
    cout<<x<<"\t"<<bessel(n,x,N)<<"\n";
  }


  return 0;
}
