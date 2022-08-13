#include <iostream>
#include <cmath>

using namespace std;


double f(double x){
  return cos(x);
}

double integralporsimpson(double a, double b,int n){
  n*=2;
  double h=(b-a)/n;
  double sum=0;

  for(int i=0;i<=n;i++){
    double x=a+i*h;
    if(i==0 || i==n){
      sum+=f(x); //si es el primero o el últmo súmelo una vez
    }
    else if(i%2==0){
      sum+=2*f(x); //si es par súmelo 2 veces
    }
    else{
      sum+=4*f(x); //si es impar súmelo 4 veces
    }
  }

  return h/3*sum;
}

int main(){
  double a=0, b=M_PI/2;
  int n=50;
  
  cout<<"La integral es "<<integralporsimpson(a,b,n)<<"\n";

  return 0;
}
