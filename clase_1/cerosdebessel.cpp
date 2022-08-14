#include <iostream>
#include <cmath>

using namespace std;

const double ErrMax=1e-7;

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

double cerosporbis(double a, double b,int n,int N){
  double m, fa, fm;
  
  fa=bessel(n,a,N);
  while(b-a >= ErrMax){
    m = (b+a)/2; fm=bessel(n,m,N);
    if(fa*fm>0){
      a=m; fa=fm;
    }
    else{
      b=m;
    }
  }

  return (a+b)/2;
}

int main(){
  int n=0;
  int N=50;
  
    cout<<"la función de Bessel "<<n<<" tiene un cero en x="<<cerosporbis(2,4,n,N)<<"\n";
 
  return 0;
}
