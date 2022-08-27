#include<iostream>
#include<cmath>

double cerosporbis(double a, double b){
  double m, fa, fm;
  double ErrMax=1e-7;
  
  fa=std::cyl_bessel_j(0,a);
  while(b-a >= ErrMax){
    m = (b+a)/2; fm=std::cyl_bessel_j(0,m);
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
  std::cout.precision(15);
  
  std::cout<<"ceros de la función de Bessel:"<<"\n";
  std::cout<<cerosporbis(2,3)<<"\n";
  std::cout<<cerosporbis(5,6)<<"\n";
  std::cout<<cerosporbis(8,9)<<"\n";
  std::cout<<cerosporbis(11,12)<<"\n";
  std::cout<<cerosporbis(14,15)<<"\n";
  
  return 0;
}
