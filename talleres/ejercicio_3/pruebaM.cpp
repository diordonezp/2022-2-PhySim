#include<iostream>
#include<cmath>
#include"vector.h"

int main(){
  std::cout.precision(15);
  matrix3D M;
  vector3D v;
  v.load(3,10,2);
  M.load(23.5,4,5.1,
         3.1415,7,90,
         0.0294,23,1);
  
  M.show();
  std::cout<<"times\n";
  v.show();
  std::cout<<"is\n";
  (M*v).show();
  

  return 0;
}