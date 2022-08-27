/*programa de solución de sist. de N ecuaciones diferenciales lineales. de la forma
  dx0/dt=f0(t,x)
  dx1/dt=f1(t,x)   x=(x0,x1,x2,...)
  dx2/dt=f2(t,x)   f=(f0(t,x),f1(t,x),f2(t,x),...)
  .
  .
  .
  Con el algoritmo rungekkuta:
  se tiene el sist dx/dt=f(t,x)
  se soluciona con pasos de tiempo dt en los que el vector x evoluciona a
  x+dx
*/
#include<iostream>
#include<cmath>
#include<vector>

//función f: dado el seleccionador i, f retorna fi la i-esima función del vector
double f(double r,std::vector<double> x,double lambda,int i){
  if(i==0){
    return x[1]; //función f0(t,x)
  }
  else if(i==1){
    return -x[1]/r-lambda*lambda*x[0]; //función f1(t,x)
  }

  return 0;
}


/*función del paso de rungekutta: dado un vector x inicial, un tiempo inicial
  y el paso de tiempo dt*/
void Unpasoderkutta(double &t,double dt,std::vector<double> &x,double lambda){
  std::vector<double> dx1(x.size()),dx2(x.size()),dx3(x.size()),dx4(x.size()),aux(x.size());

  for(int i=0;i<x.size();i++){
    dx1[i]=dt*f(t,x,lambda,i); //dx1 para todas las variables
  }
  for(int i=0;i<x.size();i++){
    aux[i]=x[i]+dx1[i]/2; //aux guarda x+dx1/2 para todas las variables, para obtener dx2
  }
  for(int i=0;i<x.size();i++){
    dx2[i]=dt*f(t+dt/2,aux,lambda,i); //dx2 para todas las variables
  }
  for(int i=0;i<x.size();i++){
    aux[i]=x[i]+dx2[i]/2; //aux guarda x+dx2/2 para todas las variables, para obtener dx2
  }
  for(int i=0;i<x.size();i++){
    dx3[i]=dt*f(t+dt/2,aux,lambda,i); //dx3 para todas las variables
  }
  for(int i=0;i<x.size();i++){
    aux[i]=x[i]+dx3[i]; //aux guarda x+dx3 para todas las variables, para obtener dx2
  }
  for(int i=0;i<x.size();i++){
    dx4[i]=dt*f(t+dt,aux,lambda,i); //dx4 para todas las variables
  }

  for(int i=0;i<x.size();i++){
    x[i]+=(dx1[i]+2*(dx2[i]+dx3[i])+dx4[i])/6; //actualización del xi para todas las variables 
  }
  t+=dt; //actualización del dt
}

/*ahora esta es una función a la que le entra r y lambda, retorna el valor de R(r) para
  el lambda deseado*/
double R(double rf,double lambda){
  //cond. iniciales y parámetros
  double
    x0=1,      //x0 es R
    x1=0,      //x1 es la derivada dR/dr de la ec. de bessel
    dr=0.0001, //paso de radio
    r=0.01;    //radio inicial
  std::vector<double> x={x0,x1};
  int N=(rf-r)/dr;

  
  for(int i=0;i<=N;i++){
    Unpasoderkutta(r,dr,x,lambda);
  }
  
  return x[0];
}

int main(){
  double lambda=0.1,dlambda=0.01;

  for(lambda;lambda<=15.0;lambda+=dlambda){
    std::cout<<lambda<<"\t"<<R(1,lambda)<<"\n";
  }
  return 0;
}
