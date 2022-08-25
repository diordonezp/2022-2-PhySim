/*programa de solución de sist. de N ecuaciones diferenciales lineales. de la forma
  dx0/dt=f0(t,x)
  dx1/dt=f1(t,x)   
  dx2/dt=f2(t,x)  
  .                x=(x0,x1,x2,...,xN)
  .                f=(f0(t,x),f1(t,x),f2(t,x),...,fN(t,x))
  .
  dxN/dt=fN(t,x)
  
  Con el algoritmo rungekkuta:
  se tiene el sist dx/dt=f(t,x)
  se soluciona con pasos de tiempo dt en los que el vector x evoluciona a
  x+dx
*/
#include<iostream>
#include<cmath>
#include<vector>

//función f: dado el seleccionador i, f retorna fi la i-esima función del vector
double f(double t,std::vector<double> x,double beta,double gamma,int i){
  if(i==0){
    return -beta*x[0]*x[1]; //función f0(t,x)
  }
  else if(i==1){
    return beta*x[0]*x[1]-gamma*x[1]; //función f1(t,x)
  }
  else if(i==2){
    return gamma*x[1]; //función f2(t,x)
  }
  
  return 0;
}


/*función del paso de rungekutta: dado un vector x inicial, un tiempo inicial
  y el paso de tiempo dt*/
void Unpasoderkutta(double &t,double dt,std::vector<double> &x,double beta, double gamma){
  std::vector<double> dx1(x.size()),dx2(x.size()),dx3(x.size()),dx4(x.size()),aux(x.size());

  for(int i=0;i<x.size();i++){
    dx1[i]=dt*f(t,x,beta,gamma,i); //dx1 para todas las variables
  }
  for(int i=0;i<x.size();i++){
    aux[i]=x[i]+dx1[i]/2; //aux guarda x+dx1/2 para todas las variables, para obtener dx2
  }
  for(int i=0;i<x.size();i++){
    dx2[i]=dt*f(t+dt/2,aux,beta,gamma,i); //dx2 para todas las variables
  }
  for(int i=0;i<x.size();i++){
    aux[i]=x[i]+dx2[i]/2; //aux guarda x+dx2/2 para todas las variables, para obtener dx2
  }
  for(int i=0;i<x.size();i++){
    dx3[i]=dt*f(t+dt/2,aux,beta,gamma,i); //dx3 para todas las variables
  }
  for(int i=0;i<x.size();i++){
    aux[i]=x[i]+dx3[i]; //aux guarda x+dx3 para todas las variables, para obtener dx2
  }
  for(int i=0;i<x.size();i++){
    dx4[i]=dt*f(t+dt,aux,beta,gamma,i); //dx4 para todas las variables
  }

  for(int i=0;i<x.size();i++){
    x[i]+=(dx1[i]+2*(dx2[i]+dx3[i])+dx4[i])/6; //actualización del xi para todas las variables 
  }
  t+=dt; //actualización del dt
}

double susceptibles(double beta){
  //cond. iniciales y parámetros
  double
    x0=0.999,   //esta variable es s: los susceptibles
    x1=0.001,   //esta variable es i: los infectados
    x2=0,       //esta variable es r: los recuperados/muertos
    dt=0.001,   //paso de tiempo
    t=0,        //tiempo inicial
    T=100,      //tiempo final
    gamma=0.08; //prob de sanarse
  std::vector<double> x={x0,x1,x2};

  for(t;t<T;){
    Unpasoderkutta(t,dt,x,beta,gamma);
  }
  
  return x[0];
}

int main(void){
  double beta=0.1;
  double dbeta=0.01;
  
  for(beta;beta<1;beta+=dbeta){
    std::cout<<beta/0.08<<"\t"<<susceptibles(beta)<<"\n";
  }
}
