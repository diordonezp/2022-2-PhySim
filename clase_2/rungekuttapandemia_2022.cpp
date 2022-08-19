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
  /*else if(i==3){
    return 3; //función f3(t,x)
    }*/
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

int main(){
  //cond. iniciales y parámetros
  double
    x0=100,    //esta variable es s: los suceptibles
    x1=1,    //esta variable es i: los infectados
    x2=0,    //esta variable es r: los recuperados/muertos
    dt=0.01, //paso de tiempo
    t=0,     //tiempo inicial
    T=10,    //tiempo final
    beta=0.1,  //prob de infección entre s y i
    gamma=1; //prob de sanarse
  std::vector<double> x={x0,x1,x2};

  for(t;t<T;){
    for(int j=0;j<x.size();j++){
      std::cout<<t<<"\t";
      for(int i=0;i<x.size();i++){
	std::cout<<x[i]<<"\t";
      }
      std::cout<<"\n";
      Unpasoderkutta(t,dt,x,beta,gamma);
    }
  }
  
  return 0;
}
