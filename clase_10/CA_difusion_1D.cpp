#include<iostream>
#include<cmath>

const int Lx=1024;
const double p=0.5; //probabilidad de no cmabiar dirección
const int Q=2; //cantidad de direcciones

/*---clase laticce---*/
class LatticeGas{
private:
  double f[Lx][Q],fnew[Lx][Q];
  int V[Q]; //Q=0->derecha
  
public:
  LatticeGas(void);
  void Inicie(int N,double mu,double sigma);
  void GrafiqueRho(void);
  void Colisione();
  void Adveccione(void);
  double rho(int ix);
  double Varianza();
};

/*---implementaciones---*/
LatticeGas::LatticeGas(void){
  V[0]=1; V[1]=-1;
}

void LatticeGas::Inicie(int N,double mu,double sigma){
  for(int ix=0;ix<Lx;ix++){
    double rho=N*1.0/(sigma*std::sqrt(2*M_PI))*std::exp(-0.5*std::pow((ix-mu)/sigma,2.0));
    for(int i=0;i<Q;i++){
      f[ix][i]=rho/Q;
    }
  }
}

double LatticeGas::rho(int ix){
  double suma; int i;
  for(suma=0,i=0;i<Q;i++){
    suma+=f[ix][i];
  }
  return suma;
}

void LatticeGas::GrafiqueRho(void){
  for(int ix=0;ix<Lx;ix++)
    std::cout<<ix<<" "<<rho(ix)<<std::endl;
}

void LatticeGas::Colisione(){
  int j;
  for(int ix=0;ix<Lx;ix++){//para cada celda
    for(int i=0;i<Q;i++){//y en cada direccion
      j=(1+i)%Q;
      fnew[ix][i]=f[ix][i]+(1-p)*(f[ix][j]-f[ix][i]);
    }
  }
}

void LatticeGas::Adveccione(void){
  for(int ix=0;ix<Lx;ix++){
    for(int i=0;i<Q;i++){
      f[(ix+V[i]+Lx)%Lx][i]=fnew[ix][i]; //condiciones de frontera periódicas
    }
  }
}

double LatticeGas::Varianza(){
  int ix;double N,Xprom,sigma2;
  //calcular N
  for(N=0,ix=0;ix<Lx;ix++){
    N+=rho(ix);
  }
  //calcular Xprom
  for(Xprom=0,ix=0;ix<Lx;ix++){
    Xprom+=ix*rho(ix);
  }
  Xprom/=N;
  //calcular Sigma2
  for(sigma2=0,ix=0;ix<Lx;ix++){
    sigma2+=std::pow(ix-Xprom,2.0)*rho(ix);
  }
  sigma2/=N;

  return sigma2;
}

/*---Funcion Principal---*/
int main(){
  LatticeGas Difusion;
  
  int N=400; double mu=Lx/2,sigma=Lx/8;
  int t, tmax=400;
  
  Difusion.Inicie(N,mu,sigma);

  for(t=0;t<tmax;t++){
    std::cout<<t<<" "<<Difusion.Varianza()<<std::endl;
    Difusion.Colisione();
    Difusion.Adveccione();
  }
  
  return 0;
}
