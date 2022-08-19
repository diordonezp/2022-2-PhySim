#include<iostream>
#include<cmath>

const double g=9.8;

//estamos declaración la clase, para que el compilador sepa que hay una clase
class cuerpo;

/*clase: identidad que tiene datos, internos a los que no se puede acceder,
  y ordenes que se ejecutan para cambiar los datos, que si son públicas*/
class cuerpo{
private: /*acá van los datos privados, solo el objeto accede a ellos. si se cambia una parte del programa,
	   el resto del programa no se tiene que enterar*/
  double x,y,vx,vy,Fx,Fy,m,R;

public: /*si hay datos guardados en publicos, el usuario puede acceder a ellos.
	  Cuando uno cambia el sistema, se tiene que informar a todos. */
  void Inicie(double x0,double y0,double vx0,double vy0,double m0,double R0); //solo se declara: el pc sabe que funciones
  void Fuerza(void);
  void Muevase(double dt);
  double getx(void){return x;}; //esto es una funcion inline declarada ahí mismo por lo sencilla
  double gety(void){return y;};
};

/*acá defino la función. cuerpo:: indica que esa es una función de la clase cuerpo. en esta función se puede meter
  las variables de las que depende (x0,y0,etc) y las privadas de la clase(x,y,etc).*/
void cuerpo::Inicie(double x0,double y0,double vx0,double vy0,double m0,double R0){
  x=x0; y=y0; vx=vx0; vy=vy0; m=m0; R=R0;
}

/*acá definimos la ecuacion de la fuerza*/
void cuerpo::Fuerza(void){
  Fx=0;
  Fy=-m*g;
}

/*algoritmo del paso en dt*/
void cuerpo::Muevase(double dt){
  x+=vx*dt;     y+=vy*dt;    //cambio de las posiciones
  vx+=Fx/m*dt;  vy+=Fy/m*dt; //cambio de las velocidades
}


int main(){
  cuerpo balon; /*cuerpo es un nuevo tipo de datos. cuando definimos balon,
		  se esta creando un espacio de memoria con todos los atributos de cuerpo.
		  balon es lo que se llama un ejemplar (instance en ingles).*/
  
  double t;
  double dt=0.1;
  
  /*----------(x0,y0,vx0,vy0,m0,R0).con balon. accedo a las ordener publicas. EJ: balon.x da error*/
  balon.Inicie(0,0,5,10,1,0.15);
  
  for(t=0;t<3;t+=dt){
    std::cout<<t<<"  "<<balon.getx()<<"  "<<balon.gety()<<std::endl;
    balon.Fuerza();
    balon.Muevase(dt);
  }
  return 0;
}
