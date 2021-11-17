#include <armadillo>
#include <math.h>
#include "matplotlibcpp.h"
#include <iostream>
#include <tuple>

using namespace std;
using namespace arma;
namespace grafica = matplotlibcpp;


//                     FUNCION POTENCIA INVERSA
            
//    Parametros de entrada:
//    A= Matriz a calcular con el metodo de sistemas de Potencia Inversa
//    b = Vector inicial
//    tol= Tolerancia de la aproximacion
//    max_itr= Iteraciones maximas
//    Parametros de salida:
//    Aproximacion del metodo iterativo

tuple<double, mat, int, double>  potencia_inversa(arma::mat A, vec x_0, double tol, int max_itr){
    

    int filas = A.n_rows; //Se obtiene la cantidad de filas de la matriz
    //Se inicializan los valores para las iteraciones del metodo
    mat y_k, x_k = x_0, x_k_n; //Valores de xk y yk
    double c_k; //Variable para las iteraciones del metodo
    double error; //Variable para ir llevando la cuenta del error del metodo

    int interaciones; //Se lleva la cuenta de las iteraciones del metodo

    vector<double> resultados; //Se almacenan los resultados en un vector 
    //para mejor legibilidad de los mismos

    while (interaciones < max_itr)// Iteraraciones
    {
        y_k = arma::solve(A, x_k); //Se resuleve el valor de yk con el vector propio y la matriz
        c_k = norm(yk, "inf");//Se toma la norma del valor de yk para que se guarde en ck del metodo
        x_k_n = (1 / c_k) * y_k; //Se resuelve el valor de la iteracion xk

        error = norm(x_k_n - x_k); //Sr obtiene el erro del metodo
        resultados.push_back(error); //Se agrega el error al metodo
        x_k = x_k_n; //Se agrega el valor de xk 

        //Si el error es menor que la tolerancia
        if (error < tol)
            break; //Se roompe el ciclo y se devuelve el error

        interaciones++;
    }

    return {c_k, x_k.t(), interaciones, error};

    //Se imprime el resultado graficamente

    plt::plot(resultados); //Se usa el vector de resultados para datos
    plt::title("Metodo Potencia Inversa iteraciones vs error"); //Se crea el titulo del grafico
    plt::xlabel("Iteraciones k"); //Eje x con las iteraciones
    plt::ylabel("Error"); //Eje y con el error del metodo
    plt::show(); //Se muestra el resultado
}

int main(int argc, char const *argv[])
{
    //matriz simetrica definida positiva a la cual se le debenalcular el vector propio.
    mat A = {{ 3, -1, 0},
             { -1, 2, -1},
             { 0, -1, 3}};


    vec  x_0= {1,1,1}; //vector inicial aleatorio no nulo.
    x_0= x_0.t(); //Se transpone el valor del vector inicial
    
    potencia_inversa(A, x_0 ,10e-10,9); //Llamado a la funcion
    
    return 0;
}
