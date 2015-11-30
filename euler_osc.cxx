#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std;

void euler_forw(const int N, double* array, const double dt);
void euler_backw(const int N, double* array2, const double dt);
void output(const int N, double* array, const double dt, double* array2);


int main(){
const double dt = M_PI/10;
const int N = 20*M_PI/dt;
double* array = new double[2*N];
double* array2 = new double[2*N];

array[0] = 1.0;
array[N] = 0.0;

array2[0] = 1.0;
array2[N] = 0.0;

euler_forw(N,array,dt);
euler_backw(N,array2,dt);
output(N,array,dt,array2);

delete[] array;
delete[] array2;
        return 0;
        }

void euler_forw(const int N, double* array, const double dt){
        for(int i = 0; i<N-1; i++){
                array[i+1] = array[i] + dt*array[N+i];
                array[N+1+i] = array[N+i]-dt*array[i];
    }
}

void euler_backw(const int N, double* array2, const double dt){
        const double a = dt*dt + 1;
        for(int i =0;i<N-1; i++){
                array2[i+1] = 1/a * array2[i] + dt*array2[N+i];
                array2[N+i+1] = 1/a * array2[N+i] - dt *array2[i];
        }
}


void output(const int N, double* array, const double dt,double* array2){
        ofstream out("Harmi_gross.txt");

        for(int i = 0; i < N-1; i++){
                out << i * dt << "\t" << array[i] << "\t" << array2[i] << endl;
        }

        out.close();

} 