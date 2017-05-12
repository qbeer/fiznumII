#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const int dim = 2; // dimension of vectors
const long double earth_mass = 5.9736e24;
const long double apo = 4.055e8;
const long double vel_apo = 9.64e2;
const long double G = 6.67384e-11;

const double sigma = 10.0;
const double rho = 28.0;
const double beta = 2.66666666667;

double step = 2.4;

void calc_a( long double* param,
             long double t,
             long double* out );
             
void calc_b( long double* param,
             long double t,
             long double* out );

typedef void (*f)( long double*,
                   long double ,
                   long double* );

void rk4( long double* param,
          long double t,
          long double* out,
          f func );

int main()
{

    FILE* fout = fopen("data.dat","w");

    long double t = 0.0;

    //long double in[] = {apo, 0.0, 0.0, vel_apo};
    
    long double in[] = {apo, 0.0, 0, vel_apo};

    long double out[2*dim];

    while( t < 10000000 ){

        rk4(in, t, out, calc_a);

        for( int i=0; i<dim*2; i++ ){

            in[i] = out[i];

            fprintf(fout, "%Lf ", in[i]);

        }

        fprintf(fout, "\n");

        t += step;

    }

    return 0;

}

void calc_a( long double* param,
             long double t,
             long double* out ){

    long double length = 0.0;

    for( int i = 0; i<dim; i++ ){

        length += param[i]*param[i];

    }

    length = sqrt(length);

    long double one_over_length_cubed = 1/pow(length,3);

    for( int i=dim; i<2*dim; i++ ){

        out[i] = -G*earth_mass*param[i-dim]*one_over_length_cubed;

    }
    
    for( int i=0; i<dim; i++ ){
    
        out[i] = param[i+dim] + out[i+dim]*step;
    
    }

}

void calc_b( long double* param,
             long double t,
             long double* out ){
 
    //out[0] = sigma*(param[1]-param[0]); // xdot
    //out[1] = param[0]*(rho-param[2]) - param[1]; // ydot
    //out[2] = param[1]*param[0] - beta*param[2]; // zdot
    
    
    out[3] = sigma*(param[4] - param[3]); // xdotdot
    out[4] = param[3]*(rho-param[2])-param[4] - param[0]*param[5]; // ydotdot
    out[5] = param[3]*param[1] + param[0]*param[4] - beta*param[5]; // zdot
    
    out[0] = param[3] + out[3]*step;
    out[1] = param[4] + out[4]*step;
    out[2] = param[5] + out[5]*step;
      
}

void rk4( long double* param,
          long double t,
          long double* out,
          f func ){

    long double k1[2*dim], k2[2*dim], k3[2*dim], k4[2*dim], temp[2*dim];

    func( param , t, k1 );
    
    for( int i=0; i<2*dim; i++ ){
    
        temp[i] = param[i] + k1[i]*step/2.;
    
    }

    func( temp, t+step/2.0, k2 );

    for( int i=0; i<2*dim; i++){
       
        temp[i] = param[i] + k2[i]*step/2.;

    }

    func( temp, t+step/2.0, k3 );

    for( int i=0; i<dim; i++ ){

        temp[i] = param[i] + k3[i]*step;

    }

    func( temp, t+step, k4 );

    for( int i=0; i<2*dim; i++ ){

        out[i] = param[i] + (step/6.0)*( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] );

    }

}
