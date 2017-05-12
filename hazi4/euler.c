#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define CUBE(x) ((x)*(x)*(x))

const double G = 6.67384E-11;
const double earth = 5.972E24;
const double apo = 405500E3;
const double vel_apo = 964.;
const double M = 0.7349E23;

double gravity( double* r, int size ){

    double length = 0.0;

    for( int i=0; i<size; i++ ){

        length += r[i]*r[i];

    }

    return -G*earth*(1./CUBE(sqrt(length)));

}

double potential( double* r, int size ){

    double length = 0.0;

    for( int i=0; i<size; i++ ){

        length += r[i]*r[i];

    }

    return -G*earth*M*(1./sqrt(length));

}

double kinetic( double* v, int size ){

    double length = 0.0;

    for( int i=0; i<size; i++){

        length += v[i]*v[i];

    }

    return (M/2.)*length;

}

int main()
{

    double T = 50000000.;

    double dt = 1000.;

    double t = 0.0;

    double r[] = { apo, 0. };

    double r_l[] = { apo, 0.};

    double v[] = { 0., vel_apo };

    double v_l[] = { 0., vel_apo };

    FILE* f = fopen("output.dat","w");

    double r_dt[] = { 0.,0. };

    double v_dt[] = { 0.,0. };


    double r_l_dt[] = { 0., 0. };

    double v_l_dt[] = { 0., 0. };

    double v_l_half_dt[] = { 0., 0. };


    double E_tot = 0.0;

    while( t < T ){

        for( int i=0; i<2; i++){

            r_dt[i] = r[i] + dt*v[i];
            v_dt[i] = v[i] + dt*gravity(r, 2)*r[i];

            v_l_half_dt[i] = v_l[i] + 0.5*dt*gravity(r_l, 2)*r_l[i];
            r_l_dt[i] = r_l[i] + v_l_half_dt[i]*dt;

        }

        for( int i=0; i<2; i++){

            v_l_dt[i] = v_l_half_dt[i] + gravity(r_l_dt, 2)*r_l_dt[i]*0.5*dt;

        }

        E_tot = potential(r, 2) + kinetic(v, 2);

        fprintf(f, "%lf %lf %lf %lf %lf %lf\n", r[0], r[1], r_l[0], r_l[1], E_tot, t);

        t += dt;

        for( int i=0; i<2; i++){

            r[i] = r_dt[i];

            v[i] = v_dt[i];

            r_l[i] = r_l_dt[i];

            v_l[i] = v_l_dt[i];

        }

    }

    fclose(f);

    return 0;
}
