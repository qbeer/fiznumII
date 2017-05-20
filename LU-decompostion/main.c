#include <stdio.h>
#include <stdlib.h>

#define N 3

/*

    Probably working properly.

*/

int main()
{
    double** A;

    A = (double**)malloc(sizeof(double*)*N);

    for(int i=0;i<N;i++){

        A[i] = (double*)malloc(sizeof(double)*N);

    }

    A[0][0] = 1.;
    A[0][1] = 2.;
    A[0][2] = 3.;
    A[1][0] = 12.;
    A[1][1] = 32.;
    A[1][2] = 19.;
    A[2][0] = 12.;
    A[2][1] = 6.;
    A[2][2] = 1.;

    double** L;
    double** U;

    L = (double**)malloc(sizeof(double*)*N);
    U = (double**)malloc(sizeof(double*)*N);


    for(int i=0; i<N; i++){

        L[i] = (double*)malloc(sizeof(double)*(i+1));
        U[i] = (double*)malloc(sizeof(double)*(N-i));

    }
    
    for(int i=0;i<N;i++){

        L[i][i] = 1.;

    }

    // Crout's algorithm

    for(int j=0;j<N;j++){

        for(int i=0;i<=j;i++){

            double temp = 0.;

            for(int k=0;k<i;k++){

                temp += L[i][k]*U[k][j];
            }

            U[i][j] = A[i][j] - temp;

        }

        for(int i=j+1;i<N;i++){

            double temp = 0.;

            for(int k=0;k<j;k++){

                temp += L[i][k]*U[k][j];

            }

            L[i][j] = (1./U[j][j])*(A[i][j] - temp );

        }

    }

    printf("\nA matrix is: \n");

    for(int i=0;i<N;i++){

        for(int j=0;j<N;j++){

            printf("%lf ", A[i][j]);

        }

        printf("\n");

    }

    printf("\nL matrix is: \n");

    for(int i=0;i<N;i++){

        for(int j=0;j<=i;j++){

            printf("%lf ", L[i][j]);

        }

        printf("\n");

    }

    printf("\nU matrix is: \n");

    for(int i=0;i<N;i++){

        for(int j=i;j<N;j++){

            printf("%lf ", U[i][j]);

        }

        printf("\n");

    }
    
    printf("\nLU is: \n");
    
    for(int i=0;i<N;i++){
        
        for(int j=0;j<N;j++){
            
            double sol = 0.0;
            
            for(int k=0;k<=i;k++){
                
                sol += L[i][k]*U[k][j]; 
                
            }
            
            printf("%lf ", sol);
            
        }
        
        printf("\n");
        
    }
    
    printf("\ndet(A) is: ");
    
    double temp = 1.;
    
    for(int i=0;i<N;i++){
        
        temp *= U[i][i];
        
    }
    
    printf("%lf\n", temp);

    return 0;
}
