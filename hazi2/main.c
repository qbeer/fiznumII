#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// SVD FUNCTION
void dgesv_(const int *N, const int *nrhs, double *A, const int *lda, int
	*ipiv, double *b, const int *ldb, int *info);

void dgels_(const char *trans, const int *M, const int *N, const int *nrhs, double *A, const int *lda, double *b, const int *ldb, double *work, const int * lwork, int *info);

// This is just for cosmetics. I don't want the ATTENTION
// text to appear several times.
unsigned int flag = 0;

// I cheched this on WikiPedia, this is approximately the standard
// machine epsilon value for doubleing point numbers.
const double eps = 10e-8;

// It is not efficient on small matrices and vectors but it is more efficient on bigger ones.
unsigned int DEF_SIZE = 20;

// This is how I get the information about row and column count in a vector typed matrix.
typedef struct{
    int length;
    int row;
    int column;
    double* matrix;
}mtrx;

// The matrix and vector read in function. My program can accept matrices as well as vectors for input and
// multiply them if it's possible
void matrix_in(FILE*, mtrx*);
// Tested

// The multiplication function.
void multiply(mtrx* A, mtrx* b);
// Tested

// Normalize the matrix and the solution vector.
void normalize(mtrx*,mtrx*);
// Tested

// Switch rows. The row numbering starts from 0 !
void switchRows(mtrx*,int,int);
// Tested

// Printing the matrix.
void printOut(mtrx*);
// Tested

// Printing the matrix and the solution vector as well.
// This is for debugging purposes only.
void printOutBoth(mtrx*, mtrx*);
// Tested

// Set the diagonals to 1. and normalize the other elements.
void diagonal(mtrx*,mtrx*,mtrx*);
// Tested

// Eliminates the numbers below the diagonal line! After it is set
// to ones.
void eliminate(mtrx*, mtrx*, int);
// Tested

// Eliminate above the diagonal line
void back_eliminate(mtrx*, mtrx*);
// Tested, but it was messed up!

// Reorder the solution of the linear equation
void reOrder(mtrx*, mtrx*);
// Tested

// Basic linear solver
void lin_solve(mtrx*,mtrx*,mtrx*);
//Tested

// Get the transponse of a matrix.
void transponse(mtrx*);
// Tested

// Get the inverse of the provided matrix.
void get_inverse(mtrx*,double**);
// Tested

// Calculate the inverse matrix with lapack
void lapack_inverse(mtrx* BASE, double**);
// Tested

// Compere the inverse matrices.
void compare(double**,double**,int);


int main(int argc, char* argv[])
{

    clock_t t0,t1,t2,t3;

    if(argc>2){
        printf("There should be an argument. \n");
        exit(1);
    }

    // Load matrix from file.
    mtrx* A = (mtrx*)malloc(sizeof(mtrx));

    FILE* fp = fopen(argv[1],"r");
    if(!fp){
        printf("An error occured loading %s \n", argv[1]);
        exit(1);
    }

    A->matrix = (double*)malloc(DEF_SIZE*sizeof(double));
    matrix_in(fp,A);
    // Matrix loaded.

    fclose(fp);

    //printOut(A);

    // Check wheter the matrix is square or not.
    if(A->column!=A->row){
        printf("The loaded matrix is not a square-matrix. Exiting...\n");
        exit(1);
    }

    printf("\n The inverse of the matrix calculated by me: \n\n");

    double **a = (double**)malloc(sizeof(double*)*A->row);
    for(int i=0;i<A->row;i++){
        a[i] = (double*)malloc(sizeof(double)*A->column);
    }
    t0 = clock();
    get_inverse(A,a);
    t1 = clock();
    printf("\nMy calculations took: %ld cpu clocks\n",(t1-t0));

    /*
    // Solution vector. Need to modify it in the code currently.
    mtrx* sol = (mtrx*)malloc(sizeof(mtrx));
    sol->matrix = (double*)malloc(sizeof(double)*A->row);
    sol->row = A->row;
    sol->column = 1;
    sol->length = A->row;
    // Order 'vector' to store row switches.
    mtrx* order = (mtrx*)malloc(sizeof(mtrx));
    order->matrix = (double*)malloc(sizeof(double)*A->row);
    order->row = A->row;
    order->column = 1;
    order->length = A->row;
    // For testing.
    for(int i=0;i<sol->row;i++){
        sol->matrix[i] = 0;
        order->matrix[i] = i;
    }
    sol->matrix[0] = 1.;
    */

    //printOut(sol);

    //lin_solve(A,sol,order);

    //printOut(sol);

    //printOutBoth(A,sol);

    // Reload the matrix! Because A is messed up at the moment.
    mtrx* B;
    B = (mtrx*)malloc(sizeof(mtrx));
    DEF_SIZE = 20;
    B->matrix = (double*)malloc(sizeof(double)*DEF_SIZE);
    fp = fopen(argv[1],"r");
    matrix_in(fp,B);

    printf("\n\n And the lapack inverse is:\n\n");

    double **b = (double**)malloc(sizeof(double*)*B->row);
    for(int i=0;i<B->row;i++){
        b[i] = (double*)malloc(sizeof(double)*B->column);
    }
    t2 = clock();
    lapack_inverse(B,b);
    t3 = clock();
    printf("\nLapack's calculations took: %ld cpu clocks \n " , (t3-t2));

    compare(a,b,B->row);

    printf("\n");

    return 0;

}

void matrix_in(FILE* fp, mtrx* matrix)
{

    char* line = (char*)malloc(10000);
    size_t chars;

    matrix->row = 0;
    matrix->column = 0;
    matrix->length = 0;

    // To handle spaces before numbers on the beginning of a row I needed to add this.
    int col = 0;

    while(getline(&line, &chars,fp) != -1){
        if(strspn(line,"-0123456789. ")){
        char delims[] = " \t";
        char* word = NULL;
        word = strtok(line,delims);
        if(strspn(word,"-0123456789.")){
        matrix->matrix[matrix->length++] = atof(word);
        matrix->column++;
        }
        while((word=strtok(NULL,delims))!=NULL){
            if(strspn(word,"-0123456789.")){
                if(matrix->length<DEF_SIZE){
                    matrix->matrix[matrix->length++] = atof(word);
                    matrix->column++;
                }else{
                    DEF_SIZE*=2;
                    matrix->matrix = (double*)realloc(matrix->matrix,sizeof(double)*DEF_SIZE);
                    matrix->matrix[matrix->length++] = atof(word);
                    matrix->column++;
                }
            }
        }
        if(matrix->column!=0){
            matrix->row++;
            col = matrix->column;
            matrix->column = 0;
        }
        }
    }

    matrix->column = col;

}

void multiply(mtrx* A, mtrx* b)
{

    printf("\n\nThe result of the multiplication is.\n\n");
    double result;
    for(int i=0;i<A->row;i++){
        for(int j=0;j<b->column;j++){
            result = 0.;
            for(int m=0;m<A->column;m++){
                result += A->matrix[i*A->column + m] * b->matrix[j*b->row + m];
            }
            printf("%lf ",result);
        }
        printf("\n");
    }
    printf("\n\n");

}

void normalize(mtrx* matrix, mtrx* vector)
{
    // Normalization of the matrix.
    for(int i=0;i<matrix->row;i++)
    {
        double grts = 0.;
        for(int j=0;j<matrix->column;j++){
            if(grts<fabs(matrix->matrix[i*matrix->column + j])){
                grts = fabs(matrix->matrix[i*matrix->column + j]);
            }
        }
        double mult = 1/grts;
        for(int j=0;j<matrix->column;j++){
            matrix->matrix[i*matrix->column + j]*=mult;
        }
        vector->matrix[i]*=mult;
    }

}

void switchRows(mtrx* matrix,int row_1,int row_2)
{
    double temp[matrix->column];
    for(int i=0;i<matrix->column;i++){
        temp[i] = matrix->matrix[row_1*matrix->column + i];
        matrix->matrix[row_1*matrix->column + i] = matrix->matrix[row_2*matrix->column + i];
        matrix->matrix[row_2*matrix->column + i] = temp[i];
    }
}

void printOut(mtrx* matrix)
{
    printf("\n");
    for(int i=0;i<matrix->row;i++){
        for(int j=0;j<matrix->column;j++){
            printf("%lf ",matrix->matrix[i*matrix->column + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printOutBoth(mtrx* matrix, mtrx* vector)
{
    printf("\n");
    for(int i=0;i<matrix->row;i++){
        for(int j=0;j<matrix->column;j++){
            printf("%lf ", matrix->matrix[i*matrix->column + j]);
        }
        printf("\t %lf", vector->matrix[i]);
        printf("\n");
    }
    printf("\n");
}

void diagonal(mtrx* matrix, mtrx* vector, mtrx* order)
{
    for(int i=0;i<matrix->row;i++){
        if(matrix->matrix[i*matrix->column + i]!=1 && matrix->matrix[i*matrix->column + i]!=0){
            double mult = 1/matrix->matrix[i*matrix->column + i];
            matrix->matrix[i*matrix->column + i] = 1.;
            vector->matrix[i] *= mult;
            for(int j=0;j<matrix->column;j++){
                if(j!=i){
                    matrix->matrix[i*matrix->column + j] *= mult;
                }
            }
        }else if(matrix->matrix[i*matrix->column + i] == 0.){
            // Try to switch it with the biggest double below
            // So it will all come to switching rows after the
            // biggest is found.
            //printf("\nHere I am!\n");
            int ind = i;
            double grts = fabs(matrix->matrix[ind*matrix->column+i]);
            //printOut(matrix);
            //printf("%d %d\n",ind, i);
            //printf("%lf\n",grts);
            for(int k=i+1;k<matrix->row;k++){
                double next = matrix->matrix[k*matrix->column + i];
                if(fabs(next)>grts){
                    ind = k;
                }
            }
            // The biggest index is found. Now it's time to switch
            // rows i and ind if ind != i. If it is equal to it then
            // the matrix is singular the program exits with 1.
            if(ind!=i){
                switchRows(matrix,i,ind);
                switchRows(vector,i,ind);
                // Now that the diagonal element is not zero I need to set it to 1. here!
                // As well as the solution vector.
                // I need to make a function out of this, due to code duplication.
                double mult = 1/matrix->matrix[i*matrix->column + i];
                matrix->matrix[i*matrix->column + i] = 1.;
                vector->matrix[i] *= mult;
                for(int j=0;j<matrix->column;j++){
                    if(j!=i){
                        matrix->matrix[i*matrix->column + j] *= mult;
                    }
                }
                switchRows(order,i,ind);
            }else{
                printf("The matrix is singular. Exiting...\n");
                exit(1);
            }
        }

        eliminate(matrix,vector,i);
    }
}

void eliminate(mtrx* matrix, mtrx* vector, int ind)
{

    for(int i=ind+1;i<matrix->row;i++){
        if(fabs(matrix->matrix[i*matrix->column+ind])>eps){
            double mult = matrix->matrix[ind*matrix->column + ind]/matrix->matrix[i*matrix->column + ind];
            vector->matrix[i] *= mult;
            for(int j=0;j<matrix->column;j++){
                matrix->matrix[i*matrix->column + j] *= mult;
                matrix->matrix[i*matrix->column + j] -= matrix->matrix[ind*matrix->column + j];
            }
            vector->matrix[i] -= vector->matrix[ind];
        }else{
            // The values below are so small that I regard them
            // as true zeroes. This will ruin the solution but
            // is going to give an approximation to the problem.
            // ATTENTION. If the value is exactly zero I do nothing
            // that is actually the goal.
            if(matrix->matrix[i*matrix->column + ind] != 0.){
                if(flag==0){
                    printf("Numerical singularity occcured.\nATTENTION: The solution is only going to be an approximation!\n\n");
                    flag++;
                }
                matrix->matrix[i*matrix->column + ind] = 0.;
            }
        }
    }

}

void back_eliminate(mtrx* matrix, mtrx* vector)
{
    for(int i=1;i<=matrix->row-1;i++){
        int j=i-1;
        while(j>=0){
            if(fabs(matrix->matrix[j*matrix->column + i])>eps){
                double mult = matrix->matrix[i*matrix->column + i]/matrix->matrix[j*matrix->column + i];
                vector->matrix[j] *= mult;
                for(int k=0;k<matrix->column;k++){
                    matrix->matrix[j*matrix->column + k] *= mult;
                }
                vector->matrix[j] -= vector->matrix[i];
                for(int k=0;k<matrix->column;k++){
                    matrix->matrix[j*matrix->column + k] -= matrix->matrix[i*matrix->row + k];
                }
            }else{
                // Approximation again.
                if(flag==0 && matrix->matrix[j*matrix->column + i] != 0.){
                    printf("Numerical singularity occcured.\nATTENTION: The solution is only going to be an approximation!\n\n");
                    flag++;
                }
                matrix->matrix[j*matrix->column + i] = 0.;
            }
            j--;
        }
    }
}

void reOrder(mtrx* vector, mtrx* order)
{
    // It is not the minimum of rearranging but I don't think this
    // is too important for efficiency of this code.
    double max;
    for(int i=0;i<vector->row-1;i++){
        max = order->matrix[i];
        for(int j=i+1;j<vector->row;j++){
            if(order->matrix[j]<max){
                double temp = order->matrix[i];
                order->matrix[i] = order->matrix[j];
                order->matrix[j] = temp;
                temp = vector->matrix[i];
                vector->matrix[i] = vector->matrix[j];
                vector->matrix[j] = temp;
            }
        }
    }
}

void lin_solve(mtrx* matrix, mtrx* vector, mtrx* order)
{
    //printOutBoth(matrix,vector);
    //normalize(matrix,vector);
    //printOutBoth(matrix,vector);
    diagonal(matrix,vector,order);
    //printOutBoth(matrix,vector);
    back_eliminate(matrix,vector);
    //printOutBoth(matrix,vector);
    diagonal(matrix,vector,order);
    //printOutBoth(matrix,vector);
    reOrder(vector,order);
}

void transponse(mtrx* matrix){
    for(int i=0;i<matrix->row;i++){
        for(int j=i;j<matrix->column;j++){
            double temp = matrix->matrix[i*matrix->row + j];
            matrix->matrix[i*matrix->row + j] = matrix->matrix[j*matrix->row + i];
            matrix->matrix[j*matrix->row + i] = temp;
        }
    }
}

void get_inverse(mtrx* BASE,double** a){
    // I apply the Gauss elimination for all base vectors!
    // a will be the inverted matrix
    mtrx* A = (mtrx*)malloc(sizeof(mtrx*));
    A->matrix = (double*)malloc(sizeof(double)*BASE->length);
    A->row = BASE->row;
    A->column = BASE ->column;
    mtrx* sol = (mtrx*)malloc(sizeof(mtrx*));
    sol->matrix = (double*)malloc(sizeof(double)*BASE->row);
    sol->row = BASE->row;
    sol->column = 1;
    mtrx* order = (mtrx*)malloc(sizeof(mtrx*));
    order->matrix = (double*)malloc(sizeof(double)*BASE->row);
    order->row = BASE->row;
    order->column = 1;
    for(int iter=0;iter<BASE->row;iter++){
        for(int i=0;i<BASE->row;i++){
            sol->matrix[i] = 0;
            order->matrix[i] = i;
            for(int j=0;j<BASE->column;j++){
                A->matrix[i*BASE->row + j] = BASE->matrix[i*BASE->row + j];
            }
        }
        sol->matrix[iter] = 1.;
        // matrix and solution vector set, here comes the gauss elimination
        lin_solve(A,sol,order);
        for(int k=0;k<sol->row;k++){
            a[k][iter] = sol->matrix[k];
        }  
    }
    for(int i=0;i<BASE->row;i++){
        for(int j=0;j<BASE->column;j++){
            printf("%lf ", a[i][j]);
        }
        printf("\n");
    }
}

void lapack_inverse(mtrx* BASE, double** a){
    // I apply the Gauss elimination for all base vectors!
    // a will be the inverted matrix
    int info;
    mtrx* A = (mtrx*)malloc(sizeof(mtrx*));
    A->matrix = (double*)malloc(sizeof(double)*BASE->length);
    A->row = BASE->row;
    A->column = BASE ->column;
    mtrx* sol = (mtrx*)malloc(sizeof(mtrx*));
    sol->matrix = (double*)malloc(sizeof(double)*BASE->row);
    sol->row = BASE->row;
    sol->column = 1;
    for(int iter=0;iter<BASE->row;iter++){
        for(int i=0;i<BASE->row;i++){
            sol->matrix[i] = 0;
            for(int j=0;j<BASE->column;j++){
                A->matrix[i*BASE->row + j] = BASE->matrix[i*BASE->row + j];
            }
        }
        sol->matrix[iter] = 1.;
        transponse(A);
        // matrix and solution vector set, here comes the gauss elimination
        //FOR THE SVD
        int N = BASE->column;
        int nrhs = 1;
        int lda = BASE->row;
        int ipiv[BASE->column];
        int ldb = BASE->row;

        dgesv_(&N,&nrhs,A->matrix,&lda,ipiv,sol->matrix,&ldb,&info);

        if(info==0){
            // Copy the solution to the inverse
            for(int k=0;k<sol->row;k++){
                a[k][iter] = sol->matrix[k];
            } 
        }else{
            printf("SVD failed.\n");
        }
    }
    if(info==0){
        for(int i=0;i<BASE->row;i++){
            for(int j=0;j<BASE->column;j++){
                printf("%lf ", a[i][j]);
            }
            printf("\n");
        }
    }
}

void compare(double** a,double** b,int size){
    printf("\nThe absolute avarage difference between my calculations and Lapack's is: \n\n");
    float abs_diff = 0.;
    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            abs_diff = fabs(a[i][j]-b[i][j]);
        }
    }
    printf("%.10lf is the absolute value of difference avaraged by elements.",abs_diff/(float)size*size);
    printf("\n");
}
