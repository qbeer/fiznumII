#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

unsigned int flag = 0;

const double eps = 10e-10;

unsigned int DEF_SIZE = 20;

typedef struct{
    int length;
    int row;
    int column;
    double* matrix;
}mtrx;

void matrix_in(FILE*, mtrx*);
void multiply(mtrx* A, mtrx* b);
void normalize(mtrx*,mtrx*);
void switchRows(mtrx*,int,int);
void printOut(mtrx*);
void printOutBoth(mtrx*, mtrx*);
void eliminate(mtrx*, mtrx*, int);
void back_eliminate(mtrx*, mtrx*);
void lin_solve(mtrx*,mtrx*);
void get_inverse(mtrx*,double**);
void fill_van_der_monde(mtrx*,mtrx*,int,int);
void trans_times_fill(mtrx*,mtrx*);
void factor_error(mtrx*, mtrx*);
void print_to_file(FILE*,mtrx*);

int main(int argc, char* argv[])
{

    if(argc<4){
        printf("\tYou should provide the filename,\n\t\tnumber of independent variables \n\t\t\tand the order of the desired fit.\n");
        exit(1);
    }

    mtrx* data = (mtrx*)malloc(sizeof(mtrx));
    FILE* fp = fopen(argv[1],"r");

    if(fp == NULL){
        printf("Error opening file %s\n",argv[1]);
    }

    // Always need to allocate some space before read in!!!
    data->matrix = (double*)malloc(DEF_SIZE*sizeof(double));
    matrix_in(fp,data);
    // reset DEF_SIZE for future readins
    DEF_SIZE = 20;
    // DATA READ IN CORRECTLY
    fclose(fp);

    int independent_val = atoi(argv[2]);
    int order = atoi(argv[3]);

    mtrx* sol = (mtrx*)malloc(sizeof(mtrx));
    mtrx* error = (mtrx*)malloc(sizeof(mtrx));
    sol->row = data->row;
    error->row = data->row;
    sol->column = 1;
    error->column = 1;
    sol->length = data->row;
    error->length = data->row;
    sol->matrix = (double*)malloc(sizeof(double)*data->row);
    error->matrix = (double*)malloc(sizeof(double)*data->row);

    // Fill the solution vector

    int loc = 0;
    for(int i=0;i<data->row;i++){
        // Here I suppose that the last element is some kind of error
        // and the one before that is the value of y_k !!!
        sol->matrix[i] = data->matrix[i*data->column + data->column-2];
        error->matrix[i] = data->matrix[i*data->column + data->column-1];
    }
    // Values filled correctly.

    factor_error(sol,error);

    mtrx* X = (mtrx*)malloc(sizeof(mtrx));
    // Data row is the number of rows in the data file
    // the independent values provided.
    X->row = data->row;
    // Column size should be the independt data with
    // times the order of the polinomial + 1 for the constant
    // value.
    X->column = independent_val*order + 1;
    // Number of elements
    X->length = X->column*X->row;
    X->matrix = (double*)malloc(X->length*sizeof(double));

    fill_van_der_monde(data,X,independent_val, order);
    // Values filled correctly.

    factor_error(X,error);

    // I don't need the data anymore.
    free(data->matrix);
    free(data);

    // I'll use the data struct to hold the value of the X^T*X
    mtrx* XTX = (mtrx*)malloc(sizeof(mtrx));
    XTX->column = X->column;
    XTX->row = X->column;
    XTX->length = XTX->row*XTX->column;
    XTX->matrix = (double*)malloc(sizeof(double)*XTX->length);

    trans_times_fill(X,XTX);

    mtrx* XTy = (mtrx*)malloc(sizeof(mtrx));
    XTy->column = 1;
    XTy->row = X->column;
    XTy->length = X->column;
    XTy->matrix = (double*)malloc(sizeof(double)*XTy->length);

    for(int i=0;i<XTy->row;i++){
        for(int j=0;j<XTy->column;j++){
            XTy->matrix[i*XTy->column+j] = 0.;
            for(int k=0;k<sol->row;k++){
                XTy->matrix[i] += X->matrix[k*X->column + i]*sol->matrix[k];
            }
        }
    }

    lin_solve(XTX,XTy);

    printf("\n\n THE SOLUTION IS: \n\n");
    printOut(XTy);

    printf("\n");

    const char* o_file = argv[1];
    const char* ext0_5 = "_";
    const char* ext1 = argv[2];
    const char* ext1_5 = "_independent_val_";
    const char* ext2 = argv[3];
    const char* ext3 = "_order_fit_vector.dat";

    char* o_file_name;
    o_file_name = malloc(strlen(o_file)+strlen(ext1)+strlen(ext2)+strlen(ext3)+strlen(ext1_5)+strlen(ext0_5)+1);
    strcpy(o_file_name,o_file);
    strcat(o_file_name,ext0_5);
    strcat(o_file_name,ext1);
    strcat(o_file_name,ext1_5);
    strcat(o_file_name,ext2);
    strcat(o_file_name,ext3);

    fp = fopen(o_file_name,"w");

    print_to_file(fp,XTy);

    fclose(fp);

    char* o_file_name_2;
    ext3 = "_order_fit_data.dat";

    o_file_name_2 = malloc(strlen(o_file)+strlen(ext1)+strlen(ext2)+strlen(ext3)+strlen(ext1_5)+strlen(ext0_5)+1);
    strcpy(o_file_name_2,o_file);
    strcat(o_file_name_2,ext0_5);
    strcat(o_file_name_2,ext1);
    strcat(o_file_name_2,ext1_5);
    strcat(o_file_name_2,ext2);
    strcat(o_file_name_2,ext3);

    fp = fopen(o_file_name_2,"w");

    for(int i=0;i<error->row;i++){
        error->matrix[i] = 1/(error->matrix[i]);
    }

    // defactor error
    factor_error(X,error);
    factor_error(sol,error);

    for(int i=0;i<X->row;i++){
        fprintf(fp,"%lf\t",sol->matrix[i]);
        double result = 0.0;
        for(int j=0;j<X->column;j++){
            result += X->matrix[i*X->column + j]*XTy->matrix[j];
        }
        fprintf(fp,"%lf\n",result);
    }

    fclose(fp);

    free(X->matrix);
    free(X);
    free(XTy->matrix);
    free(XTy);
    free(sol->matrix);
    free(sol);
    free(error->matrix);
    free(error);
    free(XTX->matrix);
    free(XTX);

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

    while(fgets(line, 10000,fp)){
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

void diagonal(mtrx* matrix, mtrx* vector)
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
            if(fabs(matrix->matrix[i*matrix->column + ind]) < eps){
                if(flag==0){
                    printf("Numerical singularity occcured.\nATTENTION: The solution is only going to be an approximation!\n\n");
                    flag++;
                    //exit(1);
                }
                //matrix->matrix[i*matrix->column + ind] = 0.;
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
                if(flag==0 && fabs(matrix->matrix[j*matrix->column + i]) < eps){
                    printf("Numerical singularity occcured.\nATTENTION: The solution is only going to be an approximation!\n\n");
                    flag++;
                    //exit(1);
                }
                //matrix->matrix[j*matrix->column + i] = 0.;
            }
            j--;
        }
    }
}

void lin_solve(mtrx* matrix, mtrx* vector)
{
    normalize(matrix,vector);
    diagonal(matrix,vector);
    back_eliminate(matrix,vector);
    diagonal(matrix,vector);
    printOutBoth(matrix,vector);
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
        lin_solve(A,sol);
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

void fill_van_der_monde(mtrx* data, mtrx* X,int ind_val, int order){
    for(int j=0;j<data->row;j++){
        X->matrix[j*X->column + 0] = 1.;
    }
    for(int i=0;i<data->row;i++){
        int loc = 1;
        for(int j=0;j<ind_val;j++){
            for(int k=1;k<=order;k++){
                X->matrix[i*X->column + loc] = pow(data->matrix[i*data->column + j],k);
                loc++;
            }
        }
    }
}

void factor_error(mtrx* mtr, mtrx* error){
    for(int i=0;i<mtr->row;i++){
        for(int j=0;j<mtr->column;j++){
            mtr->matrix[i*mtr->column+j]*=error->matrix[i];
        }
    }
}

void trans_times_fill(mtrx* X,mtrx* out){
    for(int i=0;i<out->row;i++){
        for(int j=0;j<out->column;j++){
            out->matrix[i*out->column+j] = 0.;
            for(int l=0;l<X->row;l++){
                out->matrix[i*out->column+j] += X->matrix[l*X->column + i]*X->matrix[l*X->column+j];
            }
        }
    }
}

void print_to_file(FILE* fp,mtrx* mat){
    for(int i=0;i<mat->row;i++){
        for(int j=0;j<mat->column;j++){
            fprintf(fp,"%lf ",mat->matrix[i*mat->column + j]);
        }
        fprintf(fp,"\n");
    }
}
