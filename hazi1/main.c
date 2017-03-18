/*
 Első házi feladat. Mátrix szorzás! Teljes megoldás.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// It is not efficient on small matrices and vectors but it is more efficient on bigger ones.
int DEF_SIZE = 20;

// This is how I get the information about row and column count in a vector typed matrix.
typedef struct{
    int length;
    int row;
    int column;
    float* matrix;
}mtrx;

// The matrix and vector read in function. My program can accept matrices as well as vectors for input and
// multiply them if it's possible
void matrix_in(FILE*, mtrx*);

// The multiplication function.
void multiply(mtrx* A, mtrx* b);

int main(int argc, char* argv[])
{

    if(argc<3){
        printf("There should be at least 2 arguments.\n");
        exit(1);
    }

    mtrx* A = (mtrx*)malloc(sizeof(mtrx));
    mtrx* b = (mtrx*)malloc(sizeof(mtrx));

    FILE* fp = fopen(argv[1],"r");
    if(!fp){
        printf("An error occured loading %s \n", argv[1]);
        exit(1);
    }

    A->matrix = (float*)malloc(DEF_SIZE*sizeof(float));
    matrix_in(fp,A);
    
    fclose(fp);
    
    if(!(fp=fopen(argv[2],"r"))){
        printf("An error occured loading %s \n", argv[2]);
        exit(1);
    }

    DEF_SIZE = 10;
    b->matrix = (float*)malloc(DEF_SIZE*sizeof(float));
    matrix_in(fp,b);

    if(b->row!=A->column){
        printf("Incorrect dimensions.\n");
        exit(1);
    }

    // I don't save the result! I just print it out. This way I don't need to store it!
    // I could write it to a file right away.
    multiply(A,b);

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
        char delims[] = " ";
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
                    matrix->matrix = (float*)realloc(matrix->matrix,sizeof(float)*DEF_SIZE);
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

void multiply(mtrx* A, mtrx* b){
    
    printf("\n\nThe result of the multiplication is.\n\n");
    float result;
    for(int i=0;i<A->row;i++){
        for(int j=0;j<b->column;j++){
            result = 0.;
            for(int m=0;m<A->column;m++){
                result += A->matrix[i*A->column + m] * b->matrix[j*b->row + m];
            }
            printf("%f ",result);
        }
        printf("\n");
    }
    printf("\n\n");

}
