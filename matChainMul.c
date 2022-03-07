#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdbool.h>

void matMul (
    unsigned int l,     // Rows of m1
    unsigned int m,     // Cols of m1 == Rows of m2
    unsigned int n,     // Cols of l2
    int** matrix_a,     // m1
    int** matrix_b,     // m2
    int** matMulProduct // result 
) {

    for ( unsigned int i=0; i<l; i++ ) {                                 // for every row in m1/product matrix
        for ( unsigned int k=0; k<n; k++ ) {                             // for every col in m2/product matrix
            matMulProduct[i][k] = 0;                                     // set the corresponding product location to 0
                                                                         // NOTE: not memory error here
            for ( unsigned int j=0; j<m; j++ ) {                         // for every position that need multiplication
                matMulProduct[i][k] += matrix_a[i][j] * matrix_b[j][k];  // When the lookup throws an error, it must be due to ununitialized matrix a and b. They should be equal to a value from matrices.
            }
        }
    }
}

unsigned int cost (
    unsigned int matrixCount,
    unsigned int* rowSizes,
    unsigned int* colSizes
) {
    if ( matrixCount==1 ) { return 0; } // Base case, no multplication to be done
    else {

        unsigned int numPossibleSplits = matrixCount-1; 

        unsigned int costs[numPossibleSplits];
        for ( unsigned int split=0; split<numPossibleSplits; split++ ) {

            unsigned int l = rowSizes[0];
            assert ( colSizes[split] == rowSizes[split+1] );
            unsigned int m = colSizes[split];
            unsigned int n = colSizes[matrixCount-1];

            costs[split] =
                cost( split+1, rowSizes, colSizes ) +                            // cost of left subchain
                l * m * n +                                                      // cost of multiplying the two chains
                cost( matrixCount-split-1, rowSizes+split+1, colSizes+split+1 ); // cost of right subchain
        }

        unsigned int minCost = costs[0];
        for ( unsigned int split=1; split<numPossibleSplits; split++ ) {
            if ( costs[split]<minCost ) {
                minCost = costs[split];
            }
        }

        return minCost;
    }
}

void matChainMul(
    unsigned int matrixCount,
    unsigned int* rowSizes,
    unsigned int* colSizes,
    int*** matrices,
    int** matChainMulProduct
) {  
    // copied from cost function
    if (matrixCount==1) {                                // base case.
        // printf("new MATRIX #####\n");
        for (unsigned int i=0; i<rowSizes[0]; i++) {     //initialize and populate the rows of the product matrix 
            for (unsigned int j =0; j<colSizes[0]; j++) {
                matChainMulProduct[i][j] = matrices[0][i][j];
            }
        }                                
        
    } else {
        unsigned int numPossibleSplits = matrixCount-1; // number of possible splits
        
        // declare and initialize an array for the cost of splitting at any point
        unsigned int costs[numPossibleSplits];
        
        // populate the costs array
        for ( unsigned int split=0; split<numPossibleSplits; split++ ) {
            unsigned int l = rowSizes[0];                    // l is the number of rows at the start slit'th matrix section
            assert ( colSizes[split] == rowSizes[split+1] ); // verify that the rows and colomns of adjacent matrices line-up
            unsigned int m = colSizes[split];                // m is the number of cols in the split'th matrix and rows in the split+1'th matrix
            unsigned int n = colSizes[matrixCount-1];        // n is the number of cols at the end of the split+1'th matrix line

            costs[split] =
                cost( split+1, rowSizes, colSizes ) +                            // cost of left subchain
                l * m * n +                                                      // cost of multiplying the two chains
                cost( matrixCount-split-1, rowSizes+split+1, colSizes+split+1 ); // cost of right subchain
        }

        // With a populated costs array, we seek the best cost and split location
        unsigned int minCost = costs[0];
        int bestSplit = 0;
        for ( unsigned int split=1; split<numPossibleSplits; split++ ) {
            if ( costs[split]<minCost ) {
                minCost = costs[split];
                bestSplit = split;
            }
        }
        
        int** leftSubchainProduct = (int**) malloc(sizeof(int*) * rowSizes[0]);
        for (unsigned int i=0; i<rowSizes[0]; i++) {     //initialize and populate the rows of the product matrix 
            leftSubchainProduct[i] = (int*) malloc(sizeof(int)*colSizes[bestSplit]);
        }
        matChainMul(bestSplit+1, rowSizes, colSizes, matrices, leftSubchainProduct);

        int** rightSubchainProduct = (int**) malloc(sizeof(int*) * rowSizes[(bestSplit+1)]); 
        for (unsigned int i=0; i<rowSizes[bestSplit+1]; i++) {     //initialize and populate the rows of the product matrix 
            rightSubchainProduct[i] = (int*) malloc(sizeof(int)*colSizes[matrixCount-1]);
        }
        matChainMul(matrixCount-(bestSplit+1), rowSizes+(bestSplit+1), colSizes+(bestSplit+1), matrices+(bestSplit+1), rightSubchainProduct);

        unsigned int l = rowSizes[0];               // number of rows in matrix a
        unsigned int m = colSizes[bestSplit];       // number of cols in matrix a        
        unsigned int n = colSizes[matrixCount-1];   // number of cols in matrix b
        
        matMul(l, m, n, leftSubchainProduct, rightSubchainProduct, matChainMulProduct);
        
        /*
        printf("PROD MATRIX ###### l=%d m=%d n=%d \n", l, m, n);
        for (unsigned int i=0; i<rowSizes[0]; i++) {     //initialize and populate the rows of the product matrix 
            for (unsigned int j=0; j<colSizes[matrixCount-1]; j++) {
                printf("%d ", matChainMulProduct[i][j]);
            }
            printf("\n");
        }
        */
        
        // Freeing memory
        
        for (unsigned int i=0; i<rowSizes[0]; i++) {     //initialize and populate the rows of the product matrix 
            free(leftSubchainProduct[i]);
        }
        free(leftSubchainProduct);
        
        for (unsigned int i=0; i<rowSizes[bestSplit+1]; i++) {     //initialize and populate the rows of the product matrix 
            free(rightSubchainProduct[i]);
        }
        free(rightSubchainProduct);
    }
}

int main(int argc, char* argv[]) {
    unsigned int matrixCount;
    unsigned int* rowSizes;
    unsigned int* colSizes;
    int*** matrices;
    int** productMatrix;

    // Open File
    FILE* fp = fopen(argv[1], "r");
    if (!fp) {
        perror("fopen failed");
        exit(EXIT_FAILURE);
    }

    // Save l1 to matrixCount
    if (!fscanf(fp, "%d\n", &matrixCount)) {
        perror("reading the number of matrices failed");
        exit(EXIT_FAILURE);
    }

    // intitialize rowSize, colSize, matrices
    rowSizes = calloc( matrixCount, sizeof(int) );
    colSizes = calloc( matrixCount, sizeof(int) );
    matrices = (int***) malloc(sizeof(int**)*matrixCount);

    // Fill rowSizes, colSizes, matrices
    for (unsigned int matIndex=0; matIndex<matrixCount; matIndex++) {

        // Save rows, cols to rowSizes, colSizes
        unsigned int rows, cols;
        if (!fscanf(fp, "%d %d\n", &rows, &cols)) {
            perror("reading the dimensions of matrix failed");
            exit(EXIT_FAILURE);
        }
        rowSizes[matIndex] = rows;
        colSizes[matIndex] = cols;

        // Initialize this matrix
        matrices[matIndex] = (int**) malloc(sizeof(int*)* rowSizes[matIndex]);

        // populate this matrix
        for (unsigned int i = 0; i < rows; i++) {
            matrices[matIndex][i] = (int*) malloc(sizeof(int)*cols);
            for (unsigned int j = 0; j < cols; j++) {
                int elem;
                fscanf(fp, "%d", &elem);
                matrices[matIndex][i][j] = elem;
            }
        }
    }

    // Initialize productMatrix: rows = rowSizes[0], cols = colSizes[matrixCount-1]
    productMatrix = (int**) malloc(sizeof(int*)*rowSizes[0]);
    
    // Initialize the rows of productMatrix
    for (unsigned int k = 0; k < rowSizes[0]; k++) {
        productMatrix[k] = (int*) malloc(sizeof(int)*(colSizes[matrixCount-1]));
    }

    // cunduct chain multiplication on matrices
    matChainMul(matrixCount, rowSizes, colSizes, matrices, productMatrix);
    

    //print cost(matrices), THEN the rest
    // printf("SOLUTION\n");
    int total_cost = cost(matrixCount, rowSizes, colSizes);
    printf("%d\n", total_cost);

    for (unsigned int i=0; i<rowSizes[0]; i++) {
        for (unsigned int j=0; j<colSizes[matrixCount-1]; j++) {
            printf("%d ", productMatrix[i][j]);
        }
        printf("\n");
    }

    // free matrices 
    for (unsigned int m = 0; m < matrixCount; m++) {
        for (unsigned int i = 0; i < rowSizes[m]; i++) {
            free(matrices[m][i]);
        }
        free(matrices[m]);
    }
    free(matrices);

    // free productMatrix
    for (unsigned int i = 0; i < rowSizes[0]; i++) {
        free(productMatrix[i]);
    }
    free(productMatrix);

    // free memory for everything
    free(colSizes);
    free(rowSizes);

    // close file
    fclose(fp);
    exit(EXIT_SUCCESS);
}