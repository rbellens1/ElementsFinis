#include"fem.h"
#include <stdlib.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#ifndef NORENUMBER 




typedef struct {
    int index;
    double value;
} NodeValue;

int compareNodeValueX(const void *a, const void *b) {
    NodeValue *nodeA = (NodeValue *)a;
    NodeValue *nodeB = (NodeValue *)b;
    return (nodeA->value > nodeB->value) - (nodeA->value < nodeB->value);
}

int compareNodeValueY(const void *a, const void *b) {
    NodeValue *nodeA = (NodeValue *)a;
    NodeValue *nodeB = (NodeValue *)b;
    return (nodeA->value > nodeB->value) - (nodeA->value < nodeB->value);
}

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i;
    int nNodes = theMesh->nodes->nNodes;
    NodeValue *nodeValues = malloc(nNodes * sizeof(NodeValue));
    
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[i] = i;
            free(nodeValues);
            break;
// 
// A modifier :-)
// debut
//
        case FEM_XNUM : 
            for (i = 0; i < nNodes; i++) {
                nodeValues[i].index = i;
                nodeValues[i].value = theMesh->nodes->X[i];
            }
            qsort(nodeValues, nNodes, sizeof(NodeValue), compareNodeValueX);
            for (i = 0; i < nNodes; i++) 
                theMesh->nodes->number[nodeValues[i].index] = i;
            free(nodeValues);
            break;

        case FEM_YNUM : 
            for (i = 0; i < nNodes; i++) {
                nodeValues[i].index = i;
                nodeValues[i].value = theMesh->nodes->Y[i];
            }
            qsort(nodeValues, nNodes, sizeof(NodeValue), compareNodeValueY);
            for (i = 0; i < nNodes; i++) 
                theMesh->nodes->number[nodeValues[i].index] = i;
            free(nodeValues);
            break;            
// 
// end
//
        default : Error("Unexpected renumbering option"); }
}

#endif
#ifndef NOBAND 

int femMeshComputeBand(femMesh *theMesh)
{
    int myBand = 0;
    int i, j, k;
    int nLocalNodes = theMesh->nLocalNode;
    int nElements = theMesh->nElem;
    int *elem = theMesh->elem;
    int *number = theMesh->nodes->number;

    for (i = 0; i < nElements; i++) {
        for (j = 0; j < nLocalNodes; j++) {
            for (k = 0; k < nLocalNodes; k++) {
                int node1 = number[elem[i * nLocalNodes + j]];
                int node2 = number[elem[i * nLocalNodes + k]];
                int distance = abs(node1 - node2);
                if (distance > myBand) {
                    myBand = distance;
                }
            }
        }
    }

    return(myBand + 1);
}


#endif
#ifndef NOBANDASSEMBLE


void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    // A Ecrire :-)
    int mapi;
    for (int i = 0; i < nLoc; i++)
    {
        mapi = map[i];
        for (int j = 0; j < nLoc; j++)
        {
            if (mapi <= map[j])
            {
                myBandSystem->A[mapi][map[j]] += Aloc[i * nLoc + j];
            }
            
        }
        myBandSystem->B[mapi] += Bloc[i];
    }

}


#endif
#ifndef NOBANDELIMINATE


double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;
    
    // A completer :-)


    for (k = 0; k < size; k++)
    {
        if ( fabs(A[k][k]) <= 1e-8 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);
            Error("Cannot eliminate with such a pivot");
        }
        // jend est la dernière ligne de la bande
        jend = MIN(k+band, size);
        for (i = k + 1; i < jend; i++)
        {
            factor = A[k][i] / A[k][k];
            for (j = i; j < jend; j++)
            {
                A[i][j] = A[i][j] - A[k][j] * factor;
            }
            B[i] -= B[k] * factor;
        }
    }

    // Backward substitution
    for (i = size-1; i >= 0 ; i--) {
        factor = 0;
        jend = MIN(i+band, size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i];
    }

    return(myBand->B);
}


#endif

