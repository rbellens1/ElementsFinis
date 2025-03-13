#include "fem.h"

# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femGeo* theGeometry = geoMeshCreate(filename);
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->geo  = theGeometry;
    femMesh *theMesh = theGeometry->theElements;
    if (theMesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theMesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theMesh->nodes->nNodes);
    return theProblem;
}

# endif
# ifndef NOPOISSONBOUNDARY

void femPoissonFindBoundaryNodes(femPoissonProblem *theProblem)
{
    femGeo* theGeometry = theProblem->geo;  
    femMesh* theEdges = theGeometry->theEdges; 
    int nBoundary = 0;
    
    for (int i = 0; i < theEdges->nElem; i++) {
        if (theEdges->elem[i*nBoundary] == 0) {
            nBoundary++;
        }
    }

    
    

    femDomain *theBoundary = malloc(sizeof(femDomain));
    theGeometry->nDomains++;
    theGeometry->theDomains = realloc(theGeometry->theDomains,theGeometry->nDomains*sizeof(femDomain*));
    theGeometry->theDomains[theGeometry->nDomains-1] = theBoundary;
    theBoundary->nElem = nBoundary;
    theBoundary->elem = malloc(nBoundary*sizeof(int));
    theBoundary->mesh = NULL;
    sprintf(theBoundary->name,"Boundary");
 
    for (int i = 0; i < theEdges->nElem; i++) {
        if (theEdges->elem[i*nBoundary] == 0) {
            theBoundary->elem[i] = i; 
        }
    }

}
    
# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    geoMeshFree(theProblem->geo);
    free(theProblem);
    
}
    
# endif
# ifndef NOPOISSONLOCAL

void femPoissonLocal(femPoissonProblem *theProblem, const int iElem, int *map, double *x, double *y)
{
    femMesh *theMesh = theProblem->geo->theElements;
    femDiscrete *theSpace = theProblem->space;
    int nodes = theMesh->nLocalNode;
    for (int i = 0; i < nodes; i++) {
        int j = theMesh->elem[iElem*nodes+i];
        x[i] = theMesh->nodes->X[j];
        y[i] = theMesh->nodes->Y[j]; 
        map[i] = j;
    }
}

# endif
# ifndef NOPOISSONSOLVE

void femPoissonSolve(femPoissonProblem *theProblem)
{

    femMesh *theMesh = theProblem->geo->theElements;
    femDomain *theBoundary = geoGetDomain(theProblem->geo,"Boundary");
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
 
    int nLocalNode = theMesh->nLocalNode;

    double **A = theSystem->A;
    double *B = theSystem->B;

    for(int iElem = 0; iElem < theMesh->nElem; iElem++){
        double *x = malloc(nLocalNode*sizeof(double));
        double *y = malloc(nLocalNode*sizeof(double));
        int *map = malloc(nLocalNode*sizeof(int));

        femPoissonLocal(theProblem,iElem,map,x,y);

        for(int integral = 0; integral < theRule->n; integral++){
            double weight = theRule->weight[integral];
            double xsi = theRule->xsi[integral];
            double eta = theRule->eta[integral];

            double *phi = malloc(theSpace->n * sizeof(double));
            femDiscretePhi2(theSpace,xsi,eta,phi);

            double *dphidxsi = malloc(theSpace->n * sizeof(double));
            double *dphideta = malloc(theSpace->n * sizeof(double));
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);

            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;

            for(int j = 0; j<theSpace->n; j++){
                dxdxsi += dphidxsi[j]*x[j];
                dxdeta += dphideta[j]*x[j];
                dydxsi += dphidxsi[j]*y[j];
                dydeta += dphideta[j]*y[j];
            }

            double Jacob = dxdxsi*dydeta - dxdeta*dydxsi;

            if(Jacob < 0) {
                int node = theMesh->elem[iElem*nLocalNode];
                theMesh->elem[iElem*nLocalNode] = theMesh->elem[nLocalNode*iElem+2];
                theMesh->elem[nLocalNode*iElem+2] = node;
            }

            Jacob = fabs(Jacob);

            double *dphidx = malloc(theSpace->n * sizeof(double));
            double *dphidy = malloc(theSpace->n * sizeof(double));

            for(int i = 0; i<theSpace->n; i++){
                dphidx[i] = (dphidxsi[i]*dydeta - dphideta[i]*dydxsi)/Jacob;
                dphidy[i] = (dphideta[i]*dxdxsi - dphidxsi[i]*dxdeta)/Jacob;
            }

            for(int i = 0; i<theSpace->n; i++){
                B[map[i]] += weight*phi[i]*Jacob;
                for(int j = 0; j<theSpace->n; j++){
                    A[map[i]][map[j]] += weight*(dphidx[i]*dphidx[j] + dphidy[i]*dphidy[j])*Jacob;
                }
            }
            free(dphidxsi);
            free(dphideta);
            free(phi);
            free(dphidx);
            free(dphidy);
        }

        free(x);
        free(y);
        free(map);
    }
    femMesh *theEdges = theProblem->geo->theEdges;
    for(int i = 0; i < theEdges->nElem; i++){
        for(int j = 0; j < theEdges->nLocalNode; j++){
            femFullSystemConstrain(theSystem,theEdges->elem[i*theEdges->nLocalNode+j],0.0);
        }
    }
    
    femFullSystemEliminate(theSystem);
    
}

# endif



