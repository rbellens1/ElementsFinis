#include "fem.h"
#define M_PI 3.14159265358979323846



typedef struct {
    int size;
    double **A;
} syscopy;


syscopy *A_cpy;







void femElasticityAssembleElementsAXI(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]]; // coordonnée r
            y[j]    = theNodes->Y[map[j]]; // coordonnée z
        } 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            double r = 0.0;  // Coordonnée radiale au point d'intégration
            
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i];
                r += x[i]*phi[i];  // Calcul de r au point d'intégration
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            // En axisymétrie, on intègre sur 2πr
            double dV = 2.0 * M_PI * r * jac * weight;
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }
            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                           dphidy[i] * c * dphidy[j]) * dV +
                                            dphidx[i] * b * phi[j] +
                                            phi[i] * (b * dphidx[j] + a * phi[j] / r);
                                           
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                           dphidy[i] * c * dphidx[j]) * dV +
                                           phi[i] * b * dphidy[j];
                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                           dphidx[i] * c * dphidy[j]) * dV +
                                           dphidy[i] * b * phi[j];
                                           
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                           dphidx[i] * c * dphidx[j]) * dV;
                }
                B[mapY[i]] -= phi[i] * g * rho * dV;
            }
        }
    }
}




















void femElasticityAssembleElements(femProblem *theProblem){

    if (theProblem->iCase == AXISYM)
    {
        femElasticityAssembleElementsAXI(theProblem);
        return;
    }


    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                            dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                            dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                            dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                            dphidx[i] * c * dphidx[j]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}}

}


void femElasticityAssembleNeumann(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->ruleEdge;
    femDiscrete    *theSpace = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[2],y[2],phi[2];
    int iBnd,iElem,iInteg,iEdge,i,j,d,map[2],mapU[2];
    int nLocal = 2;
    double *B  = theSystem->B;

    for(iBnd=0; iBnd < theProblem->nBoundaryConditions; iBnd++){
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        double value = theCondition->value;

        int shift=-1;
        if (type == NEUMANN_X)  shift = 0;      
        if (type == NEUMANN_Y)  shift = 1;  
        if (shift == -1) continue; 
        for(iEdge=0; iEdge < theCondition->domain->nElem; iEdge++){
            iElem = theCondition->domain->elem[iEdge];
            for (j=0; j < nLocal; j++) {
                map[j]  = theEdges->elem[iElem*nLocal+j];
                mapU[j] = 2*map[j] + shift;
                x[j]    = theNodes->X[map[j]];
                y[j]    = theNodes->Y[map[j]];} 
 
            double jac = sqrt((x[1]-x[0])*(x[1]-x[0]) + (y[1]-y[0])*(y[1]-y[0]))/2.0;
            for (iInteg=0; iInteg < theRule->n; iInteg++) {    
                double xsi    = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];  
                femDiscretePhi(theSpace,xsi,phi);
                for (i = 0; i < theSpace->n; i++) {    
                    B[mapU[i]] += jac * weight * phi[i] * value; }}}

    }
}



double *femElasticitySolve(femProblem *theProblem){
 
    //       
    // A completer :-) 
    // 

    double *residuals = theProblem->residuals;
    femFullSystemInit(theProblem->system);
    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem);
    int size = theProblem->system->size;

    A_cpy = malloc(sizeof(syscopy));
    if (A_cpy == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    A_cpy->size = size;
    A_cpy->A = (double **)  malloc(sizeof(double*) * size);

    for (int i = 0; i <size; i++) {
        A_cpy->A[i] =  malloc(sizeof(double) * size );
        memcpy(A_cpy->A[i], theProblem->system->A[i], sizeof(double) * size);
    }

    double *B = theProblem->system->B;
    for (int i = 0; i < size; i++)
    {
        residuals[i] = -B[i];
    }

    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theProblem->system,i,value); }}

    femFullSystemEliminate(theProblem->system);

    memcpy(theProblem->soluce, theProblem->system->B, sizeof(double) * size);
    
     return theProblem->soluce;
}

double * femElasticityForces(femProblem *theProblem){        
           
    //       
    // A completer :-) 
    //  


    femFullSystem  *theSystem = theProblem->system;
    //femFullSystemInit(theSystem);
    //femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem);


    double **A = A_cpy->A;
    double *soluce = theProblem->soluce;
    double *residuals = theProblem->residuals;
    //double *B = theProblem->system->B;
    int size = A_cpy->size;

    for (int i = 0; i < size; i++)
    {
        //residuals[i] = 0.0;
        for (int j = 0; j < size; j++)
        {
            residuals[i] += A[i][j] * soluce[j];
        }
    }


    for (int i = 0; i < A_cpy->size; i++)
    {
        free(A_cpy->A[i]);
        A_cpy->A[i] = NULL;
    }
    free(A_cpy->A);
    free(A_cpy);

    return theProblem->residuals;
}


