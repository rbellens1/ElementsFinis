#include "fem.h"




void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    
    int ierr;
    double r = w/4;
    int idRect = gmshModelOccAddRectangle(0.0,0.0,0.0,w,h,-1,0.0,&ierr); 
    int idDisk = gmshModelOccAddDisk(w/2.0,h/2.0,0.0,r,r,-1,NULL,0,NULL,0,&ierr); 
    int idSlit = gmshModelOccAddRectangle(w/2.0,h/2.0-r,0.0,w,2.0*r,-1,0.0,&ierr); 
    int rect[] = {2,idRect};
    int disk[] = {2,idDisk};
    int slit[] = {2,idSlit};

    gmshModelOccCut(rect,2,disk,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccCut(rect,2,slit,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccSynchronize(&ierr); 

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",11,&ierr);  
        gmshOptionSetNumber("Mesh.SmoothRatio", 21.5, &ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
 
    return;
}


double *femElasticitySolve(femProblem *theProblem)
{

    femFullSystem  *theSystem = theProblem->system; // Systeme
    femIntegration *theRule = theProblem->rule; // Points d'integration
    femDiscrete    *theSpace = theProblem->space; // Espace des fonctions de forme
    femGeo         *theGeometry = theProblem->geometry;  // Geometrie
    femNodes       *theNodes = theGeometry->theNodes; // Noeuds
    femMesh        *theMesh = theGeometry->theElements; // Maillage
    
    
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
    
    
  //
  //  A faire :-)
  //         
    //  -1- Boucle sur les elements
    //  -2- Boucle sur les noeuds du maillage
    //  -3- Boucle sur les points d'integration
    //  -4- Boucle sur les fonctions de forme
    //  -5- Boucle sur les derivees des fonctions de forme

    for(iElem = 0;iElem < theMesh->nElem; iElem++){ // Boucle sur les elements
        for (j = 0; j < nLocal; j++){ // Boucle sur les noeuds du maillage
            map[j] = theMesh->elem[nLocal*iElem+j]; // Numero du noeud j de l'element iElem
            x[j] = theMesh->nodes->X[map[j]]; // Coordonnee x du noeud j de l'element iElem
            y[j] = theMesh->nodes->Y[map[j]]; // Coordonnee y du noeud j de l'element iElem
        }
        for(iInteg = 0; iInteg < theRule->n; iInteg++){ // Boucle sur les points d'integration
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace,xsi,eta,phi); // Calcul des fonctions de forme
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta); // Calcul des derivees des fonctions de forme
            double xLocal = 0.0; // Coordonnee x du point d'integration
            double yLocal = 0.0; // Coordonnee y du point d'integration
            double dxdxsi = 0.0; // Derivee de x par rapport a xsi
            double dydxsi = 0.0; // Derivee de y par rapport a xsi
            double dxdeta = 0.0; // Derivee de x par rapport a eta
            double dydeta = 0.0; // Derivee de y par rapport a eta
            for(i = 0; i < theSpace->n; i++){ // Boucle sur les fonctions de forme
                xLocal += x[i]*phi[i]; // Coordonnee x du point d'integration
                yLocal += y[i]*phi[i]; // Coordonnee y du point d'integration
                dxdxsi += x[i]*dphidxsi[i]; // Derivee de x par rapport a xsi
                dydxsi += y[i]*dphidxsi[i]; // Derivee de y par rapport a xsi
                dxdeta += x[i]*dphideta[i]; // Derivee de x par rapport a eta
                dydeta += y[i]*dphideta[i]; // Derivee de y par rapport a eta
            } 
            double jacobian = fabs(dxdxsi*dydeta - dydxsi*dxdeta); 
            for(i = 0; i < theSpace->n; i++){  // Boucle sur les derivees des fonctions de forme
                dphidx[i] = (dphidxsi[i]*dydeta - dphideta[i]*dydxsi)/jacobian; // Derivee de phi par rapport a x
                dphidy[i] = (dphideta[i]*dxdxsi - dphidxsi[i]*dxdeta)/jacobian; // Derivee de phi par rapport a y
            }
            for(i = 0; i < theSpace->n; i++){ 
                for(j = 0; j < theSpace->n; j++){ 
                    theSystem->A[map[i]*2][map[j]*2] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jacobian * weight;                                                                                       
                    theSystem->A[map[i]*2][map[j]*2+1] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jacobian * weight;                                                                                 
                    theSystem->A[map[i]*2+1][map[j]*2] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jacobian * weight;    
                    theSystem->A[map[i]*2+1][map[j]*2+1] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jacobian * weight;
                }
                theSystem->B[map[i]*2] += (phi[i] * -rho * g + phi[i] * g)* jacobian *weight;
                theSystem->B[map[i]*2+1] += (phi[i] * -rho * g + phi[i] * g)* jacobian *weight;
            }
        }
        
    }
  
    
                
                
  
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }}
                            
    return femFullSystemEliminate(theSystem);
}
