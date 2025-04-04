/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Elasticite lineaire plane
 *  Calcul des densités de force aux noeuds contraints
 *
 *  Copyright (C) 2024 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "glfem.h"

double fun(double x, double y) 
{
    return 1;
}






int read_problem(double *E, double *nu, double *rho, double *g, double *mass, femElasticCase *type){

    
    FILE *f = fopen("../data/problem.txt", "r");
    if (f == NULL) {
        fprintf(stderr, "Error: cannot open file problem.txt\n");
        return 1;
    }
    char problemType[256];
    // Lecture du type de problème
    if (fscanf(f, "Type of problem : %[^\n]\n", problemType) != 1) {
        fprintf(stderr, "Error: cannot read problem type\n");
        return 1;
    }

    printf("prob type |%s|\n", problemType);

    // Conversion du texte en type femElasticCase
    if (strcmp(problemType, "planar stress") == 0) {
        *type = PLANAR_STRESS;
    } else if (strcmp(problemType, "planar strain") == 0) {
        *type = PLANAR_STRAIN;
    } else if (strcmp(problemType, "axisymetrie") == 0){
        *type = AXISYM;
    }else {
        fprintf(stderr, "Error: unknown problem type\n");
        return 1;
    }

    // Lecture des paramètres matériaux
    if (fscanf(f, "Young modulus : %lf\n", E) != 1) {
        fprintf(stderr, "Error: cannot read Young modulus\n");
        return 1;
    }

    if (fscanf(f, "Poisson ratio : %lf\n", nu) != 1) {
        fprintf(stderr, "Error: cannot read Poisson ratio\n");
        return 1;
    }

    if (fscanf(f, "Mass density : %lf\n", rho) != 1) {
        fprintf(stderr, "Error: cannot read density\n");
        return 1;
    }

    if (fscanf(f, "Gravity : %lf\n", g) != 1) {
        fprintf(stderr, "Error: cannot read gravity\n");
        return 1;
    }

    if (fscanf(f, "Mass : %lf\n", mass) != 1) {
        fprintf(stderr, "Error: cannot read mass\n");
        return 1;
    }

    fclose(f);
    return 0;



    return 0;
}











int main(void)
{  
    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    X : Horizontal residuals for unconstrained equations \n");
    printf("    Y : Horizontal residuals for unconstrained equations \n");
    printf("    N : Next domain highlighted\n\n\n");
      
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    
    theGeometry->elementType = FEM_TRIANGLE;

    geoMeshRead("../data/mesh.txt");
        
//
//  -2- Creation probleme 
//
    double E = 0.0;
    double nu  = 0.0;
    double rho = 0.0;
    double g   = 0.0;
    double mass = 0.0;

    femElasticCase type = PLANAR_STRESS;


    if(read_problem(&E, &nu, &rho, &g, &mass, &type)) return EXIT_FAILURE;



    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,type);

    femElasticityAddBoundaryCondition(theProblem,"Entity 52 ",NEUMANN_Y,-mass*g);

    femElasticityAddBoundaryCondition(theProblem,"DirichletConstrainRight",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"DirichletConstrainRight",DIRICHLET_Y,0.0);

    femElasticityAddBoundaryCondition(theProblem,"DirichletConstrainLeft",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"DirichletConstrainLeft",DIRICHLET_Y,0.0);

    femElasticityPrint(theProblem);

//
//  -3- Resolution du probleme et calcul des forces
//

    double *theSoluce = femElasticitySolve(theProblem);
    double *theForces = femElasticityForces(theProblem);
    double area = femElasticityIntegrate(theProblem, fun);   
   
//
//  -4- Deformation du maillage pour le plot final
//      Creation du champ de la norme du deplacement
//
    
    femNodes *theNodes = theGeometry->theNodes;
    double deformationFactor = 1e5;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    double *forcesX = malloc(theNodes->nNodes * sizeof(double));
    double *forcesY = malloc(theNodes->nNodes * sizeof(double));
    
    for (int i=0; i<theNodes->nNodes; i++){
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
        theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                                   theSoluce[2*i+1]*theSoluce[2*i+1]);
        forcesX[i] = theForces[2*i+0];
        forcesY[i] = theForces[2*i+1]; }
  
    double hMin = femMin(normDisplacement,theNodes->nNodes);  
    double hMax = femMax(normDisplacement,theNodes->nNodes);  
    printf(" ==== Minimum displacement          : %14.7e [m] \n",hMin);
    printf(" ==== Maximum displacement          : %14.7e [m] \n",hMax);

//
//  -5- Calcul de la force globaleresultante
//

    double theGlobalForce[2] = {0, 0};
    for (int i=0; i<theProblem->geometry->theNodes->nNodes; i++) {
        theGlobalForce[0] += theForces[2*i+0];
        theGlobalForce[1] += theForces[2*i+1]; }
    printf(" ==== Global horizontal force       : %14.7e [N] \n",theGlobalForce[0]);
    printf(" ==== Global vertical force         : %14.7e [N] \n",theGlobalForce[1]);
    printf(" ==== Weight                        : %14.7e [N] \n", area * rho * g);





    // Write out the solution
    int nNodes = theGeometry->theNodes->nNodes;
    femSolutionWrite(nNodes, 2, theSoluce, "../data/outfile.txt");

//
//  -6- Visualisation du maillage
//  
    
    int mode = 1; 
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];
   
 
    GLFWwindow* window = glfemInit("EPL1110 : Recovering forces on constrained nodes");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'X') == GLFW_PRESS) { mode = 2;}
        if (glfwGetKey(window,'Y') == GLFW_PRESS) { mode = 3;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        if (t-told > 0.5) {freezingButton = FALSE; }
        
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements,normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 2) {
            glfemPlotField(theGeometry->theElements,forcesX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 3) {
            glfemPlotField(theGeometry->theElements,forcesY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
         glfwSwapBuffers(window);
         glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    free(normDisplacement);
    free(forcesX);
    free(forcesY);
    femElasticityFree(theProblem) ; 
    geoFinalize();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}

 
