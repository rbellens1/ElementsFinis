/*
 *  main.c
 *  Projet 2023-2024
 *  Elasticite lineaire plane
 *
 *  Preprocesseur
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat, Miguel De Le Court
 *  All rights reserved.
 *
 */

 #include "glfem.h"

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    // Define your scroll callback functionality here
    printf("Scroll: xoffset = %f, yoffset = %f\n", xoffset, yoffset);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    // Define your callback functionality here

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        printf("Left mouse button pressed\n");
    }
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
        printf("Left mouse button released\n");
    }

}

int main(void) {
 
   //
   //  -1- Construction de la geometrie
   //
   geoInitialize();
   femGeo *theGeometry = geoGetGeometry();
   theGeometry->h = 1.0; // TODO faire varier h 
   geoMeshGenerateGeoFile("../data/mesh.geo");
 
 
   geoMeshImport();
   geoSetDomainName(3, "RouesGauche");
   geoSetDomainName(1, "RouesDroite");
   geoSetDomainName(2, "ChargeHaut");
   geoSetDomainName(0, "LibreBas");
   geoMeshWrite("../data/mesh.txt");
 
   //
   //  -2- Definition du probleme
   //
 
   double E = 50.e9; // 50 GPa carbon fiber
   double nu = 0.28; // Poisson ratio
   double rho = 1.75e3; // 1.75 g/cm3
   double gx = 0; // gravity
   double gy = -9.81; // gravity
 
   femProblem *theProblem = femElasticityCreate(theGeometry, E, nu, rho, gx, gy, PLANAR_STRAIN);
   // Les roues (fixées)
   femElasticityAddBoundaryCondition(theProblem, "RouesGauche", DIRICHLET_XY, 0.0, 0.0);
   femElasticityAddBoundaryCondition(theProblem, "RouesDroite", DIRICHLET_XY, 0.0, 0.0);

   // Charge appliquée sur le haut (ex: Poids du conducteur, 800N)
   femElasticityAddBoundaryCondition(theProblem, "ChargeHaut", NEUMANN_Y, -800.0, NAN);

   // Bord bas libre (aucune contrainte)


   femElasticityPrint(theProblem);
   femElasticityWrite(theProblem, "../data/problem.txt");
 
   //
   //  -3- Champ de la taille de référence du maillage (uniquement pour la visualisation)
   //
 
   double *meshSizeField = malloc(theGeometry->theNodes->nNodes * sizeof(double));
   femNodes *theNodes = theGeometry->theNodes;
   for (int i = 0; i < theNodes->nNodes; ++i)
     meshSizeField[i] = theGeometry->geoSize(theNodes->X[i], theNodes->Y[i]);
   double hMin = femMin(meshSizeField, theNodes->nNodes);
   double hMax = femMax(meshSizeField, theNodes->nNodes);
   printf(" ==== Global requested h : %14.7e \n", theGeometry->h);
   printf(" ==== Minimum h          : %14.7e \n", hMin);
   printf(" ==== Maximum h          : %14.7e \n", hMax);
 
   //
   //  -4- Visualisation
   //
 
   int mode = 1;
   int domain = 0;
   int freezingButton = FALSE;
   double t, told = 0;
   char theMessage[MAXNAME];
 
   GLFWwindow *window = glfemInit("EPL1110 : Project 2024-25 ");
   glfwMakeContextCurrent(window);
   glfwSetScrollCallback(window, scroll_callback);
   glfwSetMouseButtonCallback(window, mouse_button_callback);
 
   do {
     int w, h;
     glfwGetFramebufferSize(window, &w, &h);
     glfemReshapeWindows(theGeometry->theNodes, w, h);
 
     t = glfwGetTime();
     if (glfwGetKey(window, 'D') == GLFW_PRESS) {
       mode = 0;
     }
     if (glfwGetKey(window, 'V') == GLFW_PRESS) {
       mode = 1;
     }
     if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE) {
       domain++;
       freezingButton = TRUE;
       told = t;
     }
 
     if (t - told > 0.5) {
       freezingButton = FALSE;
     }
     if (mode == 1) {
       glfemPlotField(theGeometry->theElements, meshSizeField);
       glfemPlotMesh(theGeometry->theElements);
       sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
       glColor3f(1.0, 0.0, 0.0);
       glfemMessage(theMessage);
     }
     if (mode == 0) {
       domain = domain % theGeometry->nDomains;
       glfemPlotDomain(theGeometry->theDomains[domain]);
       sprintf(theMessage, "%s : %d ", theGeometry->theDomains[domain]->name, domain);
       glColor3f(1.0, 0.0, 0.0);
       glfemMessage(theMessage);
     }
 
     glfwSwapBuffers(window);
     glfwPollEvents();
   } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) != 1);
 
   // Check if the ESC key was pressed or the window was closed
 
   free(meshSizeField);
   femElasticityFree(theProblem);
   geoFree();
   glfwTerminate();
 
   exit(EXIT_SUCCESS);
   return 0;
 }
 
