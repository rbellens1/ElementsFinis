/*
 *  main.c
 *  Projet 2023-2024
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat, Miguel De Le Court
 *  All rights reserved.
 *
 */

 #include "fem.h"

 int main(void) {
   femGeo *theGeometry = geoGetGeometry();
   geoMeshRead("../data/mesh.txt");
   femProblem *theProblem = femElasticityRead(theGeometry, "../data/problem.txt");
   femElasticityPrint(theProblem);
   double *theSoluce = femElasticitySolve(theProblem);
   int nNodes = theGeometry->theNodes->nNodes;
   femSolutionWrite(nNodes, 2, theSoluce, "../data/UV.txt");
   femElasticityFree(theProblem);
   geoFree();
   return 0;
 }