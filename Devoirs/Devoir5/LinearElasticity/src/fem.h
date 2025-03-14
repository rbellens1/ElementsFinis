
/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gmshc.h" 


#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define ErrorGmsh(a)   femErrorGmsh(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1
#define MAXNAME 256

typedef enum {FEM_TRIANGLE,FEM_QUAD} femElementType;
typedef enum {DIRICHLET_X,DIRICHLET_Y} femBoundaryType;
typedef enum {PLANAR_STRESS,PLANAR_STRAIN,AXISYM} femElasticCase;


typedef struct {
    int nNodes; // Nombre de noeuds
    double *X; // Coordonnees X des noeuds
    double *Y; // Coordonnees Y des noeuds
} femNodes;

typedef struct {
    int nLocalNode; // Nombre de noeuds locaux par element
    int nElem; // Nombre d'elements dans le maillage
    int *elem; // Liste des noeuds de chaque element du maillage
    femNodes *nodes; // Noeuds du maillage
} femMesh;

typedef struct {
    femMesh *mesh; // Maillage
    int nElem; // Nombre d'elements dans le domaine
    int *elem; // Liste des elements du domaine
    char name[MAXNAME]; // Nom du domaine
} femDomain;

typedef struct {
    double LxPlate, LyPlate;
    double h;
    femElementType elementType;
    double (*geoSize)(double x, double y);
    femNodes *theNodes;
    femMesh  *theElements;
    femMesh  *theEdges;
    int nDomains;
    femDomain **theDomains;
} femGeo;

typedef struct {
    int n; // Nombre de points d'integration
    void (*x2)(double *xsi, double *eta); // Coordonnees xsi et eta des points d'integration
    void (*phi2)(double xsi, double eta, double *phi); // Fonctions de forme aux points d'integration
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta); // Derivees des fonctions de forme aux points d'integration
} femDiscrete;
    
typedef struct {
    int n; // Nombre de points d'integration
    const double *xsi; // Coordonnees xsi des points d'integration
    const double *eta; // Coordonnees eta des points d'integration
    const double *weight; // Poids des points d'integration
} femIntegration;

typedef struct {
    double *B;
    double **A;
    int size;
} femFullSystem;


typedef struct {
    femDomain* domain;
    femBoundaryType type; 
    double value;
} femBoundaryCondition;


typedef struct {
    double E,nu,rho,g;
    double A,B,C; 
    int planarStrainStress; 
    int nBoundaryConditions; 
    femBoundaryCondition **conditions;  
    int *constrainedNodes; 
    femGeo *geometry;
    femDiscrete *space;
    femIntegration *rule;
    femFullSystem *system;
} femProblem;


void                geoInitialize();
femGeo*             geoGetGeometry();
double              geoSize(double x, double y);
double              geoSizeDefault(double x, double y);
void                geoSetSizeCallback(double (*geoSize)(double x, double y));
void                geoMeshGenerate();
void                geoMeshImport();
void                geoMeshPrint();
void                geoMeshWrite(const char *filename);
void                geoMeshRead(const char *filename);
void                geoSetDomainName(int iDomain, char *name);
int                 geoGetDomain(char *name);
void                geoFinalize();

femProblem*         femElasticityCreate(femGeo* theGeometry, 
                                      double E, double nu, double rho, double g, femElasticCase iCase);
void                femElasticityFree(femProblem *theProblem);
void                femElasticityPrint(femProblem *theProblem);
void                femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value);
double*             femElasticitySolve(femProblem *theProblem);

femIntegration*     femIntegrationCreate(int n, femElementType type);
void                femIntegrationFree(femIntegration *theRule);

femDiscrete*        femDiscreteCreate(int n, femElementType type);
void                femDiscreteFree(femDiscrete* mySpace);
void                femDiscretePrint(femDiscrete* mySpace);
void                femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);

femFullSystem*      femFullSystemCreate(int size);
void                femFullSystemFree(femFullSystem* mySystem);
void                femFullSystemPrint(femFullSystem* mySystem);
void                femFullSystemInit(femFullSystem* mySystem);
void                femFullSystemAlloc(femFullSystem* mySystem, int size);
double*             femFullSystemEliminate(femFullSystem* mySystem);
void                femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value);

double              femMin(double *x, int n);
double              femMax(double *x, int n);
void                femError(char *text, int line, char *file);
void                femErrorScan(int test, int line, char *file);
void                femErrorGmsh(int test, int line, char *file);
void                femWarning(char *text, int line, char *file);


#endif
