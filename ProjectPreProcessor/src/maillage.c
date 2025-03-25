#include "fem.h"
#define M_PI		3.14159265358979323846

void chk(int ierr) {
    if (ierr) {
        printf("Error: %d\n", ierr);
        exit(ierr);
    }
}


double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
  
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;


//
//     A modifier !
//     
// Your contribution starts here ....
//
    // double hstar = h;
    // double d = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) - r0;
    // if (d < d0){
    //     double alpha = (-2*h + 2*h0)/(d0*d0*d0);
    //     double beta = (3*h - 3*h0)/(d0*d0);
    //     double gamma = 0;
    //     hstar = alpha*d*d*d + beta*d*d + gamma*d + h0;
    // }
    // d = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1)) - r1;
    // if (d < d1){
    //     double alpha = (-2*h + 2*h1)/(d1*d1*d1);
    //     double beta = (3*h - 3*h1)/(d1*d1);
    //     double gamma = 0;
    //     hstar = fmin(hstar,alpha*d*d*d + beta*d*d + gamma*d + h1);

    // }



    
     
    // return hstar;
    return h;
    
//   
// Your contribution ends here :-)
//

}


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
     
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
 
//
//  -1- Construction de la g�om�trie avec OpenCascade
//      On cr�e le rectangle
//      On cr�e les deux cercles
//      On soustrait les cercles du rectangle :-)
//
 
    int ierr;
    int idPlate = gmshModelOccAddRectangle(-w/2.0,-h/2.0-r0/6.0,0.0,2*w,r0/3,-1,0.0,&ierr); 
    ErrorGmsh(ierr);
    int idNotch_left_wheel = gmshModelOccAddDisk(x0,y0,0.0,r0,r0,-1,NULL,0,NULL,0,&ierr); 
    ErrorGmsh(ierr);


    int idNotch_right_wheel = gmshModelOccAddDisk(x0+2*w+0.65,y0,0.0,r0,r0,-1,NULL,0,NULL,0,&ierr); 
    ErrorGmsh(ierr);
    int idRotatedRectangle = gmshModelOccAddRectangle(-w/2.0+0.30,-h/2.0-r0/6.0-0.5,0.0,1.5*w,r0/3,-1,0.0,&ierr);
    ErrorGmsh(ierr);

    int idFourche = gmshModelOccAddRectangle(-w/2.0-0.90,-h-0.05,0.0,3*w,r0/3,-1,0.0,&ierr);
    ErrorGmsh(ierr);

    int idguidonrec = gmshModelOccAddRectangle(w+0.4,h-0.15,0.0,0.35,r0/4,-1,0.0,&ierr);
    ErrorGmsh(ierr);

    int idguidoncircle = gmshModelOccAddDisk(w+0.4,h-0.15+0.075,0.0,r0/6,r0/6,-1,NULL,0,NULL,0,&ierr);
    ErrorGmsh(ierr);


    int plate[] = {2,idPlate};
    int notch_left_wheel[] = {2,idNotch_left_wheel};
    int notch_right_wheel[] = {2,idNotch_right_wheel};
    int entities[] = {2, idRotatedRectangle};
    int fourche[] = {2, idFourche};
    int guidonrec[] = {2, idguidonrec};
    int guidoncircle[] = {2, idguidoncircle};
    // gmshModelOccCut(notch_right_wheel,2,hole_right_wheel,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    // ErrorGmsh(ierr);
    gmshModelOccRotate(entities, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 55.0 * M_PI / 180.0, &ierr);
    ErrorGmsh(ierr);

    gmshModelOccRotate(fourche, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 100 * M_PI / 180.0, &ierr);
    ErrorGmsh(ierr);

    gmshModelOccCut(plate,2,notch_left_wheel ,2,NULL,NULL,NULL,NULL,NULL,-1,1,0,&ierr); 
    ErrorGmsh(ierr);

    gmshModelOccCut(entities,2,plate ,2,NULL,NULL,NULL,NULL,NULL,-1,1,0,&ierr); 
    ErrorGmsh(ierr);

    gmshModelOccCut(entities,2,fourche ,2,NULL,NULL,NULL,NULL,NULL,-1,1,0,&ierr); 
    ErrorGmsh(ierr);

    gmshModelOccCut(fourche,2,notch_right_wheel ,2,NULL,NULL,NULL,NULL,NULL,-1,1,0,&ierr); 
    ErrorGmsh(ierr);

    gmshModelOccCut(fourche,2,guidonrec ,2,NULL,NULL,NULL,NULL,NULL,-1,1,0,&ierr); 
    ErrorGmsh(ierr);

    gmshModelOccCut(guidonrec,2,guidoncircle ,2,NULL,NULL,NULL,NULL,NULL,-1,1,0,&ierr); 
    ErrorGmsh(ierr);

//
//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
    int object1[] = {2, idPlate, 2, idNotch_left_wheel, 2, idNotch_right_wheel, 2, idRotatedRectangle, 2, idFourche, 2, idguidonrec, 2, idguidoncircle};
    int object2[] = {2, idguidoncircle};
    gmshModelOccFuse(object1, 14, object2, 2, NULL, NULL, NULL, NULL,NULL, -1, 1, 0, &ierr);
    ErrorGmsh(ierr);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
       
//
//  Generation de quads :-)
//
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  chk(ierr);
//    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  chk(ierr);
//    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  chk(ierr);
//    gmshModelMeshGenerate(2, &ierr);  
   
 
//
//  Plot of Fltk
//
  gmshFltkInitialize(&ierr);
  gmshFltkRun(&ierr);  
  chk(ierr);
//
    
}