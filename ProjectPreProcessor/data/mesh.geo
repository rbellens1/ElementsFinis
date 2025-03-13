h = 1.0; // Taille caractéristique du maillage

// Définition des points
Point(1) = {0, 0, 0, h};
Point(2) = {10, 0, 0, h};
Point(3) = {10, 5, 0, h};
Point(4) = {0, 5, 0, h};

// Définition des lignes
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Création du contour
Line Loop(1) = {1, 2, 3, 4};

// Création de la surface
Plane Surface(1) = {1};


// Groupes physiques pour les frontières
Physical Line("RouesGauche")  = {4};  // Fixe
Physical Line("RouesDroite")  = {2};  // Fixe
Physical Line("ChargeHaut")   = {3};  // Charge appliquée
Physical Line("LibreBas")     = {1};  // Bord libre

