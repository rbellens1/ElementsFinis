# Pré-traitement

## Auteurs

- **Etudiant:** Gabriel Zaoudi
- **Etudiant:** Romain Bellens

## Informations générales

- **Cours:** LEPL1110
- **Institution:** UCLouvain

## Description

Ce dossier contient le préprocesseur développé pour le cours LEPL1110 à l'UCLouvain axé sur les éléments finis. Le préprocesseur a pour but de générer le maillage (sauvegardé dans un fichier `mesh.txt`) et de définir le problème (sauvegardé dans un fichier `problem.txt`). //TODO problem.txt

## Prérequis

Le préprocesseur repose sur les bibliothèques suivantes :

- [CMake](https://cmake.org/)
- [Gmsh](https://gmsh.info/)
- [OpenGl](https://www.opengl.org/)
- [GLFW](https://www.glfw.org/)

## Instructions de compilation

Pour compiler le préprocesseur, exécutez les commandes suivantes dans votre terminal :

```bash
mkdir build
cd build
cmake ..
make
```

## Utilisation

Pour exécuter le préprocesseur, exécutez la commande suivante dans votre terminal :

```bash
./myFem
```


