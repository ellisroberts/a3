/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include <math.h>
#include <iostream>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  
	
	Colour ambient = ray.intersection.mat->ambient;
	Colour diffuse = ray.intersection.mat->diffuse;
	Colour specular = ray.intersection.mat->specular;
	double specular_exp = ray.intersection.mat->specular_exp;
	Vector3D N = ray.intersection.normal;
	Vector3D R = ray.dir;
	N.normalize();
	R.normalize();
	
	//std::cout << "mat:" << ambient << diffuse << specular << "\n";
	
	
	// Direction of the perfect reflected ray
	Vector3D S = R - 2*(R.dot(N)/N.dot(N)) * N;
	
	Colour Ia = ambient;
	Colour Id = (N.dot(R))*diffuse;
	Colour Is = (pow((N.dot(S)), specular_exp))*specular;
	
	ray.col = ray.col + Ia + Id + Is;
	
	//std::cout << "mat:" << N << " " << R << "\n";
	//std::cout << ray.col << "\n";
}

