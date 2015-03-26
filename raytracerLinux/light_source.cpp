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
	
	Colour ka = ray.intersection.mat->ambient;
	Colour kd = ray.intersection.mat->diffuse;
	Colour ks = ray.intersection.mat->specular;
	Vector3D N = ray.intersection.normal;
	Vector3D D = ray.dir;
	Vector3D MS = -(2*D.dot(N) * N) + D;
	Vector3D S = _pos - ray.intersection.point;
	Vector3D M = (2*S.dot(N) * N) - S;
	Vector3D C = -D;
	double s_exp = ray.intersection.mat->specular_exp;
	
	N.normalize();
	D.normalize();
	MS.normalize();
	S.normalize();
	M.normalize();
	C.normalize();
	
	Colour Ia = _col_ambient * ka;
	Colour Id = fmax(0, N.dot(S)) * (_col_diffuse * kd);
	Colour Is = (pow(fmax(0, C.dot(M)), s_exp)) * _col_specular * ks;
	
	ray.col = ray.col + Ia + Id + Is;
	ray.col.clamp();
}

