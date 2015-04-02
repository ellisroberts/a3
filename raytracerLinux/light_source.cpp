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
#include <stdlib.h>


std::vector<Ray3D> PointLight::get_shadow_rays(Ray3D& ray){
	std::vector<Ray3D> shadow_rays;
	
	Point3D R1 = ray.intersection.point - 0.0000001 * ray.dir;
	Vector3D R = get_position() - R1;
	R.normalize();
	Ray3D shadowRay(R1, R);
	
	shadow_rays.push_back(shadowRay);
	
	return shadow_rays;
}

std::vector<Ray3D> ParallelogramLight::get_shadow_rays(Ray3D& ray){
	int i;
	Vector3D R;
	Point3D P;
	std::vector<Ray3D> shadow_rays;
	std::vector<int> r;
	std::vector<int> s;	
	Vector3D x_edge = _p;
	Vector3D y_edge = _q;
	Point3D corner = _pos;
	
	double N = 12; // Number of sample
	
	// Generate Nˆ2 jittered points
	r.clear();
	s.clear();
	for(i=0; i < N; i++){
		r.push_back(((double) rand() / (RAND_MAX)) + 1);
		s.push_back(((double) rand() / (RAND_MAX)) + 1);
	}
	
	// Shuffle
	int aux;
	for(i = N-1; i >=0; i--){
		int j = rand() % (i+1);
		aux = s.at(j);
		s.at(j) = s.at(i);	
		s.at(i) = aux;
	}
	
	// Sample N ray shadows
	Point3D R1 = ray.intersection.point - 0.0000001 * ray.dir;
	for(i = 0; i < N; i++){
		
		P = corner + r.at(i) * x_edge + s.at(i) * y_edge;
		R =  P - R1;
		R.normalize();
		Ray3D shadowRay(R1, R);
	
		shadow_rays.push_back(shadowRay);
	}
	return shadow_rays;
}

void LightSource::shade( Ray3D& ray, double contribution) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  
	 
	Colour ka = ray.intersection.mat->ambient;
	Colour kd = ray.intersection.mat->diffuse;
        //we have a texture
        if (ray.intersection.isMap)
        {
            
            int row, col;
            //get texel coords
            SphereMapping sphereMap;
            sphereMap.getTexel(&row, &col, ray.intersection.text.getHeight(), ray.intersection.text.getWidth(), ray.intersection.pointObjectCoords);
            //change diffuse coefficient to that of texture
            unsigned char * red = ray.intersection.text.getRed(row,col);
            unsigned char * green = ray.intersection.text.getGreen(row,col);
            unsigned char * blue = ray.intersection.text.getBlue(row,col);
            //get values for the color components in proper type
            double RedValue = (*red)/255.0;
            double BlueValue = (*blue)/255.0;
            double GreenValue = (*green)/255.0;
            kd = kd*Colour(RedValue, GreenValue, BlueValue); 
            
        }
	Colour ks = ray.intersection.mat->specular;
	Vector3D N = ray.intersection.normal;
	Vector3D D = ray.dir;
	Vector3D MS = -(2*D.dot(N) * N) + D;
	Vector3D S = get_position() - ray.intersection.point;
	Vector3D M = (2*S.dot(N) * N) - S;
	Vector3D C = -D;
	double s_exp = ray.intersection.mat->specular_exp;
	
	N.normalize();
	D.normalize();
	MS.normalize();
	S.normalize();
	M.normalize();
	C.normalize();
	
	Colour Ia = get_ambient() * ka;
	Colour Id = fmax(0, N.dot(S)) * (get_diffuse() * kd);
	Colour Is = (pow(fmax(0, C.dot(M)), s_exp)) * get_specular() * ks;
	
	ray.col = ray.col + contribution * (Ia + Id + Is);
	ray.col.clamp();
}

