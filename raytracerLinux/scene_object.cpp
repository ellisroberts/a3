/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"
#include <stdio.h>

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	
        
	Point3D R1 = worldToModel * ray.origin;
	Vector3D R = worldToModel * ray.dir;
	
	
	// S1 (upper left vertex), S2 (upper right vertex), and S3 (lower left vertex)
	Point3D S1 = Point3D(-0.5, 0.5, 0);
	Point3D S2 = Point3D(0.5, 0.5, 0);
	Point3D S3 = Point3D(-0.5, -0.5, 0);
	
	Vector3D dS21 = S2 - S1;
	Vector3D dS31 = S3 - S1;
	Vector3D N = Vector3D(0, 0, 1);
	
	// Compute the normal in the world
	ray.intersection.normal = transNorm(modelToWorld, N);
        ray.intersection.normal.normalize();
	
	// Check if parallel
	if(N.dot(R) < 1e-25f){
		//ray.intersection.none = true;
		return false;
	}

	// Compute t value
	double t = -N.dot(R1 - S1) / N.dot(R);
	printf("t is %f\n", t);
	
	// Check if the origin is behind the square
	if(t < 0){
		//ray.intersection.none = true;
		return false;
	}
	
	// Let M be the intersection point between the ray and the surface
	Point3D M = R1 + t * R;
	ray.intersection.point = modelToWorld * M;
	
	// Get the vector between M and S1
	Vector3D dMS1 = M - S1;
	
	// Compute the position of M relative to dS21 and dS31
	double u = dMS1.dot(dS21);
        double v = dMS1.dot(dS31);

	

	// Check if the intersection point belongs to the square
	if((u >= 0.0) && (u <= dS21.dot(dS21)) && (v >= 0.0) && (v <= dS31.dot(dS31))){
                printf("squsre bitch\n");
                if (ray.intersection.none || (t < ray.intersection.t_value))
                {
                    ray.intersection.t_value = t;
		    ray.intersection.none = false;
                }
		return true;
	}else{
		//ray.intersection.none = true;
		return false;
	}
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

        //convert to model coordinates
        Point3D origin = worldToModel*ray.origin;
        Vector3D dir = worldToModel*ray.dir;

        Point3D center(0,0,0);
        Vector3D difference = origin - center;
        
        //Get the lamba value of the intersection with the sphere
        double lambda;
        double A = dir.dot(dir);
        double B = 2*dir.dot(difference);
        double C = difference.dot(difference) - 1;	
        double D = B*B-A*C;
 
        //printf("origin is %f %f %f", origin[0], origin[1], origin[2]);
        //printf("dir is %f %f %f", dir[0], dir[1], dir[2]);
        //printf("point is A is %f B is %f C is %f D is %f\n", A, B, C, D);

        //not intersection
        if (D < 0.0)
        {
            return false;
        }

        double lambda1 = -B/A + sqrt(D)/A;
        double lambda2 = -B/A - sqrt(D)/A;

        //there will be two lambda values. choose the smaller one.
        if (lambda1 >= lambda2)
        {
            lambda = lambda2;   
        }
        else
        {
            lambda = lambda1;
        }  
 
        if (ray.intersection.none || (lambda < ray.intersection.t_value))
        {
            ray.intersection.t_value = lambda;
            ray.intersection.none = false;
            //the intersection point
            ray.intersection.point = ray.origin + ray.intersection.t_value*(ray.dir);
            //the normal is the gradient
            Vector3D normal(ray.intersection.point[0], ray.intersection.point[1], ray.intersection.point[2]);
            //normalize the vectorest 
            ray.intersection.normal = normal;
            ray.intersection.point = modelToWorld*ray.intersection.point;
            ray.intersection.normal = modelToWorld.transpose()*ray.intersection.normal;
            ray.intersection.normal.normalize();
            return true;
            
        }



        
 
	return false;
}

