/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

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
	Point3D S1 = Point3D(0.5, 0.5, 0);
	Point3D S2 = Point3D(-0.5, 0.5, 0);
	Point3D S3 = Point3D(0.5, -0.5, 0);
	
	Vector3D dS21 = S2 - S1;
	Vector3D dS31 = S3 - S1;
	
	Vector3D N = dS21.cross(dS31);
	
	// Normalize
	R.normalize();
	N.normalize();
	
	// Check if parallel
	if(fabs(N.dot(R)) < 1e-25f){
		return false;
	}
	
	// Compute t value
	Vector3D dR1S1 = R1 - S1;
	double t = -N.dot(dR1S1) / N.dot(R);
	
	if(ray.intersection.t_value > t && !ray.intersection.none){
		return false;
	}
	
	// Let M be the intersection point between the ray and the surface
	Point3D M = R1 + t * R;
	
	// Get the vector between M and S1
	Vector3D dMS1 = M - S1;
	
	// Compute the position of M relative to dS21 and dS31
	double u = dMS1.dot(dS21);
    double v = dMS1.dot(dS31);
	
	// Check if the intersection point belongs to the square
	if(u >= 0.0 && u <= dS21.dot(dS21) && v >= 0.0 && v <= dS31.dot(dS31) && t > 0){
	
		// Update intersection
		ray.intersection.normal = modelToWorld * N;
		ray.intersection.normal.normalize();	
		ray.intersection.t_value = t;
		ray.intersection.point = modelToWorld * M;
		ray.intersection.none = false;
		return true;
	}else{
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
        Point3D R1 = worldToModel*ray.origin;
        Vector3D R = worldToModel*ray.dir;
		
		R.normalize();

        Point3D pc(0,0,0);
	
        //Get the lambda value of the intersection with the sphere
        double t;
        double A = R.dot(R);
        double B = 2*R.dot(R1-pc);
        double C = (R1 - pc).dot(R1 - pc) - 1;	
        double D = B*B-4*A*C;
 
        //printf("origin is %f %f %f", origin[0], origin[1], origin[2]);
        //printf("dir is %f %f %f", dir[0], dir[1], dir[2]);
        //printf("point is A is %f B is %f C is %f D is %f\n", A, B, C, D);

        // No intersection
        if (D < 0.0)
        {
            return false;
        }
		else{
			double t1 = -B/2*A + sqrt(D)/2*A;
        	double t2 = -B/2*A - sqrt(D)/2*A;
			
			//there will be two lambda values. choose the smaller one.
        	if(fmax(t1, t2) < 0 || t1 == t2){
				// Intersect in wrong direction
				return false;
			}
			else{
				// Get the minimum non negative element
				if(fmin(t1, t2) > 0){
					t = fmin(t1, t2);
				}
				else{
					t = fmax(t1, t2);
				}
			}
			
			if(ray.intersection.t_value > t && !ray.intersection.none){
				return false;
			}
	
			//the intersection point
            Point3D M = R1 + t*R;
            Vector3D N = M - pc;
			N.normalize();
			
			// Update intersection
			ray.intersection.normal = modelToWorld * N;
			ray.intersection.normal.normalize();
			ray.intersection.point = modelToWorld * M;
			ray.intersection.t_value = t;
			ray.intersection.none = false;
			return true;
		}
}