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
	Point3D center = Point3D(0, 0, 0);
	
	Vector3D dS21 = S2 - S1;
	Vector3D dS31 = S3 - S1;
	
	Vector3D N = dS21.cross(dS31);
	
	// Normalize
	N.normalize();
	//R.normalize();
	
	// Check if parallel
	if(fabs(N.dot(R)) < 1e-25f){
		return false;
	}
	
	// Compute t value
	Vector3D dCR1 = center - R1;
	double t = dCR1.dot(N)/ N.dot(R);
	
	if(ray.intersection.t_value < t && !ray.intersection.none){
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
	if(u >= 0.0 && u <= dS21.dot(dS21) && v >= 0.0 && v <= dS31.dot(dS31) && t >= 0){
	
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

        // No intersection
        if (D < 0.0)
        {
            return false;
        }
		else{
			double t1 = -B/2*A + sqrt(D)/2*A;
        	double t2 = -B/2*A - sqrt(D)/2*A;
			
			//there will be two lambda values. choose the smaller one.
        	if(fmax(t1, t2) < 0){
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
			
			if(ray.intersection.t_value < t && !ray.intersection.none){
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
                        //we have a texture map
                        ray.intersection.isMap = true;
                        ray.intersection.pointObjectCoords = M;
                        //get the texture object
                        ray.intersection.text = getTexture();
			ray.intersection.t_value = t;
			ray.intersection.none = false;
			return true;
		}
}


bool Cylinder::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {

    Point3D origin = worldToModel*ray.origin;
    Vector3D dir = worldToModel*ray.dir;
		
    dir.normalize();

    Point3D pc(0,0,0);

    Vector3D originVector = origin - pc;
    
    
    //want to intersect p = a + lambda*d with the cylinder x**2 + y**2 = 1
    
    //lambda for cylinder
    double lambda_final1= 1000;
    //lambda for plane
    double lambda_final2 = 1000;

    //to solve for equation of cylinder
    double A = dir[0]*dir[0] + dir[1]*dir[1];
    double B = origin[0]*dir[0] + origin[1]*dir[1];
    double C = origin[0]*origin[0] + origin[1]*origin[1] - 1;

    double D = B*B - A*C;

    //no intersection 
    if (D < 0.0)
    {
        return false;
    }

    if (D > 0.0){
 
        double R = sqrt(D);
        //for solving intersection with quadric
        double lambda1 = -B/A + sqrt(D)/A;
        double lambda2 = -B/A - sqrt(D)/A;

 
	//there will be two lambda values. choose the smaller one.
	double lambdaCylinder = lambda1;
        if ((lambda1 < 0)&&(lambda2 > 0))
        {
            lambdaCylinder = lambda2;
        }

        if ((lambda2 < 0)&&(lambda1 > 0))
        {
            lambdaCylinder = lambda1;
        }
        if ((lambda2 > 0)&&(lambda1 > 0))
        {
            if (lambda2 < lambda1)
            {
                lambdaCylinder = lambda2;
            }
        }

        if ((lambda2 < 0)&&(lambda1 < 0))
        {
            return false;
        }
        
	//the intersection point for the side
        Point3D intersectionCylinder = origin + lambdaCylinder*dir;
        Vector3D NormalCylinder = intersectionCylinder - pc;
        NormalCylinder[2] = 0; 
	NormalCylinder.normalize();
			
	// Update intersection if the point is within the bounds of the cylinder
        if ((intersectionCylinder[2] <= 1) & (intersectionCylinder[2] >= 0.0))
        {
            
            lambda_final1 = lambdaCylinder;
        }
    
	
        //first find the intersection the plane z = 1
        //define the plane
	
	Vector3D NormalPlane = Vector3D(0,0,1);
	
        //Normalize
	NormalPlane.normalize();
        
        //Intersection of the plane
        Point3D intersectionPlane;
	
	// If ray is not paraellel to plane, compute intersection
	if (!(fabs(NormalPlane.dot(dir)) < 1e-25f)){

	    // Compute lambda value for plane at z = 0
	    double lambdaPlane1  = (-origin[2])/dir[2];

            //compute lambdaValue for plane at z = 1
            double lambdaPlane2 = ((1-origin[2])/dir[2]);

            //printf("lambda is %f\n", lambdaPlane);
            //the intersections with both of the planes
            Point3D intersectionPlane1 = origin + lambdaPlane1*dir;
            Point3D intersectionPlane2 = origin + lambdaPlane2*dir;
        
            //we're on the tip of the cylinder and it's the closer intersection
            if (((intersectionPlane2[0]*intersectionPlane2[0] + intersectionPlane2[1]*intersectionPlane2[1]) <= 1))
            {               
                intersectionPlane = intersectionPlane2;        
                lambda_final2 = lambdaPlane2;
                
            }

            //check if we're in the bounds of the cylinder

	    if (((intersectionPlane1[0]*intersectionPlane1[0] + intersectionPlane1[1]*intersectionPlane1[1]) <= 1))
	    {    
              
                //take the lesser of the 2 plane intersections                 
                if (lambdaPlane2 < lambda_final2)
                {       
                    intersectionPlane = intersectionPlane2; 
                    lambda_final2 = lambdaPlane1;
                    //The orientation of the normal changes
                    NormalPlane = -1*NormalPlane;
                }
                
            }
	}

        if ((lambda_final1 == 1000)&&(lambda_final2 == 10000))
        {
            return false; 
        }
	
        //choose the lesser of lambda_final1 and lambda_final2
        if (lambda_final1 < lambda_final2)
        {

	    if(ray.intersection.t_value < lambda_final1 && !ray.intersection.none){
		return false;
	    }
            //we choose the side of the cylinder intersection
	    ray.intersection.normal = modelToWorld * NormalCylinder;
	    ray.intersection.normal.normalize();
	    ray.intersection.point = modelToWorld * intersectionCylinder;
	    ray.intersection.t_value = lambda_final1;
	    ray.intersection.none = false;
	    return true;       
        }

        if (lambda_final2 < lambda_final1)
        {

	    if(ray.intersection.t_value < lambda_final2 && !ray.intersection.none){
		return false;
	    }
            //we choose the side of the cylinder intersection
	    ray.intersection.normal = modelToWorld * NormalPlane;
	    ray.intersection.normal.normalize();
	    ray.intersection.point = modelToWorld * intersectionPlane;
	    ray.intersection.t_value = lambda_final1;
	    ray.intersection.none = false;
	    return true;       
        }
   }

   return false;

}
