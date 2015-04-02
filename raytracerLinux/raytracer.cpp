/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return node;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation; 	
			angle = -angle;
		} 
		else {
			node->invtrans = rotation*node->invtrans; 
		}	
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
		Vector3D up ) {
	Matrix4x4 mat; 
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat; 
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
	SceneDagNode *childPtr;

	// Applies transformation of the current node to the global
	// transformation matrices.
	_modelToWorld = _modelToWorld*node->trans;
	_worldToModel = node->invtrans*_worldToModel;
	if (node->obj) {
	
		// Perform intersection.
		if (node->obj->intersect(ray, _worldToModel, _modelToWorld)) {
			ray.intersection.mat = node->mat;
		}
	}
	
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		traverseScene(childPtr, ray);
		childPtr = childPtr->next;
	}
	
	// Removes transformation of the current node from the global
	// transformation matrices.
	_worldToModel = node->trans*_worldToModel;
	_modelToWorld = _modelToWorld*node->invtrans;
	
}

void Raytracer::computeShading( Ray3D& ray ) {
	LightListNode* curLight = _lightSource;
	int i, N;
	std::vector<Ray3D> shadow_rays;
	Ray3D shadowRay;
	for (;;) {
		if (curLight == NULL) break;
		// Each lightSource provides its own shading function.
		// Shadows
		
		shadow_rays = curLight->light->get_shadow_rays(ray);
		N = shadow_rays.size();
		
		for(i=0; i < N; i++){
			shadowRay = shadow_rays.at(i);
			
			// Check intesection
			traverseScene(_root, shadowRay);
			if(shadowRay.intersection.none){
				curLight->light->shade(ray, 1.0/N);
			}else{
				Colour ka = ray.intersection.mat->ambient;
				Colour Ia = curLight->light->get_ambient() * ka;
				ray.col = ray.col + (1.0/N) * Ia * ka;
				ray.col.clamp();
			}
		}
		//std::cout << ray.col << "\n";
		curLight = curLight->next;
	}
}

void Raytracer::initPixelBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::flushPixelBuffer( char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}

Colour Raytracer::shadeRay( Ray3D& ray ) {
	Colour curRayCol(0.0, 0.0, 0.0); 
	Colour reflRayCol(0.0, 0.0, 0.0);
	Colour refrRayCol(0.0, 0.0, 0.0);
	
	traverseScene(_root, ray);
	
	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none && (ray.maxDepth < 3)) {
		computeShading(ray);
		curRayCol = ray.col;
		
		// Spawn Reflection
		Vector3D N = ray.intersection.normal;
		Vector3D V = ray.dir;
		
		N.normalize();
		V.normalize();
		
		// Stochastic Reflection (Glossy Effect)
		double a = 0.01; // Blur degree
		
		double e1 = ((double) rand() / (RAND_MAX)) + 1;
		double e2 = ((double) rand() / (RAND_MAX)) + 1;
		
		double x = -a/2 + e1 * a;
		double y = -a/2 + e2 * a;
		
		Vector3D R = V - 2*(N.dot(V))*N;
		R.normalize();
		
		Vector3D X(1, 0, 0);
		Vector3D Y(0, 1, 0);
		
		Point3D R_perturbed_origin = ray.intersection.point - 0.0000001 * ray.dir;
		Vector3D R_perturbed = R + x*X + y*Y;
		R_perturbed.normalize();
		
		// Mirror-Like Reflection		
		//Vector3D R = V - 2*(N.dot(V))*N;
		//Point3D R1 = ray.intersection.point - 0.0000001 * ray.dir;
		
		// Define reflected ray
		Ray3D reflRay(R_perturbed_origin, R_perturbed);
		
		// Check if the ray is below the surface
		if(N.dot(R_perturbed) > 0){
			
			// Ray trace reflected ray 
			reflRay.maxDepth = ray.maxDepth + 1;
			reflRayCol = pow(ray.intersection.mat->reflection_factor, ray.maxDepth) * shadeRay(reflRay);
		}
		else{
			reflRayCol = Colour(0, 0, 0);
		}
		
		// Spawn Refraction
		double n1;
		double n2;
		if(ray.maxDepth%2 == 0){
		
			// Outside/Inside
			n1 = 1;
			n2 = ray.intersection.mat->ref_idx;

		}else{
		
			// Inside/Outside
			n1 = ray.intersection.mat->ref_idx;
			n2 = 1;
		}
		
		double n = n1/n2;
		double k;
		double c1 = -V.dot(N);
		double D = 1-n * n * (1 - c1 * c1);
		double c2 = sqrt(D);
		Vector3D Rr;
		
		// Refraction
		if(D >= 0){
			if(c1 >= 0){
				k = n * c1 - c2;
			}
			else{
				k = n * c1 + c2;
			}
			
			Rr = k * N - n * V;
			
		}
		
		// Total reflection
		else{
			Rr = 2*c1 * N - V;
		}
		Rr.normalize();
		Point3D Rr1 = ray.intersection.point + 0.0000001 * ray.dir;
		
		Ray3D refrRay(Rr1, Rr);
		refrRay.maxDepth = ray.maxDepth + 1;
		
		refrRayCol = pow(0.5, ray.maxDepth) * shadeRay(refrRay);
	}
	curRayCol = curRayCol + reflRayCol + refrRayCol;
	curRayCol = curRayCol;
	curRayCol.clamp();
	return curRayCol; 
}	

void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
		Vector3D up, double fov, char* fileName ) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	double factor = (double(height)/2)/tan(fov*M_PI/360.0);
		
	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);
	
	//number of points being sampled
	int sample = 16;
	
	// Construct a ray for each pixel.
   for (int i = 0; i < _scrHeight; i++) {
 		for (int j = 0; j < _scrWidth; j++) {
			Colour col;
   			for (float pixeli = i; pixeli < i + 1; pixeli = pixeli + .25)  {
				for (float pixelj = j; pixelj < j + 1; pixelj += .25) {
		
					// Sets up ray origin and direction in view space, 
					// image plane is at z = -1.
					Point3D origin(0, 0, 0);
					Point3D imagePlane;
					imagePlane[0] = (-double(width)/2 + 0.5 + pixelj)/factor;
					imagePlane[1] = (-double(height)/2 + 0.5 + pixeli)/factor;
					imagePlane[2] = -1;
			
					// shadeRay(ray) to generate pixel colour. 	
					//position of pixel in world coordinates
					Point3D pixelPointWorld;
			
					//direction of pixel relative to camera origin
					Vector3D pixelDirection;
            
					//direction of the casted ray
					Vector3D direction;
				
					pixelPointWorld = viewToWorld*imagePlane;
					pixelDirection = pixelPointWorld-eye;
					Ray3D ray(eye, pixelDirection);
					ray.maxDepth = 0;
					ray.dir.normalize();
			
					//sum the colours from the different ray components
					col = col + shadeRay(ray); 
				}
			}	
				//store the result of numerical integration into the pixel
			_rbuffer[i*width+j] = int(col[0]*255/sample);
			_gbuffer[i*width+j] = int(col[1]*255/sample);
			_bbuffer[i*width+j] = int(col[2]*255/sample);
		}
	}
	flushPixelBuffer(fileName);
}

int main(int argc, char* argv[])
{	
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;
	int width = 640/4; 
	int height = 480/4; 

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
			Colour(0.628281, 0.555802, 0.366065), 
			51.2 , 0.27049, 0.3);
	
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8 , 1.627, 0.25);
			
	Material metal( Colour(0.1, 0.1, 0.1), Colour(0.6, 0.6, 0.6), 
			Colour(0.7, 0.7, 0.3), 
			51.2 , 1.627, 0.8);
			
	Material glass( Colour(0.0, 0.0, 0.0), Colour(0.588235, 0.670588, 0.729412), 
			Colour(0.9, 0.9, 0.9), 
			96 , 1.5, 0.8);
			
	Material ruby( Colour(0.1745, 0.01175, 0.01175), Colour(0.61424, 0.04136, 0.04136), 
			Colour(0.727811, 0.626959, 0.626959), 
			76.8 , 1.5, 0.3);
	
	Material obsedian( Colour(0.05375, 0.05, 0.06625), Colour(0.18275, 0.17, 0.22525), 
			Colour(0.332741, 0.328634, 0.346435), 
			38.4 , 1.5, 0.4);

	// Defines a point light source.
	raytracer.addLightSource( new ParallelogramLight(Point3D(0, 0, 5), Vector3D(2, 0, 0), Vector3D(0, 2, 0),
				Colour(0.9, 0.9, 0.9) ) );
	//raytracer.addLightSource( new ParallelogramLight(Point3D(0, 0, 5), Vector3D(2, 0, 0), Vector3D(0, 2, 0),
	//			Colour(0.9, 0.9, 0.9) ) );

	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &ruby );
	SceneDagNode* sphere2 = raytracer.addObject( new UnitSphere(), &jade );
	SceneDagNode* ground = raytracer.addObject( new UnitSquare(), &metal );
	SceneDagNode* left_plane = raytracer.addObject( new UnitSquare(), &gold );
	SceneDagNode* up_plane = raytracer.addObject( new UnitSquare(), &ruby );
	//SceneDagNode* plane4 = raytracer.addObject( new UnitSquare(), &glass );
	
	// Apply some transformations to the unit square.
	double factor1[3] = { 1, 1, 1 };
	double factor2[3] = { 6.0, 6.0, 6.0 };
	double factor3[3] = { 2, 2, 1};

	raytracer.translate(sphere, Vector3D(0, 0, -8));	
	raytracer.rotate(sphere, 'z', 45); 
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);
	
	raytracer.translate(sphere2, Vector3D(0, 0, -4));	
	raytracer.scale(sphere2, Point3D(0, 0, 0), factor3);

	raytracer.translate(ground, Vector3D(0, 0, -12));	
	raytracer.scale(ground, Point3D(0, 0, 0), factor2);
	
	raytracer.translate(left_plane, Vector3D(-5, 0, -8));	
	raytracer.rotate(left_plane, 'y', 90); 
	raytracer.scale(left_plane, Point3D(0, 0, 0), factor2);
	
	raytracer.translate(up_plane, Vector3D(0, 4, -8));	
	raytracer.rotate(up_plane, 'y', 90);
	raytracer.rotate(up_plane, 'x', 90); 
	raytracer.scale(up_plane, Point3D(0, 0, 0), factor2);

	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	raytracer.render(width, height, eye, view, up, fov, "view1.bmp");
	
	// Render it from a different point of view.
	Point3D eye2(4, 2, 1);
	Vector3D view2(-4, -2, -6);
	raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");
	
	return 0;
}

