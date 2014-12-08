#include "raytracer.h"
#include "GL/gliLight.h"
#include <algorithm>
#include <cfloat>

//NEW in PA7
#define MAX_RAY_DEPTH 10		//5
#define SHADOW_SAMPLE_SIZE 20	//10

//some helper functions
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }

inline int toInt(double x){ return int( clamp(x)*255 + .5); }


//constructor
RayTracer::RayTracer(list<object3D*>& models) : m_models(models)
{
	//We need to know the camera position
	camera = gli::getCameraPos();

	//Things for gluUnProject
	modelM = new GLdouble[16];	//store the modelview matrix in this
	proj = new GLdouble[16];		//store the projection matrix in this
	view = new GLint[4];				//store the viewport matrix in this
	nearx = new GLdouble;			//store the values of the world coordinates in these
	neary = new GLdouble;
	nearz = new GLdouble;
	farx = new GLdouble;
	fary = new GLdouble;
	farz = new GLdouble;

	glGetDoublev(GL_MODELVIEW_MATRIX, modelM);	//get the modelview, projection, and viewport
	glGetDoublev(GL_PROJECTION_MATRIX, proj);
	glGetIntegerv(GL_VIEWPORT, view);

	srand(time(NULL));

	//find out where lights are
	m_lights.clear();
	for (list<object3D*>::iterator i = models.begin(); i != models.end(); i++)
	{
		object3D* obj = *i;
		if (obj->mat_emission.normsqr() != 0){  //This is a light
			m_lights.push_back(*i);

		}
	}

	if (m_lights.empty())
	{
		cerr << "! Error: There no light defined. Make sure that some object has has non-zero emission." << endl;
		exit(0);
	}
}

inline void show_progress_bar(double ratio)
{
	// Show the load bar.
	static const int w = 50;
	int   c = (int)(ratio * w);

	cout << setw(3) << (int)(ratio * 100) << "% [";
	for (int x = 0; x<c; x++) cout << "=";
	for (int x = c; x<w; x++) cout << " ";
	cout << "]\r" << flush;
}

//render an image with "nray_per_pixel" rays created per pixel
void RayTracer::render(unsigned int img_w, unsigned int img_h, unsigned int nray_per_pixel)
{
	//create a black image of size img_w X img_h
	m_img = vector< vector< Vector3d > >(img_h, vector< Vector3d >(img_w, Vector3d(0,0,0) ));
	 
	//initialize the bar
	show_progress_bar(0);

	//generate rays
	for (unsigned int y = 0; y < img_h; y++)
	{

		for (unsigned int x = 0; x < img_w; x++)
		{
			Vector3d color;
			for (unsigned int n = 0; n < nray_per_pixel; n++)
			{
				Ray r=create_a_random_ray(x,y);

				//Intersection called in raycolor now
				Vector3d rc = raycolor(r, 0);
				color = rc + color;
			}

			m_img[y][x] = color / nray_per_pixel;
		}//end of x

		show_progress_bar(y*1.0 / img_h);

	}//end of y
	show_progress_bar(1.0f);
}

// render an image with "nray_per_pixel" rays created per pixel
// create this ray randomly in the given pixel (x,y)
Ray RayTracer::create_a_random_ray(unsigned int x, unsigned int y)
{
	Ray r;

	//TODO: implement this
	//hint: see slides on generating rays for perspective views

	//offset X, Y by tiny amount
	double xf = drand48() + (double)x;
	double yf = drand48() + (double)(view[3] - y);

	gluUnProject(xf, yf, 0, modelM, proj, view, nearx, neary, nearz);  //set pointPos
	gluUnProject(xf, yf, 1, modelM, proj, view, farx, fary, farz);

	//Just to make sure we're using correct data structures...
	Vector3d nearish = Vector3d(*nearx, *neary, *nearz);
	Vector3d farish = Vector3d(*farx, *fary, *farz);
	Vector3d cameraVec = Vector3d(camera[0], camera[1], camera[2]);

	r.o = nearish + cameraVec;

	//ray direction is normalized far - near
	Vector3d rayDir = farish - nearish;
	r.v = rayDir.normalize();

	//show rays to debug
	//all_rays.push_back(r);

	return r;
}

//returns a model and triangle that intersect ray r
//return pair<NULL,NULL> if no intersection is found
pair<object3D *, triangle *> RayTracer::intersect(Ray r, Point3d& pos)
{
	double min_dist = FLT_MAX;
	triangle * closest_T = NULL;
	object3D * closest_M = NULL;

	for (list<object3D*>::iterator i = m_models.begin(); i != m_models.end(); i++)
	{
		object3D* obj = *i;

		if (dynamic_cast<model*>(obj)) //this object is a mesh
		{
			model* mesh = dynamic_cast<model*>(obj);
			triangle * t = closest_intersect(*mesh, r);
			if (t != NULL)
			{
				Point3d x;
				intersect(*mesh, t, r, x);
				double dist = (x - r.o).normsqr();
				if (dist < min_dist)
				{
					min_dist = dist;
					closest_T = t;
					closest_M = mesh;
					pos = x;
				}
			}
		}//--------------------------------------------------------------------------
		else if (dynamic_cast<sphere*>(obj)) //this object is a sphere
		{
			sphere* ball = dynamic_cast<sphere*>(obj);
			Point3d x;
			if (intersect(*ball, r, x))
			{
				double dist = (x - r.o).normsqr();
				if (dist < min_dist)
				{
					min_dist = dist;
					closest_T = NULL;
					closest_M = ball;
					pos = x;
				}
			}
		}//--------------------------------------------------------------------------
	}

	return make_pair(closest_M, closest_T);
}

//
//returns true if the sphere intersects ray r
//return false if no intersection is found
//x is the location of intersection if true is returned.
//
bool RayTracer::intersect(sphere& s, Ray r, Point3d& x)
{
	//
	// PA7 TODO: implement this
	// intersecting a sphere
	Point3d com = s.center;
	double radius = s.radius;

	//crazy calcs
	float a = (r.v[0] * r.v[0]) + (r.v[1] * r.v[1]) + (r.v[2] * r.v[2]);
	float b = (2.0 * r.v[0] * (r.o[0] - com[0])) + (2.0 * r.v[1] * (r.o[1] - com[1])) + (2.0 * r.v[2] * (r.o[2] - com[2]));
	float c = (com[0] * com[0]) + (com[1] * com[1]) + (com[2] * com[2]) + (r.o[0] * r.o[0]) + (r.o[1] * r.o[1]) + (r.o[2] * r.o[2]) +
		(-2.0 * ((com[0] * r.o[0]) + (com[1] * r.o[1]) + (com[2] * r.o[2]))) - (radius*radius);
	
	//discriminate
	float discriminant = pow(b, 2) - (4 * a * c);

	//Ray intersects sphere
	if (discriminant > 0){
		float t = (-b - sqrt(discriminant)) / (2.0 * a);
		if (t < 0.0) return false;		//both must be above 0

		x[0] = r.o[0] + (t * r.v[0]);
		x[1] = r.o[1] + (t * r.v[1]);
		x[2] = r.o[2] + (t * r.v[2]);
		return true;
	}
	return false;  //ray doesnt intersect sphere

	/*
	a = dx*dx + dy*dy + dz*dz;

            b = 2*dx*(x0-cx) +  2*dy*(y0-cy) +  2*dz*(z0-cz);

            c = cx*cx + cy*cy + cz*cz + x0*x0 + y0*y0 + z0*z0 +

                                    -2*(cx*x0 + cy*y0 + cz*z0) - R*R;
	*/
}


//returns a triangle in model m that intersect ray r
//return NULL if no intersection is found
triangle * RayTracer::intersect(model& m, Ray r)
{
	for (unsigned int i = 0; i < m.t_size; i++)
	{
		Point3d x;
		if ( intersect(m, &m.tris[i], r, x) )
			return &m.tris[i];
	}

	return NULL;
}

//returns a triangle in model m that make closest intersection with ray r
//return NULL if no intersection is found
triangle * RayTracer::closest_intersect(model& m, Ray r)
{
	double min_dist = FLT_MAX;
	triangle * closest_T=NULL;
	for (unsigned int i = 0; i < m.t_size; i++)
	{
		Point3d x;
		if (intersect(m, &m.tris[i], r, x))
		{
			double dist = (x - r.o).normsqr();
			if (dist < min_dist)
			{
				closest_T = &m.tris[i];
				min_dist = dist;
			}
		}//end if
	}//end for i

	return closest_T;
}

// determine if there is an intersection between triangle t and ray r
// return true if there is an intersection and store the intersection in x
// return false otherwise and x is undefined in this case
bool RayTracer::intersect(model& m, triangle * t, Ray r, Point3d& x)
{
	//PA6 TODO: implement this - Same as project 6

	//Get vertices of Triangle
	Vector3d v0 = Vector3d(m.vertices[t->v[0]].p[0], m.vertices[t->v[0]].p[1], m.vertices[t->v[0]].p[2]);
	Vector3d v1 = Vector3d(m.vertices[t->v[1]].p[0], m.vertices[t->v[1]].p[1], m.vertices[t->v[1]].p[2]);
	Vector3d v2 = Vector3d(m.vertices[t->v[2]].p[0], m.vertices[t->v[2]].p[1], m.vertices[t->v[2]].p[2]);

	//Ray stuffs
	Vector3d rayPos = Vector3d(r.o[0], r.o[1], r.o[2]);
	Vector3d rayDir = Vector3d(r.v[0], r.v[1], r.v[2]);

	Vector3d point; // = new Vector3d();
	bool hits = rayIntersectsTriangle(rayPos, rayDir, v0, v1, v2, point);

	x[0] = point[0];
	x[1] = point[1];
	x[2] = point[2];
	return hits;
}

bool RayTracer::rayIntersectsTriangle(Vector3d pos, Vector3d dir,
	Vector3d v0, Vector3d v1, Vector3d v2, Vector3d& point) {

	Vector3d edge1, edge2, h, s, q;
	float a, f, u, v, t;
	edge1 = v1 - v0;
	edge2 = v2 - v0;

	h = dir % edge2;
	a = edge1 * h;

	if (a > -0.00001 && a < 0.00001) return false;  //checking for perpendicularism ?

	f = 1 / a;
	s = pos - v0;
	u = f * (s * h);

	if (u < 0.0 || u > 1.0) return false;

	q = s % edge1;
	v = f * (dir * q);

	if (v < 0.0 || u + v > 1.0) return(false);

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	t = f * (edge2 * q);

	if (t > 0.00001){ // ray intersection
		point = pos + dir * t;
		return true;
	}
	else // this means that there is a line intersection
		// but not a ray intersection
		return false;

}


//
// determine the color of ray r
// 
Vector3d RayTracer::raycolor(const Ray& r, int depth)
{
	if (depth >= MAX_RAY_DEPTH) return Vector3d(0,0,0);

	Vector3d color;
	Point3d pos;	//The position of the intersection
	pair<object3D *, triangle *> X = intersect(r, pos);

	//determine the color of this ray 
	if (X.first != NULL && X.second != NULL) //this is mesh intersection
	{
		Vector3d rc = raycolor(*((model*)X.first), X.second, r,depth, pos);
		color = rc + color;
	}
	else if (X.first != NULL) //this is spehere intersection
	{
		Vector3d rc = raycolor(*((sphere*)X.first), r, depth, pos);
		color = rc + color;
	}

	return color;
}

//
// determine the color of ray r, by analizing the intersection between t and r 
// 
Vector3d RayTracer::raycolor(sphere& s, const Ray& r, int depth, Point3d pos)
{
	Vector3d color;

	//PA7 TODO: determine the direct color
	//
	//hint: (1) what is the normal direction at the intersection between s and r
	Vector3d normal = (Vector3d(pos.get()) - Vector3d(s.center.get())).normalize();  //vector from com to intersect point
	Vector3d reflection_color;
	Vector3d refraction_color;
	Vector3d lightPos;
	Vector3d lightToPoint;
	//      (2) what is the color of this sphere?  //s.mat_color
	color = Vector3d();

	//      (3) where is the light?
	for (list<object3D*>::iterator i = m_lights.begin(); i != m_lights.end(); i++)
	{
		object3D* light = *i;
		lightPos = Vector3d(light->position.get());  //lightSource
		lightToPoint = (lightPos - Vector3d(pos.get())).normalize();	

		color = color + ((s.mat_color ^ light->mat_emission) * max((normal * lightToPoint), 0.0));
	}
	//      (4) is the intersection in shadow? (call inshadow(const Point3d& p))
	double gradient = inshadow(pos);
	if (gradient > 0.0){		//intersection is in shadow
		color = color * (1.0 - gradient);
	}


	if (s.reflectiveness.normsqr() > 0)
	{
		//PA7 TODO: determine the reflection color
		//
		//      (1) generate reflection ray 
		//      (2) determine the color of the ray (remember of increase depth by 1)
		Ray refRay;


		refRay.v = (r.v - (normal * ((r.v * normal) * 2.0)));	//formula for reflection over normal
		
		refRay.o = Point3d((Vector3d(pos.get()) + (refRay.v * .01)).get());		//push origin a little off face

		reflection_color = raycolor(refRay, depth + 1);		//recurse
	}

	if (s.transparency.normsqr() > 0)
	{
		//      (1) generate refraction ray 
		//      (2) determine the color of the ray (remember of increase depth by 1)

		double ref_ind;
		Ray tranRay;

		if ((r.o - s.center).norm() - s.radius <= FLT_EPSILON){  //leaving the object
			ref_ind = s.refractive_index / 1.0006;
		}
		else{													//entering the object
			ref_ind = 1.0006 / s.refractive_index;
		}
		
		double work = normal * r.v;
		double plz = 1.0f - ((ref_ind * ref_ind) * (1.0f - (work * work)));  //Simplify order of operations

		tranRay.v = ((normal * (ref_ind * work - sqrt(plz))) - (r.v * ref_ind)).normalize();
		//((n * (N * I) - (sqrt(1 - (pow(n,s) * (1 - pow(N * I, 2)))))) * N) - (n * I)

		tranRay.o = Point3d((Vector3d(pos.get()) + (tranRay.v * .01)).get());

		refraction_color = raycolor(tranRay, depth + 1);

	}

	//finally gather all the colors
	for (int i = 0; i < 3; i++)
	{
		color[i] = clamp(color[i] + s.reflectiveness[i] * reflection_color[i] + s.transparency[i] * refraction_color[i]);
	}

	return color;
}


//
// determine the color of ray r, by analizing the intersection between t and r 
// 
Vector3d RayTracer::raycolor(model& m, triangle * t, const Ray& r, int depth, Point3d pos)
{
	if (m.mat_emission.normsqr() != 0.0){

		return m.mat_emission;		//this model is a light
	}
	Vector3d color;
	//PA6 TODO: implement this - Same as Project 6?

	vertex v0 = m.vertices[t->v[0]];
	vertex v1 = m.vertices[t->v[1]];
	vertex v2 = m.vertices[t->v[2]];

	//TODO: implement this
	//Interpolate the normal
	Vector3d weights = getBarycentricCoordinatesAt(Vector3d(pos.get()), Vector3d(v0.p.get()), Vector3d(v1.p.get()), Vector3d(v2.p.get()), t->n);
	Vector3d interNorm = (weights[0] * v0.n) + (weights[1] * v1.n) + (weights[2] * v2.n);

	//This is the color
	//color = m.mat_color;

	//Get the light and eye vector?
	Vector3d lightPos = Vector3d(light0_position[0], light0_position[1], light0_position[2]);
	Vector3d eye = (Vector3d(camera[0], camera[1], camera[2]) - Vector3d(pos.get())).normalize();
	Vector3d fragToLight = (lightPos - Vector3d(pos.get())).normalize();

	//Get the diffusion and specular value
	float difTemp = fragToLight * interNorm;
	Vector3d diffuse = m.mat_color * max(difTemp, 0.0f);
	Vector3d lye = (fragToLight + eye).normalize();
	Vector3d specular = m.mat_specular * pow(max(lye * interNorm, 0.0), m.mat_shininess);

	color = diffuse + specular;

	//is the frag facing away from the light?
	if (difTemp <= 0){
		//diffuse = 0;
		return color;
	}
	//diffuse not less than zero
	else{
		double gradient = inshadow(pos);
		if (gradient > 0.0){
			//std::cout << color << " " << gradient << std::endl;
			color = color * (1.0 - gradient);
			//std::cout << color << std::endl;
			return color;
		}
	}

	return color;
}

//check if a point p is in shadow
//return 0 is p is not in the shadow
//a value between 0 and 1 to indicate the "probability"
//that p is in full darkness
double RayTracer::inshadow(const Point3d& p)
{
	double gradient = 0.0;
	Vector3d pos = Vector3d(p.get());
	Vector3d lightPos;
	object3D *light;
	Ray lightR;

	for (int i = 0; i < SHADOW_SAMPLE_SIZE; i++){

		//random location in the light
		list<object3D*>::iterator lightI = m_lights.begin();  //TODO make this for multiple lights
		light = *lightI;
		lightPos = Vector3d(light->position.get());

		lightPos = randomLight(lightPos);
		
		//Shoot ray from light to intersect point
		lightR.v = (lightPos - pos).normalize();
		lightR.o = Point3d((pos + (lightR.v * .0001)).get());  //Push off face

		//Test for intersection
		Point3d lightRef;
		pair<object3D *, triangle *> X = intersect(lightR, lightRef);

		//hit something
		if (X.first != NULL){
			if ((lightR.o - lightRef).norm() < (lightR.o - lightPos).norm()){
				gradient += 1.0;
			}
		}
	}
	
	return (gradient / (double)SHADOW_SAMPLE_SIZE) / 2;
}

//Helper for barycentric normals
Vector3d RayTracer::getBarycentricCoordinatesAt(Vector3d point, Vector3d v0, Vector3d v1, Vector3d v2, Vector3d normal)
{
	Vector3d bary;

	// The area of a triangle is 
	float areaABC = normal * ((v1 - v0) % (v2 - v0));
	float areaPBC = normal * ((v1 - point) % (v2 - point));
	float areaPCA = normal * ((v2 - point) % (v0 - point));

	bary[0] = areaPBC / areaABC; // alpha
	bary[1] = areaPCA / areaABC; // beta
	bary[2] = 1.0f - bary[0] - bary[1]; // gamma

	return bary;
}

//save rendered image to file
bool RayTracer::save2file(const string& filename)
{
	FILE *f = fopen(filename.c_str(), "w");         // Write image to PPM file.

	int h = m_img.size();
	if (h == 0) return true; //nothing to save...

	int w = m_img.front().size();
	
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);

	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
			fprintf(f, "%d %d %d ", toInt(m_img[i][j][0]), toInt(m_img[i][j][1]), toInt(m_img[i][j][2]));

		show_progress_bar(i*1.0 / h);
	}
	show_progress_bar(1.0f);
	fclose(f);
}

Vector3d RayTracer::randomLight(Vector3d light){
	Vector3d newLight;
	newLight[1] = light[1] - 5;
	float somVal1 = -1.0 + (3.0 * rand() / (RAND_MAX + 1.0));
	float somVal2 = -1.0 + (3.0 * rand() / (RAND_MAX + 1.0));
	newLight[0] = light[0] + (somVal1 * 5);
	newLight[2] = light[2] + (somVal2 * 5);
	//haxsar
	return newLight;
}