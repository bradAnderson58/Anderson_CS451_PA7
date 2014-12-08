//------------------------------------------------------------------------------
//  Copyright 2007-2014 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------


#pragma once

#include "objReader.h"
#include "model.h"

#include <list>
#include <float.h>
#include <algorithm>
using namespace std;

//-----------------------------------------------------------------------------
// INPUTS
string rt_filename; //an rt (raytracing) file defines models and their transformation and material properties

unsigned int image_w = 400, image_h = 400; //image size
unsigned int n_ray = 1;

//-----------------------------------------------------------------------------
// Intermediate data
list<object3D *> models;

float R=0;       //radius
Point3d COM;     //center of mass
double BOX[6];   //bounding box

//-----------------------------------------------------------------------------
//read models from file
bool readfromfile();
void computeCOM_R();
void buildBoxModel(double * box);
void printUsage(char * name);

//-----------------------------------------------------------------------------
bool parseArg(int argc, char ** argv)
{
    for(int i=1;i<argc;i++){
        if(argv[i][0]=='-')
        {
			string tmp(argv[i]);
            return false; //unknown
        }
        else{
			rt_filename=argv[i];
        }
    }

	if (rt_filename.empty())
	{
		cerr << "! Error: No RT (Raytracing) file is given." << endl;
		printUsage(argv[0]);
		return false;
	}

	if (rt_filename.find(".rt") == string::npos)
	{
		cerr << "! Error: the given file ("<< rt_filename<<") is not an RT(Raytracing) file." << endl;
		printUsage(argv[0]);
		return false;
	}

    return true;
}

void printUsage(char * name)
{
    //int offset=20;
	cerr << "Usage: " << name << " *.rt \n";
    cerr<<"\n-- Report bugs to: Jyh-Ming Lien jmlien@gmu.edu"<<endl;
}

//-----------------------------------------------------------------------------

void computeCOM_R()
{
	//compute a bbox
	double box[6] = { FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX };
	//-------------------------------------------------------------------------
	for (list<object3D*>::iterator i = models.begin(); i != models.end(); i++)
	{
		object3D* obj = *i;
		
		//this obj is a mesh
		if (dynamic_cast<model*>(obj))
		{
			model* mesh = dynamic_cast<model*>(obj);
			for (unsigned int j = 0; j < mesh->v_size; j++){ //this is a sphere
				Point3d& p = mesh->vertices[j].p;
				if (p[0]<box[0]) box[0] = p[0];
				if (p[0]>box[1]) box[1] = p[0];
				if (p[1]<box[2]) box[2] = p[1];
				if (p[1]>box[3]) box[3] = p[1];
				if (p[2]<box[4]) box[4] = p[2];
				if (p[2]>box[5]) box[5] = p[2];
			}//j
		}
		else if (dynamic_cast<sphere*>(obj)) //this is a ball
		{
			sphere* ball = dynamic_cast<sphere*>(obj);
			double sbox[6] = { ball->center[0] - ball->radius, ball->center[0] + ball->radius, 
				               ball->center[1] - ball->radius, ball->center[1] + ball->radius, 
				               ball->center[2] - ball->radius, ball->center[2] + ball->radius };

			if (sbox[0] < box[0]) box[0] = sbox[0];
			if (sbox[1] > box[1]) box[1] = sbox[1];
			if (sbox[2] < box[2]) box[2] = sbox[2];
			if (sbox[3] > box[3]) box[3] = sbox[3];
			if (sbox[4] < box[4]) box[4] = sbox[4];
			if (sbox[5] > box[5]) box[5] = sbox[5];
		}
	}//i

	//-------------------------------------------------------------------------
	// compute center of mass and R...
	COM.set((box[1] + box[0]) / 2, (box[3] + box[2]) / 2, (box[5] + box[4]) / 2);

	//scale the twice box first
	box[0] = max((box[0] - COM[0]) * 10 + COM[0], box[0] - 30);
	box[1] = min((box[1] - COM[0]) * 10 + COM[0], box[1] + 30);
	box[3] = box[3] + 0.01;
	box[4] = max((box[4] - COM[2]) * 10 + COM[2], box[4] - 30);
	box[5] = min((box[5] - COM[2]) * 10 + COM[2], box[5] + 30);

	COM.set((box[1] + box[0]) / 2, (box[3] + box[2]) / 2, (box[5] + box[4]) / 2);

	//-------------------------------------------------------------------------
	R = (float)(COM - Point3d(box[0], box[2], box[4])).norm();

	//build box walls
	buildBoxModel(box);

	//remember the box
	for (int i = 0; i < 6; i++) BOX[i] = box[i];
}

inline void buildBoxWall(const Point3d& a, const Point3d& b, const Point3d& c, const Point3d& d, const Vector3d& n, const Vector3d& color)
{
	model * m=new model();
	assert(m);

	//build left wall
	m->v_size = 4;
	m->vertices = new vertex[4];
	m->vertices[0].p = a;
	m->vertices[0].n = n;
	m->vertices[1].p = b;
	m->vertices[1].n = n;
	m->vertices[2].p = c;
	m->vertices[2].n = n;
	m->vertices[3].p = d;
	m->vertices[3].n = n;

	m->e_size = 5;
	m->edges = new edge[5];
	m->edges[0].vid[0] = 0;
	m->edges[0].vid[1] = 1;
	m->edges[1].vid[0] = 1;
	m->edges[1].vid[1] = 2;
	m->edges[2].vid[0] = 2;
	m->edges[2].vid[1] = 3;
	m->edges[3].vid[0] = 3;
	m->edges[3].vid[1] = 0;
	m->edges[4].vid[0] = 2;
	m->edges[4].vid[1] = 0;

	m->t_size = 2;
	m->tris = new triangle[2];
	m->tris[0].v[0] = 0;
	m->tris[0].v[1] = 1;
	m->tris[0].v[2] = 2;
	m->tris[0].n = n;
	m->tris[1].v[0] = 2;
	m->tris[1].v[1] = 3;
	m->tris[1].v[2] = 0;
	m->tris[1].n = n;

	m->mat_color = color;
	m->mat_specular = color;

	models.push_back(m);
}

//create a model with the front open...
void buildBoxModel(double * box)
{
	Point3d llb(box[0], box[2], box[4]), lrb(box[1], box[2], box[4]), urb(box[1], box[3], box[4]), ulb(box[0], box[3], box[4]),
	    	llf(box[0], box[2], box[5]), lrf(box[1], box[2], box[5]), urf(box[1], box[3], box[5]), ulf(box[0], box[3], box[5]);

	//build left wall (yellow)
	buildBoxWall(llf, llb, ulb, ulf, Vector3d(1,0,0), Vector3d(1,1,0));
	
	//build right wall (green)
	buildBoxWall(lrf, urf, urb, lrb, Vector3d(-1, 0, 0), Vector3d(0, 1, 0));

	//build back wall (white)
	buildBoxWall(llb, lrb, urb, ulb, Vector3d(0, 0, 1), Vector3d(1, 1, 1));

	//build floor (white)
	buildBoxWall(llf, lrf, lrb, llb, Vector3d(0, 1, 0), Vector3d(1, 1, 1));

	//build ceiling (blue)
	buildBoxWall(ulf, ulb, urb, urf, Vector3d(0, -1, 0), Vector3d(0, 0, 1));
}


//-----------------------------------------------------------------------------
//
//
//
//  parsing objects
//
//
//-----------------------------------------------------------------------------

//convert a string to a list of tokens
inline list<string> tokenize(char * tmp, const char * ignore)
{
	list<string> tokens;
	char * tok = strtok(tmp, ignore);
	while (tok != NULL)
	{
		tokens.push_back(tok);
		tok = strtok(NULL, ignore);
	}
	return tokens;
}


inline bool findImageSize(istream& in)
{
	const int size = 1024;
	char * tmp = new char[size];
	in.getline(tmp, size);//read lines from a file into strings
	list<string> tok = tokenize(tmp, " =\t[]()<>,");
	if (tok.size() == 0) return false;			//empty line
	if (tok.front()[0] == '#') return false; 	//commented out

	for (list<string>::iterator i = tok.begin(); i != tok.end(); i++){
		string& t = *i;
		if (t[0] == '#') break;
		if (t == "width") image_w = (unsigned int)atoi((++i)->c_str());
		else if (t == "height") image_h = (unsigned int)atoi((++i)->c_str());
		else if (t == "n") n_ray = (unsigned int)atoi((++i)->c_str());
	}
	delete[] tmp;
	return true;
}

inline bool createModel(istream& in)
{
	const int size = 1024;
	char * tmp = new char[size];
	in.getline(tmp, size);//read lines from a file into strings
	list<string> tok = tokenize(tmp, " =\t[]()<>,");
	if (tok.size() == 0) return false;			//empty line
	if (tok.front()[0] == '#') return false; 	//commented out

	Vector3d m_local_rot(0, 0, 0);
	Vector3d m_pos(0, 0, 0);
	double s=1; //scale
	string filename;
	Vector3d mat_color = Vector3d(1, 1, 1);
	Vector3d mat_specular = Vector3d(1, 1, 1);
	Vector3d mat_emission = Vector3d(0, 0, 0);
	float    mat_shininess = 128;

	for (list<string>::iterator i = tok.begin(); i != tok.end(); i++){
		string& t = *i;
		if (t[0] == '#') break;
		if (t == "tx") m_pos[0] = (float)atof((++i)->c_str());
		else if (t == "ty") m_pos[1] = (float)atof((++i)->c_str());
		else if (t == "tz") m_pos[2] = (float)atof((++i)->c_str());
		else if (t == "rx") m_local_rot[0] = (float)atof((++i)->c_str());
		else if (t == "ry") m_local_rot[1] = (float)atof((++i)->c_str());
		else if (t == "rz") m_local_rot[2] = (float)atof((++i)->c_str());
		else if (t == "scale") s = (float)atof((++i)->c_str());
		//difuss color
		else if (t == "cr") mat_color[0] = (float)atof((++i)->c_str());
		else if (t == "cg") mat_color[1] = (float)atof((++i)->c_str());
		else if (t == "cb") mat_color[2] = (float)atof((++i)->c_str());
		//specular color
		else if (t == "sr") mat_specular[0] = (float)atof((++i)->c_str());
		else if (t == "sg") mat_specular[1] = (float)atof((++i)->c_str());
		else if (t == "sb") mat_specular[2] = (float)atof((++i)->c_str());
		//emission color
		else if (t == "er") mat_emission[0] = (float)atof((++i)->c_str());
		else if (t == "eg") mat_emission[1] = (float)atof((++i)->c_str());
		else if (t == "eb") mat_emission[2] = (float)atof((++i)->c_str());
		else if (t == "shininess") mat_shininess = (float)atof((++i)->c_str());
		else filename = t;
	}

	delete[] tmp;

	Quaternion rot = Quaternion(m_local_rot.get()); //rz*ry*rx;

	model * m=new model();
	assert(m);
	if (!m->build(filename)) return false;

	m->transform(m_pos, rot.getMatrix(), s);
	m->mat_color = mat_color;
	m->mat_specular = mat_specular;
	m->mat_emission = mat_emission;
	m->mat_shininess = mat_shininess;

	m->position = m_pos;

	models.push_back(m);

	return true;
}


inline bool createSphere(istream& in)
{
	const int size = 1024;
	char * tmp = new char[size];
	in.getline(tmp, size);//read lines from a file into strings
	list<string> tok = tokenize(tmp, " =\t[]()<>,");
	if (tok.size() == 0) return false;			//empty line
	if (tok.front()[0] == '#') return false; 	//commented out

	Vector3d m_pos(0, 0, 0);
	double radius = 1; //radius

	Vector3d mat_color = Vector3d(1, 1, 1);
	Vector3d mat_specular = Vector3d(1, 1, 1);
	Vector3d mat_emission = Vector3d(0, 0, 0);
	double    mat_shininess = 128;

	Vector3d transparency = Vector3d(0, 0, 0);  //how transparent this object is
	Vector3d reflectiveness = Vector3d(0, 0, 0); //how mirrow like this object is
	double refractive_index = 1;

	for (list<string>::iterator i = tok.begin(); i != tok.end(); i++){
		string& t = *i;
		if (t[0] == '#') break;
		if (t == "tx") m_pos[0] = (float)atof((++i)->c_str());
		else if (t == "ty") m_pos[1] = (float)atof((++i)->c_str());
		else if (t == "tz") m_pos[2] = (float)atof((++i)->c_str());
		else if (t == "r") radius = (float)atof((++i)->c_str());
		//difuss color
		else if (t == "cr") mat_color[0] = (float)atof((++i)->c_str());
		else if (t == "cg") mat_color[1] = (float)atof((++i)->c_str());
		else if (t == "cb") mat_color[2] = (float)atof((++i)->c_str());
		//specular color
		else if (t == "sr") mat_specular[0] = (float)atof((++i)->c_str());
		else if (t == "sg") mat_specular[1] = (float)atof((++i)->c_str());
		else if (t == "sb") mat_specular[2] = (float)atof((++i)->c_str());
		//emission color
		else if (t == "er") mat_emission[0] = (float)atof((++i)->c_str());
		else if (t == "eg") mat_emission[1] = (float)atof((++i)->c_str());
		else if (t == "eb") mat_emission[2] = (float)atof((++i)->c_str());
		else if (t == "shininess") mat_shininess = (float)atof((++i)->c_str());
		//transparency
		else if (t == "Tr") transparency[0] = (float)atof((++i)->c_str());
		else if (t == "Tg") transparency[1] = (float)atof((++i)->c_str());
		else if (t == "Tb") transparency[2] = (float)atof((++i)->c_str());
		//relectiveness (how mirrow-like this object is)
		else if (t == "Rr") reflectiveness[0] = (float)atof((++i)->c_str());
		else if (t == "Rg") reflectiveness[1] = (float)atof((++i)->c_str());
		else if (t == "Rb") reflectiveness[2] = (float)atof((++i)->c_str());
		//refractive index
		else if (t == "RI") refractive_index = (float)atof((++i)->c_str());
	}

	delete[] tmp;

	sphere * m = new sphere();
	assert(m);
	m->center = m_pos;
	m->radius = radius;
	m->mat_color = mat_color;
	m->mat_specular = mat_specular;
	m->mat_emission = mat_emission;
	m->mat_shininess = (float)mat_shininess;
	m->transparency = transparency;
	m->reflectiveness = reflectiveness;
	m->refractive_index = (float)refractive_index;
	
	//position
	m->position = m_pos;
	
	models.push_back(m);

	return true;
}


inline bool initialize(istream& in)
{
	const int size = 1024;
	char * tmp = new char[size];

	while (!in.eof())
    {
		in.getline(tmp, size);
		list<string> tok = tokenize(tmp, " \t\r[]()<>,=");
		if (tok.empty()) continue;

		string label = tok.front();
		tok.pop_front();
		if (label[0] == '#') continue;
		
		//initialize image
		if (label == "image" || label == "Image" || label == "IMAGE")
		{
			findImageSize(in);
		}//end if obst

		//initialize meshes
		else if (label == "model" || label == "Model" || label == "MODEL")
		{
			int m_size = atoi(tok.front().c_str());
			for (int i = 0; i<m_size; i++)
			{
				createModel(in);
			}//end outer for
		}//end if obst


		//initialize sphere
		else if (label == "sphere" || label == "Sphere" || label == "SPHERE")
		{
			int sphere_size = atoi(tok.front().c_str());
			for (int i = 0; i<sphere_size; i++)
			{
				createSphere(in);
			}//end outer for
		}//end if obst

		//Unknown
		else{
			cerr << "! Warning: Unknown value: (" << tmp <<")"<< endl;
			continue;
		}
	}

	delete[] tmp;

	return true;
}

inline bool initialize(const string& filename)
{
	ifstream fin(filename.c_str());

	if (fin.good() == false) return false;

	if (initialize(fin) == false) return false;

	fin.close();

	computeCOM_R();

	return true;
}

//-----------------------------------------------------------------------------
//
//
//
//  Open GL stuff below
//
//
//-----------------------------------------------------------------------------

#include <draw.h>



