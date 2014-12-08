//------------------------------------------------------------------------------
//  Copyright 2007-2009 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _BF_MODEL_H_
#define _BF_MODEL_H_

#include <Point.h>
#include <Vector.h>
#include <Matrix.h>
#include <Quaternion.h>
using namespace mathtool;

#include <string>
#include <cassert>
#include <set>
using namespace std;

#include "objReader.h"

//shorten some notations
typedef unsigned int uint;

struct object3D
{
	object3D()
	{
		for (int i = 0; i < 3; i++)
		{
			mat_color[i] = 1;
			mat_specular[i] = 1;
			mat_emission[i] = 0;
		}

		mat_shininess = 128;
	}

	virtual ~object3D(){}

	//material properties
	Vector3d mat_color;
	Vector3d mat_specular;
	Vector3d mat_emission;
	float    mat_shininess;
	
	//This might be helpful!!
	Vector3d position;
};

//a sphere
struct sphere : public object3D
{
	sphere(){ radius = 1; refractive_index = 1; }

	Point3d center;
	double radius;

	//
	Vector3d transparency;  //how transparent this object is
	double refractive_index; //we assume that the refractive index of air is 1
	                         //diamond is 2.419
	                         //amber is about 1.55
	                         //water is 1.333 and ice is about 1.31
	Vector3d reflectiveness; //how mirrow like this object is
};


//a triangle of the model
struct triangle
{
    triangle(){}

    uint v[3]; //vertex id
    uint e[3]; //edge id
	Vector2d texcoord[3]; //texture coordinates

    Vector3d n; //normal

    //backups
    Vector3d bk_n;
};

//a vertex of the model
struct vertex
{
    vertex(){ concave=false; }
    Point3d p;  //position
	Vector3d n; //per vertex normal

    list<uint> m_f;
    list<uint> m_e; //a list of edges

    //backups
    Point3d bk_p; 

    //if concave, set to true
    bool concave;
};

//an edge of the model
struct edge
{
    edge(){ type='x'; vid[0]=vid[1]=UINT_MAX; }
    uint vid[2];
    vector<uint> fid;

    Vector3d v;       //parallel vector
    Vector3d in_n[2]; //inface normals

    //backups
    Vector3d bk_v;       //parallel vector
    Vector3d bk_in_n[2]; //inface normals

    //type, c-convex, r-reflex, p-plane, b-border
    char type;
};

struct model : public object3D
{
    //initialization
    model()
	{ 
        v_size=e_size=t_size=0;
        vertices=NULL;
        edges=NULL;
        tris=NULL; 

        for(int i=0;i<3;i++) for(int j=0;j<3;j++) current_rot[i][j]=0;
        current_rot[0][0]=current_rot[1][1]=current_rot[2][2]=1;
    }

    ~model(){}

    void destroy()
	{
        delete [] tris;     tris=NULL;
        delete [] edges;    edges=NULL;
        delete [] vertices; vertices=NULL;
        v_size=e_size=t_size=0;
    }

    //build this model
    bool build(const string & name);

	//transform this model
	void transform(const Vector3d& T, const Matrix3x3& M, double S);

    //rotate points
    void rotate(const Matrix2x2& m);
    void rotate(const Matrix3x3& M);

    //scale the model
    void scale(double s);

    //negate point/facets ...
    void negate();

    //reverse facets ...
    void reverse();
    
    //data
    vertex   * vertices;  //vertices
    triangle * tris;      //triangles
    edge     * edges;     //edges
    uint v_size, e_size, t_size;

    //current orientation
    double   current_rot[3][3];
};



#endif //_BF_MODEL_H_
