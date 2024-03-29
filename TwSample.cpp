//  ---------------------------------------------------------------------------
//
//  @file       TwSimpleGLUT.cpp
//  @brief      A simple example that uses AntTweakBar with OpenGL and GLUT.
//
//              AntTweakBar: http://anttweakbar.sourceforge.net/doc
//              OpenGL:      http://www.opengl.org
//              GLUT:        http://opengl.org/resources/libraries/glut
//  
//  @author     Philippe Decaudin
//  @date       2006/05/20
//
//  Modified by Guoning Chen
//  @date       2018/08/09
//
//  Modified by Binoy Dalal
//  @date       2018/09/18
//  ---------------------------------------------------------------------------


#include <AntTweakBar.h>

#include <stdlib.h>
#include <stdio.h>
#include <vector>


#define _USE_MATH_DEFINES
#include <math.h>

#if defined(_WIN32) || defined(_WIN64)
//  MiniGLUT.h is provided to avoid the need of having GLUT installed to 
//  recompile this example. Do not use it in your own programs, better
//  install and use the actual GLUT library SDK.
#   define USE_MINI_GLUT
#endif

#if defined(USE_MINI_GLUT)
#   include "./src/MiniGLUT.h"
#elif defined(_MACOSX)
#   include <GLUT/glut.h>
#else
#   include <GL/glut.h>
#endif

// Main window ID
int MainWindow;

// This example displays one of the following shapes
//typedef enum { SHAPE_TEAPOT=1, SHAPE_TORUS, SHAPE_CONE, BUNNY } Shape;

#define NUM_SHAPES 3
//Shape g_CurrentShape = SHAPE_TORUS;
int g_CurrentShape = 0;

// Shapes scale
float g_Zoom = 1.0f;
// Shape orientation (stored as a quaternion)
float g_Rotation[] = { 0.0f, 0.0f, 0.0f, 1.0f };
// Auto rotate
int g_AutoRotate = 0;
int g_RotateTime = 0;
float g_RotateStart[] = { 0.0f, 0.0f, 0.0f, 1.0f };
// Shapes material
float g_MatAmbient[] = { 0.5f, 0.0f, 0.0f, 1.0f };
float g_MatDiffuse[] = { 1.0f, 1.0f, 0.0f, 1.0f };
// Light parameter
float g_LightMultiplier = 1.0f;
float g_LightDirection[] = { -0.57735f, -0.57735f, -0.57735f };
double s_min, s_max, a_max, a_min, vx_max, vx_min, vy_max, vy_min;
float max_sv1, max_sv2, max_sv3, min_sv1, min_sv2, min_sv3;
float max_vx1, max_vx2, max_vx3, min_vx1, min_vx2, min_vx3;
float max_vy1, max_vy2, max_vy3, min_vy1, min_vy2, min_vy3;
float max_vz1, max_vz2, max_vz3, min_vz1, min_vz2, min_vz3;
float abs_s_min = 0., abs_s_max = 0.;
// pointers to max and min values in any dataset
float *max_ptr, *min_ptr;
// White threshold
float g_WhiteThreshold = 0.5;

TwBar *bar = NULL; // Pointer to the tweak bar

void(*colorFunction)(float, float[]); // pointer to color function of choice

void(*vectorFunction)(float, float, float, float&, float&, float&, float&); // pointer to vector field of choice

// some sample color definitions:
// this order must list order of the colors

const GLfloat Colors[7][3] =
{
	{ 1., 0., 0. },		// red
	{ 1., 1., 0. },		// yellow
	{ 0., 1., 0. },		// green
	{ 0., 1., 1. },		// cyan
	{ 0., 0., 1. },		// blue
	{ 1., 0., 1. },		// magenta
	{ 1., 1., 1. },		// white
};
#define NUM_COLORS 5
int whichColor = 0;
int whichPlot = 0;

// the stroke characters 'X' 'Y' 'Z' :

static float xx[] = {
	0.f, 1.f, 0.f, 1.f
};

static float xy[] = {
	-.5f, .5f, .5f, -.5f
};

static int xorder[] = {
	1, 2, -3, 4
};


static float yx[] = {
	0.f, 0.f, -.5f, .5f
};

static float yy[] = {
	0.f, .6f, 1.f, 1.f
};

static int yorder[] = {
	1, 2, 3, -2, 4
};


static float zx[] = {
	1.f, 0.f, 1.f, 0.f, .25f, .75f
};

static float zy[] = {
	.5f, .5f, -.5f, -.5f, 0.f, 0.f
};

static int zorder[] = {
	1, 2, 3, 4, -5, 6
};

// size of the box:

const float BOXSIZE = { 2.f };


// fraction of the length to use as height of the characters:

const float LENFRAC = 0.10f;


// fraction of length to use as start location of the characters:

const float BASEFRAC = 1.10f;

// line width for the axes:

const GLfloat AXES_WIDTH = { 3. };


GLuint	BoxList = 100;		// object display list
GLuint	AxesList = 101;		// list to hold the axes

int g_Axes = 1;   // Toggle Axes
int g_Box = 0;    // Toggle Box

int g_arrows = 0; // Toggle Arrow heads
int g_streamlines = 0; // Toggle streamlines
int g_streamribbon = 0; // Toggle streamribbon
float g_probeX = 0.; // Probe X coordinate
float g_probeY = 0.; // Probe Y coordinate
float g_probeZ = 0.; // Probe Z coordinate
float g_step = 0.1; //  step size when integrating
float g_linedist = 0.1; // line dist for computing stream bunch
float g_rdist = 0.1; // width of streamribbon
int whichIntegrator = 0; // euler or Rk2

#include "Skeleton.h"
Polyhedron *poly = NULL;

typedef struct vecNode {
	float x, y, z;
	float vx, vy, vz;
	float magnitude;
};

const int NX3d = 8, NY3d = 8, NZ3d = 8;

vecNode vec_field1[NX3d][NY3d][NZ3d];
vecNode vec_field2[NX3d][NY3d][NZ3d];
vecNode vec_field3[NX3d][NY3d][NZ3d];
int currentField;

//streamline length
int g_streamLength = 200;

typedef struct streampoint {
	float nextX, nextY, nextZ;
	float magnitude;
};

std::vector<streampoint> streamline;
std::vector<streampoint> streamr1;
std::vector<streampoint> streamr2;
std::vector<std::vector<streampoint>> streamlines;

const int MINUS = { 0 };
const int PLUS = { 1 };

#define X 0
#define Y 1
#define Z 2
#define WINGS 0.10

/* x, y, z, axes: */
static float axx[3] = { 1., 0., 0. };
static float ayy[3] = { 0., 1., 0. };
static float azz[3] = { 0., 0., 1. };

const float A = sqrt(3);
const float B = sqrt(2);

void Display(void);
void Reshape(int width, int height);
void Arrow(float[3], float[3]);
void cross(float[3], float[3], float[3]);
float dot(float[3], float[3]);
float unit(float[3], float[3]);

void Arrow(float tail[3], float head[3])
{
	float u[3], v[3], w[3]; /* arrow coordinate system */
	float d; /* wing distance */
	float x, y, z; /* point to plot */
	float mag; /* magnitude of major direction */
	float f; /* fabs of magnitude */
	int axis; /* which axis is the major */
	/* set w direction in u-v-w coordinate system: */
	w[0] = head[0] - tail[0];
	w[1] = head[1] - tail[1];
	w[2] = head[2] - tail[2];
	/* determine major direction: */
	axis = X;
	mag = fabs(w[0]);
	if ((f = fabs(w[1])) > mag)
	{
		axis = Y;
		mag = f;
	}
	if ((f = fabs(w[2])) > mag)
	{
		axis = Z;
		mag = f;
	}
	/* set size of wings and turn w into a unit vector: */
	d = WINGS * unit(w, w);
	/* draw the shaft of the arrow: */
	glBegin(GL_LINE_STRIP);
	glVertex3fv(tail);
	glVertex3fv(head);
	glEnd();
	/* draw two sets of wings in the non-major directions: */
	if (axis != X)
	{
		cross(w, axx, v);
		(void)unit(v, v);
		cross(v, w, u);
		x = head[0] + d * (u[0] - w[0]);
		y = head[1] + d * (u[1] - w[1]);
		z = head[2] + d * (u[2] - w[2]);
		glBegin(GL_LINE_STRIP);
		glVertex3fv(head);
		glVertex3f(x, y, z);
		glEnd();
		x = head[0] + d * (-u[0] - w[0]);
		y = head[1] + d * (-u[1] - w[1]);
		z = head[2] + d * (-u[2] - w[2]);
		glBegin(GL_LINE_STRIP);
		glVertex3fv(head);
		glVertex3f(x, y, z);
		glEnd();
	}
	if (axis != Y)
	{
		cross(w, ayy, v);
		(void)unit(v, v);
		cross(v, w, u);
		x = head[0] + d * (u[0] - w[0]);
		y = head[1] + d * (u[1] - w[1]);
		z = head[2] + d * (u[2] - w[2]);
		glBegin(GL_LINE_STRIP);
		glVertex3fv(head);
		glVertex3f(x, y, z);
		glEnd();
		x = head[0] + d * (-u[0] - w[0]);
		y = head[1] + d * (-u[1] - w[1]);
		z = head[2] + d * (-u[2] - w[2]);
		glBegin(GL_LINE_STRIP);
		glVertex3fv(head);
		glVertex3f(x, y, z);
		glEnd();
	}
	if (axis != Z)
	{
		cross(w, azz, v);
		(void)unit(v, v);
		cross(v, w, u);
		x = head[0] + d * (u[0] - w[0]);
		y = head[1] + d * (u[1] - w[1]);
		z = head[2] + d * (u[2] - w[2]);
		glBegin(GL_LINE_STRIP);
		glVertex3fv(head);
		glVertex3f(x, y, z);
		glEnd();
		x = head[0] + d * (-u[0] - w[0]);
		y = head[1] + d * (-u[1] - w[1]);
		z = head[2] + d * (-u[2] - w[2]);
		glBegin(GL_LINE_STRIP);
		glVertex3fv(head);
		glVertex3f(x, y, z);
		glEnd();
	}
	/* done: */
}

///calculate the dot production of two vectors
float dot(float v1[3], float v2[3])
{
	return(v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}

///calculate the cross production of two vectors
void cross(float v1[3], float v2[3], float vout[3])
{
	float tmp[3];
	tmp[0] = v1[1] * v2[2] - v2[1] * v1[2];
	tmp[1] = v2[0] * v1[2] - v1[0] * v2[2];
	tmp[2] = v1[0] * v2[1] - v2[0] * v1[1];
	vout[0] = tmp[0];
	vout[1] = tmp[1];
	vout[2] = tmp[2];
}

///Normalize vector
float unit(float vin[3], float vout[3])
{
	float dist, f;
	dist = vin[0] * vin[0] + vin[1] * vin[1] + vin[2] * vin[2];
	if (dist > 0.0)
	{
		dist = sqrt(dist);
		f = 1. / dist;
		vout[0] = f * vin[0];
		vout[1] = f * vin[1];
		vout[2] = f * vin[2];
	}
	else
	{
		vout[0] = vin[0];
		vout[1] = vin[1];
		vout[2] = vin[2];
	}
	return(dist);
}

void get_vector_field_1(float x, float y, float z, float &vx, float &vy, float &vz, float &magnitude) {
	vx = -3 + 6.*x - 4.*x*(y + 1.) - 4.*z;
	vy = 12.*x - 4.*x*x - 12.*z + 4.*z*z;
	vz = 3. + 4.*x - 4.*x*(y + 1.) - 6.*z + 4.*(y + 1.)*z;
	magnitude = sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2));
}

void get_vector_field_2(float x, float y, float z, float &vx, float &vy, float &vz, float &magnitude) {
	vx = (A * sin(z)) + cos(y);
	vy = (B * sin(x)) + (A * cos(z));
	vz = sin(y) + (B * cos(x));
	magnitude = sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2));
}

void get_vector_field_3(float x, float y, float z, float &vx, float &vy, float &vz, float &magnitude) {
	vx = -1. * y;
	vy = -1. * z;
	vz = x;
	magnitude = sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2));
}

void gen_vec_fields() {
	float ix = 2. / NX3d;
	float iy = 2. / NY3d;
	float iz = 2. / NZ3d;
	for (int i = 0; i < NX3d; i++) {
		for (int j = 0; j < NY3d; j++) {
			for (int k = 0; k < NZ3d; k++) {
				vecNode node;
				node.x = (ix * i) - 1;
				node.y = (iy * j) - 1;
				node.z = (iz * k) - 1;
				get_vector_field_1(node.x, node.y, node.z, node.vx, node.vy, node.vz, node.magnitude);
				vec_field1[i][j][k] = node;
				get_vector_field_2(node.x, node.y, node.z, node.vx, node.vy, node.vz, node.magnitude);
				vec_field2[i][j][k] = node;
				get_vector_field_3(node.x, node.y, node.z, node.vx, node.vy, node.vz, node.magnitude);
				vec_field3[i][j][k] = node;
			}
		}
	}
}

void draw_cube() {
	glColor3f(1., 1., 1.);
	glBegin(GL_LINE_STRIP);
	glVertex3f(-1., -1., 1.);
	glVertex3f(1., -1., 1.);
	glVertex3f(1., 1., 1.);
	glVertex3f(-1., 1., 1.);
	glEnd();
	glBegin(GL_LINE_STRIP);
	glVertex3f(1., -1., 1.);
	glVertex3f(1., -1., -1.);
	glVertex3f(1., 1., -1.);
	glVertex3f(1., 1., 1.);
	glEnd();
	glBegin(GL_LINE_STRIP);
	glVertex3f(1., -1., -1.);
	glVertex3f(1., 1., -1.);
	glVertex3f(-1., 1., -1.);
	glVertex3f(-1., -1., -1.);
	glEnd();
	glBegin(GL_LINE_STRIP);
	glVertex3f(-1., -1., -1.);
	glVertex3f(-1., 1., -1.);
	glVertex3f(-1., 1., 1.);
	glVertex3f(-1., -1., 1.);
	glEnd();
	glBegin(GL_LINE_STRIP);
	glVertex3f(-1., 1., 1.);
	glVertex3f(1., 1., 1.);
	glVertex3f(1., 1., -1.);
	glVertex3f(-1., 1., -1.);
	glEnd();
	glBegin(GL_LINE_STRIP);
	glVertex3f(-1., -1., 1.);
	glVertex3f(1., -1., 1.);
	glVertex3f(1., -1., -1.);
	glVertex3f(-1., -1., -1.);
	glVertex3f(-1., -1., 1.);
	glEnd();
}

// Routine to set a quaternion from a rotation axis and angle
// ( input axis = float[3] angle = float  output: quat = float[4] )
void SetQuaternionFromAxisAngle(const float *axis, float angle, float *quat)
{
	float sina2, norm;
	sina2 = (float)sin(0.5f * angle);
	norm = (float)sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
	quat[0] = sina2 * axis[0] / norm;
	quat[1] = sina2 * axis[1] / norm;
	quat[2] = sina2 * axis[2] / norm;
	quat[3] = (float)cos(0.5f * angle);
}

void draw_3d_arrows_field1() {
	for (int i = 0; i < NX3d; i++) {
		for (int j = 0; j < NY3d; j++) {
			for (int k = 0; k < NZ3d; k++) {
				vecNode *node = &vec_field1[i][j][k];
				glPushMatrix();
				float x = node->x;
				float y = node->y;
				float z = node->z;
				float vx = node->vx;
				float vy = node->vy;
				float vz = node->vz;
				float magnitude = node->magnitude;
				vx /= magnitude;
				vy /= magnitude;
				vz /= magnitude;
				float arrow_head[3] = { x, y, z };
				float arrow_direct[3] = { vx, vy, vz };
				float rgb[3];
				colorFunction(magnitude, rgb);
				glColor3fv(rgb);
				glTranslatef(x, y, z);
				glScalef(0.1, 0.1, 0.1);
				Arrow(arrow_head, arrow_direct);
				glScalef(1., 1., 1.);
				glPopMatrix();
			}
		}
	}
}

void draw_3d_arrows_field2() {
	for (int i = 0; i < NX3d; i++) {
		for (int j = 0; j < NY3d; j++) {
			for (int k = 0; k < NZ3d; k++) {
				vecNode *node = &vec_field2[i][j][k];
				glPushMatrix();
				float x = node->x;
				float y = node->y;
				float z = node->z;
				float vx = node->vx;
				float vy = node->vy;
				float vz = node->vz;
				float magnitude = node->magnitude;
				vx /= magnitude;
				vy /= magnitude;
				vz /= magnitude;
				float arrow_head[3] = { x, y, z };
				float arrow_direct[3] = { vx, vy, vz };
				float rgb[3];
				colorFunction(magnitude, rgb);
				glColor3fv(rgb);
				glTranslatef(x, y, z);
				glScalef(0.1, 0.1, 0.1);
				Arrow(arrow_head, arrow_direct);
				glPopMatrix();
			}
		}
	}
}

void draw_3d_arrows_field3() {
	for (int i = 0; i < NX3d; i++) {
		for (int j = 0; j < NY3d; j++) {
			for (int k = 0; k < NZ3d; k++) {
				vecNode *node = &vec_field3[i][j][k];
				glPushMatrix();
				float x = node->x;
				float y = node->y;
				float z = node->z;
				float vx = node->vx;
				float vy = node->vy;
				float vz = node->vz;
				float magnitude = node->magnitude;
				vx /= magnitude;
				vy /= magnitude;
				vz /= magnitude;
				float arrow_head[3] = { x, y, z };
				float arrow_direct[3] = { vx, vy, vz };
				float rgb[3];
				colorFunction(magnitude, rgb);
				glColor3fv(rgb);
				glTranslatef(x, y, z);
				glScalef(0.1, 0.1, 0.1);
				Arrow(arrow_head, arrow_direct);
				glPopMatrix();
			}
		}
	}
}


// Routine to convert a quaternion to a 4x4 matrix
// ( input: quat = float[4]  output: mat = float[4*4] )
void ConvertQuaternionToMatrix(const float *quat, float *mat)
{
	float yy2 = 2.0f * quat[1] * quat[1];
	float xy2 = 2.0f * quat[0] * quat[1];
	float xz2 = 2.0f * quat[0] * quat[2];
	float yz2 = 2.0f * quat[1] * quat[2];
	float zz2 = 2.0f * quat[2] * quat[2];
	float wz2 = 2.0f * quat[3] * quat[2];
	float wy2 = 2.0f * quat[3] * quat[1];
	float wx2 = 2.0f * quat[3] * quat[0];
	float xx2 = 2.0f * quat[0] * quat[0];
	mat[0 * 4 + 0] = -yy2 - zz2 + 1.0f;
	mat[0 * 4 + 1] = xy2 + wz2;
	mat[0 * 4 + 2] = xz2 - wy2;
	mat[0 * 4 + 3] = 0;
	mat[1 * 4 + 0] = xy2 - wz2;
	mat[1 * 4 + 1] = -xx2 - zz2 + 1.0f;
	mat[1 * 4 + 2] = yz2 + wx2;
	mat[1 * 4 + 3] = 0;
	mat[2 * 4 + 0] = xz2 + wy2;
	mat[2 * 4 + 1] = yz2 - wx2;
	mat[2 * 4 + 2] = -xx2 - yy2 + 1.0f;
	mat[2 * 4 + 3] = 0;
	mat[3 * 4 + 0] = mat[3 * 4 + 1] = mat[3 * 4 + 2] = 0;
	mat[3 * 4 + 3] = 1;
}


// Routine to multiply 2 quaternions (ie, compose rotations)
// ( input q1 = float[4] q2 = float[4]  output: qout = float[4] )
void MultiplyQuaternions(const float *q1, const float *q2, float *qout)
{
	float qr[4];
	qr[0] = q1[3] * q2[0] + q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1];
	qr[1] = q1[3] * q2[1] + q1[1] * q2[3] + q1[2] * q2[0] - q1[0] * q2[2];
	qr[2] = q1[3] * q2[2] + q1[2] * q2[3] + q1[0] * q2[1] - q1[1] * q2[0];
	qr[3] = q1[3] * q2[3] - (q1[0] * q2[0] + q1[1] * q2[1] + q1[2] * q2[2]);
	qout[0] = qr[0]; qout[1] = qr[1]; qout[2] = qr[2]; qout[3] = qr[3];
}

// Color space conversion function
// HSV to RGB
void
HsvRgb(float hsv[3], float rgb[3])
{
	float h, s, v;			// hue, sat, value
	float r, g, b;			// red, green, blue
	float i, f, p, q, t;		// interim values


								// guarantee valid input:

	h = hsv[0] / 60.;
	while (h >= 6.)	h -= 6.;
	while (h < 0.) 	h += 6.;

	s = hsv[1];
	if (s < 0.)
		s = 0.;
	if (s > 1.)
		s = 1.;

	v = hsv[2];
	if (v < 0.)
		v = 0.;
	if (v > 1.)
		v = 1.;


	// if sat==0, then is a gray:

	if (s == 0.0)
	{
		rgb[0] = rgb[1] = rgb[2] = v;
		return;
	}


	// get an rgb from the hue itself:

	i = floor(h);
	f = h - i;
	p = v * (1. - s);
	q = v * (1. - s * f);
	t = v * (1. - (s * (1. - f)));

	switch ((int)i)
	{
	case 0:
		r = v;	g = t;	b = p;
		break;

	case 1:
		r = q;	g = v;	b = p;
		break;

	case 2:
		r = p;	g = v;	b = t;
		break;

	case 3:
		r = p;	g = q;	b = v;
		break;

	case 4:
		r = t;	g = p;	b = v;
		break;

	case 5:
		r = v;	g = p;	b = q;
		break;
	}


	rgb[0] = r;
	rgb[1] = g;
	rgb[2] = b;
}

// Return elapsed time in milliseconds
int GetTimeMs()
{
#if !defined(_WIN32)
	return glutGet(GLUT_ELAPSED_TIME);
#else
	// glutGet(GLUT_ELAPSED_TIME) seems buggy on Windows
	return (int)GetTickCount();
#endif
}

//
//	Draw a set of 3D axes:
//	(length is the axis length in world coordinates)
//

void Axes(float length)
{
	int i, j;			// counters
	float fact;			// character scale factor
	float base;			// character start location

	glEnable(GL_COLOR_MATERIAL);
	glBegin(GL_LINE_STRIP);
	glColor3f(1, 0, 0);
	glVertex3f(length, 0., 0.);
	glVertex3f(0., 0., 0.);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glColor3f(0, 1, 0);
	glVertex3f(0., 0., 0.);
	glVertex3f(0., length, 0.);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glColor3f(0, 0, 1);
	glVertex3f(0., 0., 0.);
	glVertex3f(0., 0., length);
	glEnd();

	fact = LENFRAC * length;
	base = BASEFRAC * length;

	glColor3f(1, 0, 0);
	glBegin(GL_LINE_STRIP);
	for (i = 0; i < 4; i++)
	{
		j = xorder[i];
		if (j < 0)
		{

			glEnd();
			glBegin(GL_LINE_STRIP);
			j = -j;
		}
		j--;
		glVertex3f(base + fact * xx[j], fact*xy[j], 0.0);
	}
	glEnd();

	glColor3f(0, 1, 0);
	glBegin(GL_LINE_STRIP);
	for (i = 0; i < 5; i++)
	{
		j = yorder[i];
		if (j < 0)
		{

			glEnd();
			glBegin(GL_LINE_STRIP);
			j = -j;
		}
		j--;
		glVertex3f(fact*yx[j], base + fact * yy[j], 0.0);
	}
	glEnd();

	glColor3f(0, 0, 1);
	glBegin(GL_LINE_STRIP);
	for (i = 0; i < 6; i++)
	{
		j = zorder[i];
		if (j < 0)
		{

			glEnd();
			glBegin(GL_LINE_STRIP);
			j = -j;
		}
		j--;
		glVertex3f(0.0, fact*zy[j], base + fact * zx[j]);
	}
	glEnd();

}

void InitAxesLists(void)
{
	float dx = BOXSIZE / 2.;
	float dy = BOXSIZE / 2.;
	float dz = BOXSIZE / 2.;

	// create the object:

	//BoxList = glGenLists(1);
	glNewList(BoxList, GL_COMPILE);

	glBegin(GL_QUADS);

	glColor3f(0., 0., 1.);
	glNormal3f(0., 0., 1.);
	glVertex3f(-dx, -dy, dz);
	glVertex3f(dx, -dy, dz);
	glVertex3f(dx, dy, dz);
	glVertex3f(-dx, dy, dz);

	glNormal3f(0., 0., -1.);
	glTexCoord2f(0., 0.);
	glVertex3f(-dx, -dy, -dz);
	glTexCoord2f(0., 1.);
	glVertex3f(-dx, dy, -dz);
	glTexCoord2f(1., 1.);
	glVertex3f(dx, dy, -dz);
	glTexCoord2f(1., 0.);
	glVertex3f(dx, -dy, -dz);

	glColor3f(1., 0., 0.);
	glNormal3f(1., 0., 0.);
	glVertex3f(dx, -dy, dz);
	glVertex3f(dx, -dy, -dz);
	glVertex3f(dx, dy, -dz);
	glVertex3f(dx, dy, dz);

	glNormal3f(-1., 0., 0.);
	glVertex3f(-dx, -dy, dz);
	glVertex3f(-dx, dy, dz);
	glVertex3f(-dx, dy, -dz);
	glVertex3f(-dx, -dy, -dz);

	glColor3f(0., 1., 0.);
	glNormal3f(0., 1., 0.);
	glVertex3f(-dx, dy, dz);
	glVertex3f(dx, dy, dz);
	glVertex3f(dx, dy, -dz);
	glVertex3f(-dx, dy, -dz);

	glNormal3f(0., -1., 0.);
	glVertex3f(-dx, -dy, dz);
	glVertex3f(-dx, -dy, -dz);
	glVertex3f(dx, -dy, -dz);
	glVertex3f(dx, -dy, dz);

	glEnd();

	glEndList();


	// create the axes:

	//AxesList = glGenLists(1);
	glNewList(AxesList, GL_COMPILE);
	glLineWidth(AXES_WIDTH);
	Axes(1.5);
	glLineWidth(1.);
	glEndList();
}

void Rainbow_color(float s, float rgb[3])
{
	float t = (s - *min_ptr) / (*max_ptr - *min_ptr);
	// make sure t is between 0 and 1, if not, rgb should be black
	if (t < 0 || t>1) {
		rgb[0] = rgb[1] = rgb[2] = 0.;
		return;
	}
	float hsv[3] = { 1. };
	// map the scalar value linearly to the hue channel of the HSV
	hsv[0] = (1.0 - t) * 240;
	hsv[1] = hsv[2] = 1.; // set the saturation and value as 1
	// Call the HSV to RGB conversion function
	HsvRgb(hsv, rgb);
}

void BWR_Divergent(float s, float rgb[3]) {
	float t = (s - *min_ptr) / (*max_ptr - *min_ptr);
	float hsv[4];
	hsv[2] = hsv[3] = 1;
	if (t <= g_WhiteThreshold) {
		hsv[0] = 240;
		hsv[1] = 1 - ((1 / g_WhiteThreshold) * t);
	}
	else {
		hsv[0] = 0;
		hsv[1] = ((1 / g_WhiteThreshold) * t) - 1;
	}
	HsvRgb(hsv, rgb);
}

void HeatMap(float s, float rgb[3]) {
	float t = (s - *min_ptr) / (*max_ptr - *min_ptr);
	if (t <= 0) {
		rgb[0] = rgb[1] = rgb[2] = 0.; //This is the coldest hence black
		return;
	}
	if (t >= 1) {
		rgb[0] = rgb[1] = rgb[2] = 1.; // This is the hottest hence white
		return;
	}
	// We follow the sequential R + G + B scale to get the heatmap color scheme.
	// Since red needs to be full first before going to green, we therefore split the normalized scalar space into 3 equal parts,
	// if the value is in the 1st part, then we know that we only need a shade of red. The intensity of red will depend on the 
	// value of the scalar space.
	// Since the space has been divided into 3 parts, no color value can be bigger than (1./3.) therefore we multiply the value 
	// by 3 to scale it back to the normal color range of 0 - 1. 
	// If a color lies in the 2nd part, then we maximize red and carryforward the spillover to set green.
	// We do similarly for blue
	rgb[0] = min((3 * max(t, 0)), 1); // 3 * (min(t, 1/3))
	rgb[1] = min((3 * max(t - (1. / 3.), 0)), 1); // 3 * (min(t-1/3, 1/3))
	rgb[2] = min((3 * max(t - (2. / 3.), 0)), 1); // 3 * (min(t-2/3, 1/3))
}

void Discrete(float s, float rgb[3]) {
	int t = floor((s - *min_ptr) / (*max_ptr - *min_ptr) * 10);
	if (t >= 0 && t <= 5) {
		rgb[0] = Colors[t][0];
		rgb[1] = Colors[t][1];
		rgb[2] = Colors[t][2];
		return;
	}
	if (t < 0) {
		rgb[0] = rgb[1] = rgb[2] = 0.;
		return;
	}
	rgb[0] = rgb[1] = rgb[2] = 1.;
}

void NonLinear(float s, float rgb[3]) {
	float t = (sqrt(s) - sqrt(*min_ptr)) / (sqrt(*max_ptr) - sqrt(*min_ptr));
	float hsv[4];
	hsv[1] = hsv[2] = hsv[3] = 1.;
	hsv[0] = (1 - t) * 240;
	HsvRgb(hsv, rgb);
}

void setColorFunction() {
	switch (whichColor) {
	case 0:
		colorFunction = &Rainbow_color;
		break;
	case 1:
		colorFunction = &BWR_Divergent;
		break;
	case 2:
		colorFunction = &HeatMap;
		break;
	case 3:
		colorFunction = &Discrete;
		break;
	case 4:
		colorFunction = &NonLinear;
		break;
	}
}

void draw_streamlines() {
	float rgb[3];
	for (int i = 0; i < (int) streamlines.size(); i++) {
		std::vector<streampoint> streamline = streamlines[i];
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j < streamline.size(); j++) {
			streampoint *point = &streamline[j];
			colorFunction(point->magnitude, rgb);
			glColor3fv(rgb);
			glVertex3f(point->nextX, point->nextY, point->nextZ);
		}
		glEnd();
	}
}

void draw_streamribbon() {
	float rgb[3];
	for (int i = 0; i < (int) streamr1.size() - 1; i++) {
		glBegin(GL_QUADS);
		streampoint *p1 = &streamr1[i];
		streampoint *p2 = &streamr2[i];
		colorFunction(p1->magnitude, rgb);
		glColor3fv(rgb);
		glVertex3f(p1->nextX, p1->nextY, p1->nextZ);
		colorFunction(p2->magnitude, rgb);
		glColor3fv(rgb);
		glVertex3f(p2->nextX, p2->nextY, p2->nextZ);
		++p1;
		++p2;
		colorFunction(p2->magnitude, rgb);
		glColor3fv(rgb);
		glVertex3f(p2->nextX, p2->nextY, p2->nextZ);
		colorFunction(p1->magnitude, rgb);
		glColor3fv(rgb);
		glVertex3f(p1->nextX, p1->nextY, p1->nextZ);
		glEnd();
	}
}

// Callback function called by GLUT to render screen
void Display(void)
{
	float v[4]; // will be used to set light parameters
	float mat[4 * 4]; // rotation matrix

	// Clear frame buffer
	glClearColor(0, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glEnable(GL_NORMALIZE);

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	//glShadeModel(GL_FLAT);

	// Set light
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	v[0] = v[1] = v[2] = g_LightMultiplier * 0.4f; v[3] = 1.0f;
	glLightfv(GL_LIGHT0, GL_AMBIENT, v);
	v[0] = v[1] = v[2] = g_LightMultiplier * 0.8f; v[3] = 1.0f;
	glLightfv(GL_LIGHT0, GL_DIFFUSE, v);
	v[0] = -g_LightDirection[0]; v[1] = -g_LightDirection[1]; v[2] = -g_LightDirection[2]; v[3] = 0.0f;
	glLightfv(GL_LIGHT0, GL_POSITION, v);

	// Set material
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, g_MatAmbient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, g_MatDiffuse);

	// Rotate and draw shape
	glPushMatrix();
	glTranslatef(0.5f, -0.3f, 0.0f);
	if (g_AutoRotate)
	{
		float axis[3] = { 0, 1, 0 };
		float angle = (float)(GetTimeMs() - g_RotateTime) / 1000.0f;
		float quat[4];
		SetQuaternionFromAxisAngle(axis, angle, quat);
		MultiplyQuaternions(g_RotateStart, quat, g_Rotation);
	}
	ConvertQuaternionToMatrix(g_Rotation, mat);
	glMultMatrixf(mat);
	glScalef(g_Zoom, g_Zoom, g_Zoom);

	//glCallList(g_CurrentShape);

	setColorFunction();
	draw_cube();

	if (g_arrows) {
		if (currentField == 1) {
			max_ptr = &max_sv1;
			min_ptr = &min_sv1;
			draw_3d_arrows_field1();
		}
		else if (currentField == 2) {
			max_ptr = &max_sv2;
			min_ptr = &min_sv2;
			draw_3d_arrows_field2();			
		}
		else {
			max_ptr = &max_sv3;
			min_ptr = &min_sv3;
			draw_3d_arrows_field3();
		}
	}

	if (g_streamlines) {
		max_ptr = &abs_s_max;
		min_ptr = &abs_s_min;
		draw_streamlines();
	}

	if (g_streamribbon) {
		max_ptr = &abs_s_max;
		min_ptr = &abs_s_min;
		draw_streamribbon();
	}

		// Draw axes
	if (g_Axes)
		glCallList(AxesList);


	glPopMatrix();

	// Draw tweak bars
	TwDraw();

	// Present frame buffer
	glutSwapBuffers();

	// Recall Display at next frame
	glutPostRedisplay();
}

// Callback function called by GLUT when window size changes
void Reshape(int width, int height)
{
	// Set OpenGL viewport and camera
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(40, (double)width / height, 1, 10);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);
	glTranslatef(0, 0.6f, -1);

	// Send the new window size to AntTweakBar
	TwWindowSize(width, height);
}


// Function called at exit
void Terminate(void)
{
	//glDeleteLists(SHAPE_TEAPOT, NUM_SHAPES);

	TwTerminate();
}


//  Callback function called when the 'AutoRotate' variable value of the tweak bar has changed
void TW_CALL SetAutoRotateCB(const void *value, void *clientData)
{
	(void)clientData; // unused

	g_AutoRotate = *(const int *)value; // copy value to g_AutoRotate
	if (g_AutoRotate != 0)
	{
		// init rotation
		g_RotateTime = GetTimeMs();
		g_RotateStart[0] = g_Rotation[0];
		g_RotateStart[1] = g_Rotation[1];
		g_RotateStart[2] = g_Rotation[2];
		g_RotateStart[3] = g_Rotation[3];

		// make Rotation variable read-only
		TwDefine(" TweakBar/ObjRotation readonly ");
	}
	else
		// make Rotation variable read-write
		TwDefine(" TweakBar/ObjRotation readwrite ");
}


//  Callback function called by the tweak bar to get the 'AutoRotate' value
void TW_CALL GetAutoRotateCB(void *value, void *clientData)
{
	(void)clientData; // unused
	*(int *)value = g_AutoRotate; // copy g_AutoRotate to value
}


//  Callback function called when the 'Axes' variable value of the tweak bar has changed
void TW_CALL SetAxesCB(const void *value, void *clientData)
{
	(void)clientData; // unused

	g_Axes = *(const int *)value; // copy value to g_Axes
}

//  Callback function called by the tweak bar to get the 'Axes' value
void TW_CALL GetAxesCB(void *value, void *clientData)
{
	(void)clientData; // unused
	*(int *)value = g_Axes; // copy g_AutoRotate to value
}

double getMagnitude(float x, float y, float z) {
	float magnitude, vx, vy, vz;
	vectorFunction(x, y, z, vx, vy, vz, magnitude);
	return magnitude;
}

double edist(float x1, float x2, float y1, float y2, float z1, float z2) {
	return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2));
}

double getDistance(float x, float y, float z) {
	float distance = 1.;
	// return shortest distance of this point from any all points in the streamline
	for (int i = 0; i < (int)streamline.size(); i++) {
		streampoint *p = &streamline[i];
		float temp = edist(x, p->nextX, y, p->nextY, z, p->nextZ);
		if (temp < distance)
			distance = temp;
	}
	return distance;
}

double getSeparation(float x, float y, float z) {
	float separation = 1.;
	// return the shortest distance of this point on this streamline from the distance on all other streamlines
	// we don't care which streamline it's closest to, we just want to terminate computation if it is too close to any of the other
	for (int i = 0; i < (int)streamlines.size(); i++) {
		std::vector<streampoint> stream = streamlines[i];
		for (int j = 0; j < (int)stream.size(); j++) {
			streampoint *p = &stream[j];	
			float temp = edist(x, p->nextX, y, p->nextY, z, p->nextZ);
			if (temp < separation)
				separation = temp;
		}
	}
	return separation;
}

bool checkConditions(float next_x, float next_y, float next_z, int n_steps, bool computeRibbon) {
	bool conditions, condition1, condition2, condition3, condition4, condition5;	
	// next value should be within bounds
	condition1 = (next_x >= -1. && next_x <= 1.) && (next_y >= -1. && next_y <= 1.) && (next_z >= -1. && next_z <= 1.);
	// next value is not a fixed point
	condition2 = getMagnitude(next_x, next_y, next_z) > 0.000001;
	// streamline isn't curving back on itself - check if the euclidean distance is greater than threshold
	// condition3 = getDistance(next_x, next_y, next_z) > 0.001;
	// streamline length in steps has been reached
	condition4 = n_steps <= g_streamLength;
	// streamline is too close to other streamlines - don't need to check this for streamribbons
	// condition5 = getSeparation(next_x, next_y, next_z) > 0.05 || computeRibbon;
	// check if all conditions are satisfied
	conditions = condition1 && condition2 && condition4;
	
	return conditions;
}

bool computeStreamline(float seedX, float seedY, float seedZ) {
	int n_steps = 0;
	float next_x, next_y, next_z;
	next_x = seedX;
	next_y = seedY;
	next_z = seedZ;
	float vx, vy, vz, magnitude;
	streamline.clear();
	streampoint point;
	if (whichIntegrator == 0) {
		// Euler
		// Sn = Sn-1 + v(Sn-1).dt
		while (checkConditions(next_x, next_y, next_z, n_steps, false)) {
			vectorFunction(next_x, next_y, next_z, vx, vy, vz, magnitude);
			point.nextX = next_x;
			point.nextY = next_y;
			point.nextZ = next_z;
			point.magnitude = magnitude;
			if (abs_s_max < magnitude)
				abs_s_max = magnitude;
			if (abs_s_min > magnitude)
				abs_s_min = magnitude;
			streamline.push_back(point);
			n_steps++;
			vx /= magnitude;
			vy /= magnitude;
			vz /= magnitude;
			vx *= g_step;
			vy *= g_step;
			vz *= g_step;
			next_x += vx;
			next_y += vy;
			next_z += vz;
		}
	}
	else {
		// RK - 2
		float vx2, vy2, vz2, vxavg, vyavg, vzavg;
		while (checkConditions(next_x, next_y, next_z, n_steps, false)) {
			vectorFunction(next_x, next_y, next_z, vx, vy, vz, magnitude);
			point.nextX = next_x;
			point.nextY = next_y;
			point.nextZ = next_z;
			point.magnitude = magnitude;
			if (abs_s_max < magnitude)
				abs_s_max = magnitude;
			if (abs_s_min > magnitude)
				abs_s_min = magnitude;
			streamline.push_back(point);
			n_steps++;
			vx /= magnitude;
			vy /= magnitude;
			vz /= magnitude;
			vx *= g_step;
			vy *= g_step;
			vz *= g_step;
			next_x += vx;
			next_y += vy;
			next_z += vz;
			vectorFunction(next_x, next_y, next_z, vx2, vy2, vz2, magnitude);
			vx2 /= magnitude;
			vy2 /= magnitude;
			vz2 /= magnitude;
			vx2 *= g_step;
			vy2 *= g_step;
			vz2 *= g_step;
			vxavg = (vx + vx2) / 2;
			vyavg = (vy + vy2) / 2;
			vzavg = (vz + vz2) / 2;
			next_x += vxavg - vx;
			next_y += vyavg - vy;
			next_z += vzavg - vz;
		}
	}
	if (n_steps == 0)
		return false;
	streamlines.push_back(streamline);
	return true;
}

void computeStreamBunch() {
	streamlines.clear();
	if (!computeStreamline(g_probeX, g_probeY, g_probeZ)) // the probe - to check if the probe point is fixed
		return;
	computeStreamline(g_probeX + g_linedist, g_probeY, g_probeZ); // first streamline
	float x_coord, y_coord, theta;
	for (int i = 40; i <= 160; i += 40) {
		theta = i * M_PI / 180;
		x_coord = sqrt(pow(g_linedist, 2) / (1 + pow(tan(theta), 2))) * (i < 90 ? 1 : -1);
		y_coord = x_coord * tan(theta);
		computeStreamline(g_probeX + x_coord, g_probeY + y_coord, g_probeZ);
		computeStreamline(g_probeX + x_coord, g_probeY - y_coord, g_probeZ);
	}
}

void computeStreamribbon() {
	int n_steps = 0;
	float next_x1, next_y1, next_z1, next_x2, next_y2, next_z2;
	next_x1 = g_probeX;
	next_y1 = g_probeY;
	next_z1 = g_probeZ;
	next_x2 = g_probeX + g_rdist;
	next_y2 = g_probeY;
	next_z2 = g_probeZ;
	float vx, vy, vz, magnitude;
	streamr1.clear();
	streamr2.clear();
	streampoint point1;
	streampoint point2;
	if (whichIntegrator == 0) {
		// Euler
		// Sn = Sn-1 + v(Sn-1).dt
		while (checkConditions(next_x1, next_y1, next_z1, n_steps, true) && checkConditions(next_x2, next_y2, next_z2, n_steps, true)) {
			n_steps++;
			vectorFunction(next_x1, next_y1, next_z1, vx, vy, vz, magnitude);
			point1.nextX = next_x1;
			point1.nextY = next_y1;
			point1.nextZ = next_z1;
			point1.magnitude = magnitude;
			if (abs_s_max < magnitude)
				abs_s_max = magnitude;
			if (abs_s_min > magnitude)
				abs_s_min = magnitude;
			streamr1.push_back(point1);
			vx /= magnitude;
			vy /= magnitude;
			vz /= magnitude;
			vx *= g_step;
			vy *= g_step;
			vz *= g_step;
			next_x1 += vx;
			next_y1 += vy;
			next_z1 += vz;
			vectorFunction(next_x2, next_y2, next_z2, vx, vy, vz, magnitude);
			point2.nextX = next_x2;
			point2.nextY = next_y2;
			point2.nextZ = next_z2;
			point2.magnitude = magnitude;
			if (abs_s_max < magnitude)
				abs_s_max = magnitude;
			if (abs_s_min > magnitude)
				abs_s_min = magnitude;
			streamr2.push_back(point2);
			vx /= magnitude;
			vy /= magnitude;
			vz /= magnitude;
			vx *= g_step;
			vy *= g_step;
			vz *= g_step;
			next_x2 += vx;
			next_y2 += vy;
			next_z2 += vz;
			float dist = edist(next_x1, next_x2, next_y1, next_y2, next_z1, next_z2);
			next_x2 = next_x1 + ((g_rdist * (next_x2 - next_x1)) / dist);
			next_y2 = next_y1 + ((g_rdist * (next_y2 - next_y1)) / dist);
			next_z2 = next_z1 + ((g_rdist * (next_z2 - next_z1)) / dist);
		}
	}
	else {
		// RK - 2 
		float vx21, vy21, vz21, vx22, vy22, vz22, vxavg, vyavg, vzavg;
		while (checkConditions(next_x1, next_y1, next_z1, n_steps, true) && checkConditions(next_x2, next_y2, next_z2, n_steps, true)) {
			n_steps++;
			vectorFunction(next_x1, next_y1, next_z1, vx, vy, vz, magnitude);
			point1.nextX = next_x1;
			point1.nextY = next_y1;
			point1.nextZ = next_z1;
			point1.magnitude = magnitude;
			if (abs_s_max < magnitude)
				abs_s_max = magnitude;
			if (abs_s_min > magnitude)
				abs_s_min = magnitude;
			streamr1.push_back(point1);
			vx /= magnitude;
			vy /= magnitude;
			vz /= magnitude;
			vx *= g_step;
			vy *= g_step;
			vz *= g_step;
			next_x1 += vx;
			next_y1 += vy;
			next_z1 += vz;
			vectorFunction(next_x1, next_y1, next_z1, vx21, vy21, vz21, magnitude);
			vx21 /= magnitude;
			vy21 /= magnitude;
			vz21 /= magnitude;
			vx21 *= g_step;
			vy21 *= g_step;
			vz21 *= g_step;
			vxavg = (vx + vx21) / 2;
			vyavg = (vy + vy21) / 2;
			vzavg = (vz + vz21) / 2;
			next_x1 += vxavg - vx;
			next_y1 += vyavg - vy;
			next_z1 += vzavg - vz;
			vectorFunction(next_x2, next_y2, next_z2, vx, vy, vz, magnitude);
			point2.nextX = next_x2;
			point2.nextY = next_y2;
			point2.nextZ = next_z2;
			point2.magnitude = magnitude;
			if (abs_s_max < magnitude)
				abs_s_max = magnitude;
			if (abs_s_min > magnitude)
				abs_s_min = magnitude;
			streamr2.push_back(point2);
			vx /= magnitude;
			vy /= magnitude;
			vz /= magnitude;
			vx *= g_step;
			vy *= g_step;
			vz *= g_step;
			next_x2 += vx;
			next_y2 += vy;
			next_z2 += vz;
			vectorFunction(next_x2, next_y2, next_z2, vx22, vy22, vz22, magnitude);
			vx22 /= magnitude;
			vy22 /= magnitude;
			vz22 /= magnitude;
			vx22 *= g_step;
			vy22 *= g_step;
			vz22 *= g_step;
			vxavg = (vx + vx22) / 2;
			vyavg = (vy + vy22) / 2;
			vzavg = (vz + vz22) / 2;
			next_x2 += vxavg - vx;
			next_y2 += vyavg - vy;
			next_z2 += vzavg - vz;
			float dist = edist(next_x1, next_x2, next_y1, next_y2, next_z1, next_z2);
			next_x2 = next_x1 + ((g_rdist * (next_x2 - next_x1)) / dist);
			next_y2 = next_y1 + ((g_rdist * (next_y2 - next_y1)) / dist);
			next_z2 = next_z1 + ((g_rdist * (next_z2 - next_z1)) / dist);
		}
	}
}

void TW_CALL setArrowCB(const void* value, void* clientData) {
	g_arrows = *(const int *)value;
	g_streamribbon = g_streamlines = 0;
}

void TW_CALL getArrowCB(void* value, void* clientData) {
	*(int *)value = g_arrows;
}

void TW_CALL setStreamlinesCB(const void* value, void* clientData) {
	g_streamlines = *(const int *)value;
	// g_streamribbon = 0;
	g_arrows = 0;
}

void TW_CALL getStreamlinesCB(void* value, void* clientData) {
	*(int *)value = g_streamlines;
}

void TW_CALL setStreamribbonCB(const void* value, void* clientData) {
	g_streamribbon = *(const int *)value;
	// g_streamlines = 0;
	g_arrows = 0;
}

void TW_CALL getStreamribbonCB(void* value, void* clientData) {
	*(int *)value = g_streamribbon;
}

void setVecFieldPointer() {
	if (currentField == 1)
		vectorFunction = &get_vector_field_1;
	else if (currentField == 2)
		vectorFunction = &get_vector_field_2;
	else if (currentField == 3)
		vectorFunction = &get_vector_field_3;
}

void TW_CALL loadNewObjCB(void *clientData)
{
	switch (g_CurrentShape) {
	case 0:
		currentField = 1;
		max_ptr = &max_sv1;
		min_ptr = &min_sv1;
		break;

	case 1:
		currentField = 2;
		max_ptr = &max_sv2;
		min_ptr = &min_sv2;
		break;

	case 2:
		currentField = 3;
		max_ptr = &max_sv3;
		min_ptr = &min_sv3;
		break;
	}
	abs_s_max = abs_s_min = 0.;
	setVecFieldPointer();
	computeStreamBunch();
	computeStreamribbon();

	g_WhiteThreshold = 0.5; // reset g_WhiteThreshold
	glutSetWindow(MainWindow);
	glutPostRedisplay();
}

void calcLimits3dVec() {
	max_sv1 = min_sv1 = vec_field1[0][0][0].magnitude;
	max_sv2 = min_sv2 = vec_field2[0][0][0].magnitude;
	max_sv3 = min_sv3 = vec_field3[0][0][0].magnitude;
	max_vx1 = min_vx1 = vec_field1[0][0][0].vx;
	max_vx2 = min_vx2 = vec_field2[0][0][0].vx;
	max_vx3 = min_vx3 = vec_field3[0][0][0].vx;
	max_vy1 = min_vy1 = vec_field1[0][0][0].vy;
	max_vy2 = min_vy2 = vec_field2[0][0][0].vy;
	max_vy3 = min_vy3 = vec_field3[0][0][0].vy;
	max_vz1 = min_vz1 = vec_field1[0][0][0].vz;
	max_vz2 = min_vz2 = vec_field2[0][0][0].vz;
	max_vz3 = min_vz3 = vec_field3[0][0][0].vz;
	for (int i = 0; i < NX3d; i++) {
		for (int j = 0; j < NY3d; j++) {
			for (int k = 0; k < NZ3d; k++) {
				if (vec_field1[i][j][k].magnitude > max_sv1)
					max_sv1 = vec_field1[i][j][k].magnitude;
				if (vec_field1[i][j][k].magnitude < min_sv1)
					min_sv1 = vec_field1[i][j][k].magnitude;
				if (vec_field2[i][j][k].magnitude > max_sv2)
					max_sv2 = vec_field2[i][j][k].magnitude;
				if (vec_field2[i][j][k].magnitude < min_sv2)
					min_sv2 = vec_field2[i][j][k].magnitude;
				if (vec_field3[i][j][k].magnitude > max_sv3)
					max_sv3 = vec_field3[i][j][k].magnitude;
				if (vec_field3[i][j][k].magnitude < min_sv3)
					min_sv3 = vec_field3[i][j][k].magnitude;

				if (vec_field1[i][j][k].magnitude > max_vx1)
					max_vx1 = vec_field1[i][j][k].magnitude;
				if (vec_field1[i][j][k].magnitude < min_vx1)
					min_vx1 = vec_field1[i][j][k].magnitude;
				if (vec_field2[i][j][k].magnitude > max_vx2)
					max_vx2 = vec_field2[i][j][k].magnitude;
				if (vec_field2[i][j][k].magnitude < min_vx2)
					min_vx2 = vec_field2[i][j][k].magnitude;
				if (vec_field3[i][j][k].magnitude > max_vx3)
					max_vx3 = vec_field3[i][j][k].magnitude;
				if (vec_field3[i][j][k].magnitude < min_vx3)
					min_vx3 = vec_field3[i][j][k].magnitude;

				if (vec_field1[i][j][k].magnitude > max_vy1)
					max_vy1 = vec_field1[i][j][k].magnitude;
				if (vec_field1[i][j][k].magnitude < min_vy1)
					min_vy1 = vec_field1[i][j][k].magnitude;
				if (vec_field2[i][j][k].magnitude > max_vy2)
					max_vy2 = vec_field2[i][j][k].magnitude;
				if (vec_field2[i][j][k].magnitude < min_vy2)
					min_vy2 = vec_field2[i][j][k].magnitude;
				if (vec_field3[i][j][k].magnitude > max_vy3)
					max_vy3 = vec_field3[i][j][k].magnitude;
				if (vec_field3[i][j][k].magnitude < min_vy3)
					min_vy3 = vec_field3[i][j][k].magnitude;

				if (vec_field1[i][j][k].magnitude > max_vz1)
					max_vz1 = vec_field1[i][j][k].magnitude;
				if (vec_field1[i][j][k].magnitude < min_vz1)
					min_vz1 = vec_field1[i][j][k].magnitude;
				if (vec_field2[i][j][k].magnitude > max_vz2)
					max_vz2 = vec_field2[i][j][k].magnitude;
				if (vec_field2[i][j][k].magnitude < min_vz2)
					min_vz2 = vec_field2[i][j][k].magnitude;
				if (vec_field3[i][j][k].magnitude > max_vz3)
					max_vz3 = vec_field3[i][j][k].magnitude;
				if (vec_field3[i][j][k].magnitude < min_vz3)
					min_vz3 = vec_field3[i][j][k].magnitude;
			}
		}
	}
}

void recompStreamline(void* clientData) {
	abs_s_max = abs_s_min = 0.;
	computeStreamBunch();
	computeStreamribbon();
}

void InitTwBar()
{
	// Create a tweak bar
	bar = TwNewBar("TweakBar");
	TwDefine(" GLOBAL help='This example shows how to integrate AntTweakBar with GLUT and OpenGL.' "); // Message added to the help bar.
	TwDefine(" TweakBar size='290 650' color='0 128 255' alpha=128 position='0 0' "); // change default tweak bar size and color
	TwDefine(" TweakBar  label='Visual Parameters'");        // change the title of the Tweakbar

	// Add callback to toggle reference axes (callback functions are defined above).
	TwAddVarCB(bar, "Axes", TW_TYPE_BOOL32, SetAxesCB, GetAxesCB, NULL,
		" label='Axes' key=a help='Toggle reference axes.' ");

	TwAddSeparator(bar, NULL, NULL);

	// Add 'g_Zoom' to 'bar': this is a modifable (RW) variable of type TW_TYPE_FLOAT. Its key shortcuts are [z] and [Z].
	TwAddVarRW(bar, "Zoom", TW_TYPE_FLOAT, &g_Zoom,
		" min=0.01 max=10 step=0.01 keyIncr=z keyDecr=Z help='Scale the object (1=original size).' ");

	// Add 'g_Rotation' to 'bar': this is a variable of type TW_TYPE_QUAT4F which defines the object's orientation
	TwAddVarRW(bar, "ObjRotation", TW_TYPE_QUAT4F, &g_Rotation,
		" label='Object rotation' opened=true help='Change the object orientation.' ");

	// Add callback to toggle auto-rotate mode (callback functions are defined above).
	TwAddVarCB(bar, "AutoRotate", TW_TYPE_BOOL32, SetAutoRotateCB, GetAutoRotateCB, NULL,
		" label='Auto-rotate' key=space help='Toggle auto-rotate mode.' ");

	// Add 'g_LightMultiplier' to 'bar': this is a variable of type TW_TYPE_FLOAT. Its key shortcuts are [+] and [-].
	TwAddVarRW(bar, "Multiplier", TW_TYPE_FLOAT, &g_LightMultiplier,
		" label='Light booster' min=0.1 max=4 step=0.02 keyIncr='+' keyDecr='-' help='Increase/decrease the light power.' ");

	// Add 'g_LightDirection' to 'bar': this is a variable of type TW_TYPE_DIR3F which defines the light direction
	TwAddVarRW(bar, "LightDir", TW_TYPE_DIR3F, &g_LightDirection,
		" label='Light direction' opened=true help='Change the light direction.' ");

	// Add 'g_MatAmbient' to 'bar': this is a variable of type TW_TYPE_COLOR3F (3 floats color, alpha is ignored)
	// and is inserted into a group named 'Material'.
	TwAddVarRW(bar, "Ambient", TW_TYPE_COLOR3F, &g_MatAmbient, " group='Material' ");

	// Add 'g_MatDiffuse' to 'bar': this is a variable of type TW_TYPE_COLOR3F (3 floats color, alpha is ignored)
	// and is inserted into group 'Material'.
	TwAddVarRW(bar, "Diffuse", TW_TYPE_COLOR3F, &g_MatDiffuse, " group='Material' ");

	TwAddSeparator(bar, " objects ", NULL);

	// Add the enum variable 'g_CurrentShape' to 'bar'
	// (before adding an enum variable, its enum type must be declared to AntTweakBar as follow)
	{
		TwEnumVal shapeEV[NUM_SHAPES] = { {0, "Field 1"}, {1, "Field 2"}, {2, "Field 3"} };

		// Create a type for the enum shapeEV
		TwType shapeType = TwDefineEnum("ShapeType", shapeEV, NUM_SHAPES);

		// add 'g_CurrentShape' to 'bar': this is a variable of type ShapeType. Its key shortcuts are [<] and [>].
		TwAddVarRW(bar, "Shape", shapeType, &g_CurrentShape, " keyIncr='<' keyDecr='>' help='Change object shape.' ");

		// add a button to reload the selected object
		TwAddButton(bar, "Update (Re-load)", loadNewObjCB, NULL, " label='Load new object after selection' ");
	}

	TwAddSeparator(bar, " others ", NULL);

	// Add the enum variable 'whichColor' to 'bar' 
	{
		TwEnumVal ColorEV[NUM_COLORS] = { {0, "Rainbow"}, {1, "Blue-White-Red"}, {2, "Heat map"}, {3, "Discrete"}, {4, "NonLinear - Extremes"} };
		// Create a type for the enum ColorEV
		TwType ColorType = TwDefineEnum("ColoType", ColorEV, NUM_COLORS);

		// add 'whichColor' to 'bar': this is a variable of type ShapeType. Its key shortcuts are [+] and [-].
		TwAddVarRW(bar, "Object colors", ColorType, &whichColor, " help='Change object color.' ");
	}

	// Add modifier for the white threshold
	TwAddVarRW(bar, "WhiteThreshold", TW_TYPE_FLOAT, &g_WhiteThreshold,
		" label = 'Adjust white threshold' min=0 max=1 step=0.01 keyIncr = 'w' keyDecr = 's' help='Increase/decrease white threshold' ");

	TwAddVarCB(bar, "toggleArrows", TW_TYPE_BOOL32, setArrowCB, getArrowCB, NULL, "label='Toggle Arrows'");
	TwAddVarCB(bar, "toggleStreamlines", TW_TYPE_BOOL32, setStreamlinesCB, getStreamlinesCB, NULL, "label='Toggle Streamlines'");
	TwAddVarCB(bar, "toggleStreamribbon", TW_TYPE_BOOL32, setStreamribbonCB, getStreamribbonCB, NULL, "label='Toggle Streamribbon'");
	TwAddVarRW(bar, "modifyStreamLength", TW_TYPE_INT32, &g_streamLength, "label='Change streamline length' min=2 max=500 step=1 help='NOTE: Large values will take longer to compute'");
	TwAddVarRW(bar, "moveProbeX", TW_TYPE_FLOAT, &g_probeX, "label='Change X coordinate' min=-1.0 max=1.0 step=0.01");
	TwAddVarRW(bar, "moveProbeY", TW_TYPE_FLOAT, &g_probeY, "label='Change Y coordinate' min=-1.0 max=1.0 step=0.01");
	TwAddVarRW(bar, "moveProbeZ", TW_TYPE_FLOAT, &g_probeZ, "label='Change Z coordinate' min=-1.0 max=1.0 step=0.01");
	TwEnumVal Integrators[2] = { {0, "Euler"}, {1, "Runge-Kutta Second Order"} };
	TwType Integrator = TwDefineEnum("Integrators", Integrators, 2);
	TwAddVarRW(bar, "changeIntegrator", Integrator, &whichIntegrator, "label='Change integration method'");
	TwAddVarRW(bar, "changeDt", TW_TYPE_FLOAT, &g_step, "label='Change step size' min=0.01 max=0.5 step=0.01");
	TwAddVarRW(bar, "changeSepDist", TW_TYPE_FLOAT, &g_linedist, "label='Change dist. between lines' min=0.1 max=0.5 step=0.01");
	TwAddButton(bar, "updateProbe", recompStreamline, NULL, "label='Update Streamline/Streamribbon'");
}

// Main
int main(int argc, char *argv[])
{
	float axis[] = { 0.7f, 0.7f, 0.0f }; // initial model rotation
	float angle = 0.8f;

	// Initialize GLUT
	glutInit(&argc, argv);
	// First parameter is the buffer - single/double
	// probably no noticeable difference between single and double
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(1024, 768);
	MainWindow = glutCreateWindow("Assignment 7 � Binoy Dalal");
	glutCreateMenu(NULL);

	// Set GLUT callbacks
	glutDisplayFunc(Display); // DisplayOG is a function pointer to the original display function - function contains info. about objects to be drawn on the screen
	glutReshapeFunc(Reshape); // Reshape function tells how resizing the window will affect the graphics on screen
	atexit(Terminate);  // Called after glutMainLoop ends

	// Initialize AntTweakBar
	TwInit(TW_OPENGL, NULL);

	// Set GLUT event callbacks
	// - Directly redirect GLUT mouse button events to AntTweakBar
	glutMouseFunc((GLUTmousebuttonfun)TwEventMouseButtonGLUT);
	// - Directly redirect GLUT mouse motion events to AntTweakBar
	glutMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
	// - Directly redirect GLUT mouse "passive" motion events to AntTweakBar (same as MouseMotion)
	glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
	// - Directly redirect GLUT key events to AntTweakBar
	glutKeyboardFunc((GLUTkeyboardfun)TwEventKeyboardGLUT);
	// - Directly redirect GLUT special key events to AntTweakBar
	glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);
	// - Send 'glutGetModifers' function pointer to AntTweakBar;
	//   required because the GLUT key event functions do not report key modifiers states.
	TwGLUTModifiersFunc(glutGetModifiers);

	gen_vec_fields();

	calcLimits3dVec();
	setColorFunction();
	max_ptr = &max_sv1;
	min_ptr = &min_sv1;
	currentField = 1;
	setVecFieldPointer();
	computeStreamBunch();
	computeStreamribbon();

	// Build a display list for the axes

	InitAxesLists();

	// Initialize the AntTweakBar interface

	InitTwBar();

	// Store time
	g_RotateTime = GetTimeMs();
	// Init rotation
	SetQuaternionFromAxisAngle(axis, angle, g_Rotation);
	SetQuaternionFromAxisAngle(axis, angle, g_RotateStart);

	// Call the GLUT main loop
	glutMainLoop(); // Constantly renders the object on screen

	return 0;
}