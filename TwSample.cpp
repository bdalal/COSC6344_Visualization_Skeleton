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
// First 3 params are x, y and z - the rotation axes
// The 4th param is w - the rotation angle in radians
// For glRotate, the rotation angle will be the first param and the axes in order will be the rest
// Yrot is when y = 1 and x = z = 0
// Xrot is when x = 1 and y = z = 0
// Xrot is when x = 1 and y = z = 0
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
double abs_s_min, abs_s_max;
// pointers to max and min values in any dataset
float *max_ptr, *min_ptr;
// White threshold
float g_WhiteThreshold = 0.5;
//iso scalar value
float g_sprime = 50;
//iso contour count
int g_ncontours = 1;
// iso surface flag
bool g_isoSurfaces = false;
// bilinear flag
bool g_bilinear = false;
// opacity value
float g_opacity = 1.0;
// enable/disable slices
bool g_enableSlices = false;
// enable/disable DVR
bool g_enableDVR = true;
// enable/disable original display function
bool g_DisplayOG = false;
// enable/disable enhanced LIC
bool g_enhanceLIC = false;
// enable/disable colorplot
bool g_colorPlot = false;
// enable/disable colored LIC
bool g_coloredLIC = false;

TwBar *bar = NULL; // Pointer to the tweak bar

unsigned short isPoly = 0; // check if we're drawing quads or triangles or 3d

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

int g_XYplane = 1; // Toggle XY cutting plane
int g_YZplane = 1; // Toggle YZ cutting plane
int g_XZplane = 1; // Toggle XZ cutting plane

float g_gradientMin = 0.; // Value of change in the 3dVis
float g_gradientMax = 0.; // Value of change in the 3dVis
float g_gradientAbsMax = 0; // Set maximum value for gradient for use in interface

int g_arrows = 0; // Toggle Arrow heads

int g_streamlines = 0; // Toggle streamlines

int g_streamribbon = 0; // Toggle streamribbon

float g_probeX = 0.;
float g_probeY = 0.;
float g_probeZ = 0.;

float g_step = 0.1; //  step size when integrating

float g_linedist = 0.1; // line dist for computing stream bunch

float g_rdist = 0.1; // width of streamribbon

int whichIntegrator = 0;

#include "Skeleton.h"
Polyhedron *poly = NULL;

// To read from dat files
typedef struct node
{
	float x, y, z, s; // represents a vertex with the scalar
};
typedef struct lineseg
{
	node n1, n2; // indices of the vertices in grid_pts
	node intersection;
};
typedef struct quad
{
	node v0, v1, v2, v3; // indices of 4 vertices in grid_pts
	lineseg e0, e1, e2, e3; // indices of 4 edges in edgeList
};
typedef struct isoSurfaceNode {
	double x, y, z; // vertex
	double T; // Temperature
	float rgb[3]; // assigned color
	float rad; // radius
	float dTdx, dTdy, dTdz; // gradient in each direction
	float grad; // total gradient
	bool draw; // flag to determine if the particular node needs to be drawn or not
};
typedef struct sources {
	double xc, yc, zc;
	double a; // temperature value of the source
};
sources Sources[] = {
	{1.f, 0.f, 0.f, 90.f},
	{-1.f, 0.3f, 0.f, 120.f},
	{0.f, 1.f, 0.f, 120.f},
	{0.f, 0.4f, 1.f, 170.f}
};
typedef struct vecNode {
	float x, y, z;
	float vx, vy, vz;
	float magnitude;
};

int NX, NY;
std::vector<lineseg> isocontours; // stores all contours for quads
std::vector<std::vector<node>> grid;
std::vector<std::vector<lineseg>> rightlines;
std::vector<std::vector<lineseg>> toplines;
std::vector<std::vector<quad>> faces;
std::vector<lineseg> isocontours_t; // stores all contours for triangles

std::vector<lineseg> isosurfacecontours; // stores all contours for quads to be used in drawing the iso surface

// TODO: let user configure the size of the vector field
const int NX3d = 8, NY3d = 8, NZ3d = 8;
const double TEMPMAX = 100., TEMPMIN = 0.;
isoSurfaceNode grid3d[NX3d][NY3d][NZ3d];

vecNode vec_field1[NX3d][NY3d][NZ3d];
vecNode vec_field2[NX3d][NY3d][NZ3d];
vecNode vec_field3[NX3d][NY3d][NZ3d];
int currentField;

int g_Xslice = NX3d / 2; // default YZ plane at origin
int g_Yslice = NY3d / 2; // default XZ plane at origin
int g_Zslice = NZ3d / 2; // default XY plane at origin

unsigned char TextureXY[NZ3d][NY3d][NX3d][4];
unsigned char TextureXZ[NY3d][NZ3d][NX3d][4];
unsigned char TextureYZ[NX3d][NZ3d][NY3d][4];

//streamline length
int g_streamLength = 200;
// texture images
const int IMG_RES = 512; // resolution of the image
unsigned char noise_tex[IMG_RES][IMG_RES][3];
unsigned char vec_img[IMG_RES][IMG_RES][3];
unsigned char lic_tex[IMG_RES][IMG_RES][3];
unsigned char lic_tex_enhanced[IMG_RES][IMG_RES][3];

typedef struct streampoint {
	float nextX, nextY, nextZ;
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

int Major; /* X, Y, or Z */
int Xside, Yside, Zside; /* which side is visible, PLUS or MINUS */

const float A = sqrt(3);
const float B = sqrt(2);

void updateDataRange(void* clientData);
void setupTwBar();
void TW_CALL setXYCB(const void* value, void* clientData);
void TW_CALL getXYCB(void* value, void* clientData);
void TW_CALL setYZCB(const void* value, void* clientData);
void TW_CALL getYZCB(void* value, void* clientData);
void TW_CALL setXZCB(const void* value, void* clientData);
void TW_CALL getXZCB(void* value, void* clientData);
void Display(void);
void Reshape(int width, int height);
void ReshapeNew(int width, int height);
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
	glPushMatrix();
	glTranslatef(head[0], head[1], head[2]);
	// glScalef(0.05, 0.05, 0.05);
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
	glPopMatrix();
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

//void gen_noise_tex() {
//	for (int x = 0; x < IMG_RES; x++) {
//		for (int y = 0; y < IMG_RES; y++) {
//			float noise = 255 * (rand() % 32768) / 32768.0;
//			noise_tex[x][y][0] = noise_tex[x][y][1] = noise_tex[x][y][2] = (unsigned char)noise;
//		}
//	}
//}
//
//void determineVisibility(float mat[16]) {
//	float nzx, nzy, nzz;
//	nzx = mat[2];
//	nzy = mat[6];
//	nzz = mat[10];
//	/* which sides of the cube are showing:
//	*/
//	/* the Xside being shown to the user is MINUS or PLUS */
//	Xside = (nzx > 0. ? PLUS : MINUS);
//	Yside = (nzy > 0. ? PLUS : MINUS);
//	Zside = (nzz > 0. ? PLUS : MINUS);
//	/* which direction needs to be composited: */
//	if (fabs(nzx) > fabs(nzy) && fabs(nzx) > fabs(nzz))
//		Major = X;
//	else if (fabs(nzy) > fabs(nzx) && fabs(nzy) > fabs(nzz))
//		Major = Y;
//	else
//		Major = Z;
//}
//
//void CompositeXY(void)
//{
//	int x, y, z, zz;
//	float alpha; /* opacity at this voxel */
//	float r, g, b; /* running color composite */
//	for (x = 0; x < NX3d; x++) {
//		for (y = 0; y < NY3d; y++) {
//			r = g = b = 0.;
//			for (zz = 0; zz < NZ3d; zz++) {
//				/* which direction to composite: */
//				if (Zside == PLUS)
//					z = zz;
//				else
//					z = (NZ3d - 1) - zz;
//				isoSurfaceNode* node;
//				node = &grid3d[x][y][z];
//				if ((node->T < s_min || node->T > s_max) || (node->grad < g_gradientMin || node->grad > g_gradientMax)) { // determine whether the value is out of the range set by the range slider
//					r = g = b = 0.;
//					alpha = 0.;
//				}
//				else {
//					r = node->rgb[0];
//					g = node->rgb[1];
//					b = node->rgb[2];
//					alpha = g_opacity;
//				}
//				unsigned char ru = (unsigned char)(255.*r + .5);
//				TextureXY[zz][y][x][0] = (unsigned char)(255.*r + .5);
//				TextureXY[zz][y][x][1] = (unsigned char)(255.*g + .5);
//				TextureXY[zz][y][x][2] = (unsigned char)(255.*b + .5);
//				TextureXY[zz][y][x][3] = (unsigned char)(255.*alpha + .5);
//			}
//		}
//	}
//}
//
//void CompositeYZ(void) {
//	int x, y, z, xx;
//	float alpha; /* opacity at this voxel */
//	float r, g, b; /* running color composite */
//	for (y = 0; y < NY3d; y++) {
//		for (z = 0; z < NZ3d; z++) {
//			r = g = b = 0.;
//			for (xx = 0; xx < NX3d; xx++) {
//				/* which direction to composite: */
//				if (Xside == PLUS)
//					x = xx;
//				else
//					x = (NX3d - 1) - xx;
//				isoSurfaceNode* node;
//				node = &grid3d[x][y][z];
//				if ((node->T < s_min || node->T > s_max) || (node->grad < g_gradientMin || node->grad > g_gradientMax)) { // determine whether the value is out of the range set by the range slider
//					r = g = b = 0.;
//					alpha = 0.;
//				}
//				else {
//					r = node->rgb[0];
//					g = node->rgb[1];
//					b = node->rgb[2];
//					alpha = g_opacity;
//				}
//				TextureYZ[xx][z][y][0] = (unsigned char)(255.*r + .5);
//				TextureYZ[xx][z][y][1] = (unsigned char)(255.*g + .5);
//				TextureYZ[xx][z][y][2] = (unsigned char)(255.*b + .5);
//				TextureYZ[xx][z][y][3] = (unsigned char)(255.*alpha + .5);
//			}
//		}
//	}
//}
//
//void CompositeXZ(void) {
//	int x, y, z, yy;
//	float alpha; /* opacity at this voxel */
//	float r, g, b; /* running color composite */
//	for (x = 0; x < NX3d; x++) {
//		for (z = 0; z < NZ3d; z++) {
//			r = g = b = 0.;
//			for (yy = 0; yy < NY3d; yy++) {
//				/* which direction to composite: */
//				if (Yside == PLUS)
//					y = yy;
//				else
//					y = (NY3d - 1) - yy;
//				isoSurfaceNode* node;
//				node = &grid3d[x][y][z];
//				if ((node->T < s_min || node->T > s_max) || (node->grad < g_gradientMin || node->grad > g_gradientMax)) { // determine whether the value is out of the range set by the range slider
//					r = g = b = 0.;
//					alpha = 0.;
//				}
//				else {
//					r = node->rgb[0];
//					g = node->rgb[1];
//					b = node->rgb[2];
//					alpha = g_opacity;
//				}
//				TextureXZ[yy][z][x][0] = (unsigned char)(255.*r + .5);
//				TextureXZ[yy][z][x][1] = (unsigned char)(255.*g + .5);
//				TextureXZ[yy][z][x][2] = (unsigned char)(255.*b + .5);
//				TextureXZ[yy][z][x][3] = (unsigned char)(255.*alpha + .5);
//			}
//		}
//	}
//}
//
//void drawTexture() {
//	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
//	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
//	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
//	GLfloat filter = GL_NEAREST;
//	if (g_bilinear)
//		filter = GL_LINEAR;
//	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filter);
//	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filter);
//	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
//	glEnable(GL_TEXTURE_2D);
//	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//	glEnable(GL_BLEND);
//	if (Major == Z) {
//		float z0, dz, zcoord;
//		int z;
//		if (Zside == PLUS) {
//			z0 = -1.;
//			dz = 2. / (float)(NZ3d - 1);
//		}
//		else {
//			z0 = 1.;
//			dz = -2. / (float)(NZ3d - 1);
//		}
//		for (z = 0, zcoord = z0; z < NZ3d; z++, zcoord += dz) {
//			glTexImage2D(GL_TEXTURE_2D, 0, 4, NX3d, NY3d, 0, GL_RGBA, GL_UNSIGNED_BYTE, &TextureXY[z][0][0][0]);
//			glBegin(GL_QUADS);
//			glTexCoord2f(0., 0.);
//			glVertex3f(-1., -1., zcoord);
//			glTexCoord2f(1., 0.);
//			glVertex3f(1., -1., zcoord);
//			glTexCoord2f(1., 1.);
//			glVertex3f(1., 1., zcoord);
//			glTexCoord2f(0., 1.);
//			glVertex3f(-1., 1., zcoord);
//			glEnd();
//		}
//	}
//	else if (Major == Y) {
//		float y0, dy, ycoord;
//		int y;
//		if (Yside == PLUS) {
//			y0 = -1.;
//			dy = 2. / (float)(NY3d - 1);
//		}
//		else {
//			y0 = 1.;
//			dy = -2. / (float)(NY3d - 1);
//		}
//		for (y = 0, ycoord = y0; y < NY3d; y++, ycoord += dy) {
//			glTexImage2D(GL_TEXTURE_2D, 0, 4, NX3d, NZ3d, 0, GL_RGBA, GL_UNSIGNED_BYTE, &TextureXZ[y][0][0][0]);
//			glBegin(GL_QUADS);
//			glTexCoord2f(0., 0.);
//			glVertex3f(-1., ycoord, -1.);
//			glTexCoord2f(1., 0.);
//			glVertex3f(1., ycoord, -1.);
//			glTexCoord2f(1., 1.);
//			glVertex3f(1., ycoord, 1.);
//			glTexCoord2f(0., 1.);
//			glVertex3f(-1., ycoord, 1.);
//			glEnd();
//		}
//	}
//	else {
//		float x0, dx, xcoord;
//		int x;
//		if (Xside == PLUS) {
//			x0 = -1.;
//			dx = 2. / (float)(NX3d - 1);
//		}
//		else {
//			x0 = 1.;
//			dx = -2. / (float)(NX3d - 1);
//		}
//		for (x = 0, xcoord = x0; x < NX3d; x++, xcoord += dx) {
//			glTexImage2D(GL_TEXTURE_2D, 0, 4, NY3d, NZ3d, 0, GL_RGBA, GL_UNSIGNED_BYTE, &TextureYZ[x][0][0][0]);
//			glBegin(GL_QUADS);
//			glTexCoord2f(0., 0.);
//			glVertex3f(xcoord, -1., -1.);
//			glTexCoord2f(1., 0.);
//			glVertex3f(xcoord, 1., -1.);
//			glTexCoord2f(1., 1.);
//			glVertex3f(xcoord, 1., 1.);
//			glTexCoord2f(0., 1.);
//			glVertex3f(xcoord, -1., 1.);
//			glEnd();
//		}
//	}
//	glDisable(GL_TEXTURE_2D);
//}

//double getTemperature(double x, double y, double z) {
//	double t = 0.;
//	for (int i = 0; i < 4; i++)
//	{
//		double dx = x - Sources[i].xc;
//		double dy = y - Sources[i].yc;
//		double dz = z - Sources[i].zc;
//		double rsqd = dx * dx + dy * dy + dz * dz;
//		t += Sources[i].a * exp(-5.*rsqd);
//	}
//
//	if (t > TEMPMAX)
//		t = TEMPMAX;
//	return t;
//}
//
//void computeGradient() {
//	for (int i = 0; i < NX3d; i++) {
//		for (int j = 0; j < NY3d; j++) {
//			for (int k = 0; k < NZ3d; k++) {
//				if (i == 0)
//					grid3d[i][j][k].dTdx = (grid3d[i + 1][j][k].T - grid3d[i][j][k].T) / (grid3d[i + 1][j][k].x - grid3d[i][j][k].x);
//				else if (i == NX3d - 1)
//					grid3d[i][j][k].dTdx = (grid3d[i][j][k].T - grid3d[i - 1][j][k].T) / (grid3d[i][j][k].x - grid3d[i - 1][j][k].x);
//				else
//					grid3d[i][j][k].dTdx = (grid3d[i + 1][j][k].T - grid3d[i - 1][j][k].T) / (grid3d[i + 1][j][k].x - grid3d[i - 1][j][k].x);
//				if (j == 0)
//					grid3d[i][j][k].dTdy = (grid3d[i][j + 1][k].T - grid3d[i][j][k].T) / (grid3d[i][j + 1][k].y - grid3d[i][j][k].y);
//				else if (j == NY3d - 1)
//					grid3d[i][j][k].dTdy = (grid3d[i][j][k].T - grid3d[i][j - 1][k].T) / (grid3d[i][j][k].y - grid3d[i][j - 1][k].y);
//				else
//					grid3d[i][j][k].dTdy = (grid3d[i][j + 1][k].T - grid3d[i][j - 1][k].T) / (grid3d[i][j + 1][k].y - grid3d[i][j - 1][k].y);
//				if (k == 0)
//					grid3d[i][j][k].dTdz = (grid3d[i][j][k + 1].T - grid3d[i][j][k].T) / (grid3d[i][j][k + 1].z - grid3d[i][j][k].z);
//				else if (k == NZ3d - 1)
//					grid3d[i][j][k].dTdz = (grid3d[i][j][k].T - grid3d[i][j][k - 1].T) / (grid3d[i][j][k].z - grid3d[i][j][k - 1].z);
//				else
//					grid3d[i][j][k].dTdz = (grid3d[i][j][k + 1].T - grid3d[i][j][k - 1].T) / (grid3d[i][j][k + 1].z - grid3d[i][j][k - 1].z);
//
//				grid3d[i][j][k].grad = sqrt(pow(grid3d[i][j][k].dTdx, 2) + pow(grid3d[i][j][k].dTdy, 2) + pow(grid3d[i][j][k].dTdz, 2));
//			}
//		}
//	}
//}
//
//void populate3dStruct() {
//	double ix = 2. / NX3d; // splitting data into 50 samples along the x axis
//	double iy = 2. / NY3d;
//	double iz = 2. / NZ3d;
//	for (int i = 0; i < NX3d; i++) {
//		for (int j = 0; j < NY3d; j++) {
//			for (int k = 0; k < NZ3d; k++) {
//				isoSurfaceNode node;
//				node.x = (ix * i) - 1;
//				node.y = (iy * j) - 1;
//				node.z = (iz * k) - 1;
//				node.T = getTemperature(node.x, node.y, node.z);
//				node.draw = true;
//				grid3d[i][j][k] = node;
//			}
//		}
//	}
//}
//
//void color3dStruct() {
//	for (int i = 0; i < NX3d; i++) {
//		for (int j = 0; j < NY3d; j++) {
//			for (int k = 0; k < NZ3d; k++) {
//				isoSurfaceNode* node;
//				node = &grid3d[i][j][k];
//				if ((node->T >= s_min && node->T <= s_max) && (node->grad >= g_gradientMin && node->grad <= g_gradientMax)) {
//					colorFunction(node->T, node->rgb);
//					node->draw = true;
//				}
//				else
//					node->draw = false;
//			}
//		}
//	}
//}

//void Load_data_on_uniformGrids(const char *name)
//{
//	int i, j;
//	FILE* fp = fopen(name, "r");
//	if (fp == NULL) return;
//	fscanf(fp, "%d %d\n", &NX, &NY);
//	grid.clear();
//	for (i = 0; i < NY; i++) {
//		std::vector<node> tmpv;
//		for (j = 0; j < NX; j++) {
//			node tmp;
//			fscanf(fp, "%f, %f, %f, %f \n", &tmp.x, &tmp.y, &tmp.z, &tmp.s);
//			tmpv.push_back(tmp);
//		}
//		grid.push_back(tmpv);
//	}
//	fclose(fp);
//}
//
//void build_edge_list() {
//	int i, j;
//	rightlines.clear();
//	for (i = 0; i < NY; i++) {
//		std::vector<lineseg> tmpl;
//		for (j = 0; j < NX - 1; j++) {
//			node n1 = grid[i][j];
//			node n2 = grid[i][j + 1];
//			lineseg rightedge;
//			rightedge.n1 = n1;
//			rightedge.n2 = n2;
//			rightedge.intersection;
//			tmpl.push_back(rightedge);
//		}
//		rightlines.push_back(tmpl);
//	}
//	toplines.clear();
//	for (i = 0; i < NY - 1; i++) {
//		std::vector<lineseg> tmpl;
//		for (j = 0; j < NX; j++) {
//			node n1 = grid[i][j];
//			node n2 = grid[i + 1][j];
//			lineseg topedge;
//			topedge.n1 = n1;
//			topedge.n2 = n2;
//			topedge.intersection;
//			tmpl.push_back(topedge);
//		}
//		toplines.push_back(tmpl);
//	}
//}
//
//void build_face_list() {
//	int i, j;
//	faces.clear();
//	for (i = 0; i < NY - 1; i++) {
//		std::vector<quad> tmpf;
//		for (j = 0; j < NX - 1; j++) {
//			quad face;
//			node v0 = grid[i][j];
//			node v1 = grid[i][j + 1];
//			node v2 = grid[i + 1][j + 1];
//			node v3 = grid[i + 1][j];
//			lineseg ebr = rightlines[i][j]; // bottom right edge
//			lineseg etr = rightlines[i + 1][j]; // top right edge
//			lineseg elt = toplines[i][j]; // left top edge
//			lineseg ert = toplines[i][j + 1]; // right top edge
//			face.v0 = v0;
//			face.v1 = v1;
//			face.v2 = v2;
//			face.v3 = v3;
//			face.e0 = ebr;
//			face.e1 = ert;
//			face.e2 = etr;
//			face.e3 = elt;
//			tmpf.push_back(face);
//		}
//		faces.push_back(tmpf);
//	}
//}
//
//void copyVertices(isoSurfaceNode* s1, isoSurfaceNode* s2, node* n1, node* n2) {
//	n1->x = s1->x;
//	n1->y = s1->y;
//	n1->z = s1->z;
//	n2->x = s2->x;
//	n2->y = s2->y;
//	n2->z = s2->z;
//}
//
//void computeIsoSurfaces(isoSurfaceNode* n[]) {
//	bool intersects[4] = { false, false, false, false };
//	node intersections[4];
//	int ctr = 0;
//	for (int l = 0; l < 4; l++) {
//		isoSurfaceNode *s1, *s2;
//		switch (l) {
//		case 0:
//			s1 = n[0];
//			s2 = n[1];
//			break;
//		case 1:
//			s1 = n[1];
//			s2 = n[2];
//			break;
//		case 2:
//			s1 = n[2];
//			s2 = n[3];
//			break;
//		case 3:
//			s1 = n[3];
//			s2 = n[0];
//			break;
//		}
//		if (s1->T == s2->T && s1->T != g_sprime)
//			continue;
//		if (s1->T == s2->T && s1->T == g_sprime) {
//			node node1, node2;
//			copyVertices(s1, s2, &node1, &node2);
//			lineseg contour;
//			contour.n1 = node1;
//			contour.n2 = node2;
//			isosurfacecontours.push_back(contour);
//			ctr++;
//			continue;
//		}
//		float tprime;
//		tprime = (g_sprime - s1->T) / (s2->T - s1->T);
//		if (tprime >= 0 && tprime <= 1) {
//			node intersect;
//			ctr++;
//			intersect.x = ((1 - tprime) * s1->x) + (tprime * s2->x);
//			intersect.y = ((1 - tprime) * s1->y) + (tprime * s2->y);
//			intersect.z = ((1 - tprime) * s1->z) + (tprime * s2->z);
//			intersections[l] = intersect;
//			intersects[l] = true;
//		}
//	}
//	if (ctr == 0)
//		return;
//	if (ctr == 1 || ctr == 3)
//		char *c = "Something is very wrong";
//	if (ctr == 2) {
//		lineseg contour;
//		if (intersects[0]) {
//			contour.n1 = intersections[0];
//			if (intersects[1])
//				contour.n2 = intersections[1];
//			else if (intersects[2])
//				contour.n2 = intersections[2];
//			else
//				contour.n2 = intersections[3];
//		}
//		else if (intersects[1]) {
//			contour.n1 = intersections[1];
//			if (intersects[2])
//				contour.n2 = intersections[2];
//			else
//				contour.n2 = intersections[3];
//		}
//		else {
//			contour.n1 = intersections[2];
//			contour.n2 = intersections[3];
//		}
//		isosurfacecontours.push_back(contour);
//		return;
//	}
//	if (ctr == 4) {
//		float s0 = n[0]->T;
//		float s1 = n[1]->T;
//		float s2 = n[2]->T;
//		float s3 = n[3]->T;
//		float m = (s0 + s1 + s2 + s3) / 4;
//		lineseg contour1, contour2;
//		if (g_sprime <= m) {
//			if (s0 <= m) {
//				contour1.n1 = intersections[0];
//				contour1.n2 = intersections[3];
//				contour2.n1 = intersections[1];
//				contour2.n2 = intersections[2];
//			}
//			else {
//				contour1.n1 = intersections[0];
//				contour1.n2 = intersections[1];
//				contour2.n1 = intersections[3];
//				contour2.n2 = intersections[2];
//			}
//		}
//		else {
//			if (s0 <= m) {
//				contour1.n1 = intersections[0];
//				contour1.n2 = intersections[1];
//				contour2.n1 = intersections[3];
//				contour2.n2 = intersections[2];
//			}
//			else {
//				contour1.n1 = intersections[0];
//				contour1.n2 = intersections[3];
//				contour2.n1 = intersections[1];
//				contour2.n2 = intersections[2];
//			}
//		}
//		isosurfacecontours.push_back(contour1);
//		isosurfacecontours.push_back(contour2);
//	}
//}
//
//void computeIsoSurfaces(float g_sprime) {
//	if (g_sprime < s_min || g_sprime > s_max)
//		return;
//	isosurfacecontours.clear();
//	int i, j, k;
//	for (k = 0; k < NZ3d; k++) {
//		for (i = 0; i < NX3d - 1; i++) {
//			for (j = 0; j < NY3d - 1; j++) {
//				// Process quad whose corner is at [i, j, k] in XY plane
//				isoSurfaceNode* n[4] = { &grid3d[i][j][k], &grid3d[i + 1][j][k], &grid3d[i + 1][j + 1][k], &grid3d[i][j + 1][k] };
//				computeIsoSurfaces(n);
//			}
//		}
//	}
//	for (i = 0; i < NX3d; i++) {
//		for (k = 0; k < NZ3d - 1; k++) {
//			for (j = 0; j < NY3d - 1; j++) {
//				// Process quad whose corner is at [i, j, k] in YZ plane
//				isoSurfaceNode* n[4] = { &grid3d[i][j][k], &grid3d[i][j + 1][k], &grid3d[i][j + 1][k + 1], &grid3d[i][j][k + 1] };
//				computeIsoSurfaces(n);
//			}
//		}
//	}
//	for (j = 0; j < NY3d; j++) {
//		for (k = 0; k < NZ3d - 1; k++) {
//			for (i = 0; i < NX3d - 1; i++) {
//				// Process quad whose corner is at [i, j, k] in XZ plane
//				isoSurfaceNode* n[4] = { &grid3d[i][j][k], &grid3d[i + 1][j][k], &grid3d[i + 1][j][k + 1], &grid3d[i][j][k + 1] };
//				computeIsoSurfaces(n);
//			}
//		}
//	}
//}
//
//void computeContours(float g_sprime) {
//	if (g_sprime < s_min || g_sprime > s_max)
//		return;
//	// for each face compute the intersection
//	int i, j, k;
//	for (i = 0; i < faces.size(); i++) {
//		std::vector<quad> tmpf = faces[i];
//		for (j = 0; j < tmpf.size(); j++) { // get edges for every face in a row
//			quad face = tmpf[j];
//			lineseg* e[4] = { &face.e0, &face.e1, &face.e2, &face.e3 };
//			bool intersects[4] = { false, false, false, false };
//			int ctr = 0;
//			for (k = 0; k < 4; k++) {
//				// compute intersection for every edge and store on edge itself
//				node n1 = e[k]->n1;
//				node n2 = e[k]->n2;
//				if (n1.s == n2.s && n1.s != g_sprime)
//					continue;
//				if (n1.s == n2.s && n1.s == g_sprime) {
//					lineseg contour;
//					contour.n1 = n1;
//					contour.n2 = n2;
//					isocontours.push_back(contour);
//					ctr++; // there's an intersection
//					continue;
//				}
//				float tprime, xprime, yprime, zprime;
//				tprime = (g_sprime - n1.s) / (n2.s - n1.s);
//				if (tprime >= 0 && tprime <= 1) {
//					ctr++; // there's an intersection so increment the counter
//					xprime = ((1 - tprime) * n1.x) + (tprime * n2.x);
//					yprime = ((1 - tprime) * n1.y) + (tprime * n2.y);
//					zprime = ((1 - tprime) * n1.z) + (tprime * n2.z);
//					e[k]->intersection.x = xprime;
//					e[k]->intersection.y = yprime;
//					e[k]->intersection.z = zprime;
//					e[k]->intersection.s = g_sprime;
//					intersects[k] = true;
//				}
//			}
//			if (ctr == 0)
//				continue;
//			if (ctr == 1 || ctr == 3) {
//				char* c = "Something is very wrong";
//			}
//			// connect the intersections based on the count
//			if (ctr == 2) { // if there are 2 intersections
//				lineseg contour;
//				if (intersects[0]) {
//					contour.n1 = e[0]->intersection;
//					if (intersects[1])
//						contour.n2 = e[1]->intersection;
//					else if (intersects[2])
//						contour.n2 = e[2]->intersection;
//					else
//						contour.n2 = e[3]->intersection;
//				}
//				else if (intersects[1]) {
//					contour.n1 = e[1]->intersection;
//					if (intersects[2])
//						contour.n2 = e[2]->intersection;
//					else
//						contour.n2 = e[3]->intersection;
//				}
//				else {
//					contour.n1 = e[2]->intersection;
//					contour.n2 = e[3]->intersection;
//				}
//				isocontours.push_back(contour);
//				continue;
//			}
//			if (ctr == 4) { // If there are 4 intersections. Condition hit for dataset2 at g_sprime = 49.000
//				float s0 = face.v0.s;
//				float s1 = face.v1.s;
//				float s2 = face.v2.s;
//				float s3 = face.v3.s;
//				float m = (s0 + s1 + s2 + s3) / 4;
//				lineseg contour1, contour2;
//				if (g_sprime <= m) {
//					if (s0 <= m) {
//						contour1.n1 = e[0]->intersection;
//						contour1.n2 = e[3]->intersection;
//						contour2.n1 = e[1]->intersection;
//						contour2.n2 = e[2]->intersection;
//					}
//					else {
//						contour1.n1 = e[0]->intersection;
//						contour1.n2 = e[1]->intersection;
//						contour2.n1 = e[3]->intersection;
//						contour2.n2 = e[2]->intersection;
//					}
//				}
//				else {
//					if (s0 <= m) {
//						contour1.n1 = e[0]->intersection;
//						contour1.n2 = e[1]->intersection;
//						contour2.n1 = e[3]->intersection;
//						contour2.n2 = e[2]->intersection;
//					}
//					else {
//						contour1.n1 = e[0]->intersection;
//						contour1.n2 = e[3]->intersection;
//						contour2.n1 = e[1]->intersection;
//						contour2.n2 = e[2]->intersection;
//					}
//				}
//				isocontours.push_back(contour1);
//				isocontours.push_back(contour2);
//			}
//		}
//	}
//}
//
//void computeContoursTriangles(float g_sprime) {
//	if (g_sprime < s_min || g_sprime > s_max)
//		return;
//	for (int i = 0; i < poly->ntris; i++) {
//		Triangle* tmpt = poly->tlist[i];
//		Edge **e = tmpt->edges;
//		int ctr = 0;
//		std::vector<node> intersections;
//		for (int j = 0; j < 3; j++) {
//			Vertex** v = e[j]->verts;
//			Vertex* v0 = v[0];
//			Vertex* v1 = v[1];
//			if (v0->s == v1->s && v0->s != g_sprime)
//				continue;
//			if (v1->s == v1->s && v0->s == g_sprime) {
//				lineseg contour;
//				node n1, n2;
//				n1.x = v0->x;
//				n2.x = v1->x;
//				n1.y = v0->y;
//				n2.y = v1->y;
//				n1.z = v0->z;
//				n2.z = v1->z;
//				n1.s = v0->s;
//				n2.s = v1->s;
//				contour.n1 = n1;
//				contour.n2 = n2;
//				ctr++;
//				isocontours_t.push_back(contour);
//				break;
//			}
//			float tprime, xprime, yprime, zprime;
//			node intersection;
//			tprime = (g_sprime - v0->s) / (v1->s - v0->s);
//			if (tprime >= 0 && tprime <= 1) {
//				xprime = ((1 - tprime) * v0->x) + (tprime * v1->x);
//				yprime = ((1 - tprime) * v0->y) + (tprime * v1->y);
//				zprime = ((1 - tprime) * v0->z) + (tprime * v1->z);
//				intersection.x = xprime;
//				intersection.y = yprime;
//				intersection.z = zprime;
//				intersection.s = g_sprime;
//				intersections.push_back(intersection);
//				ctr++;
//			}
//		}
//		if (ctr == 1 || ctr == 3) {
//			char* c = "Something is very wrong";
//		}
//		if (ctr == 2) {
//			lineseg contour;
//			node n1, n2;
//			n1.x = intersections[0].x;
//			n1.y = intersections[0].y;
//			n1.z = intersections[0].z;
//			n1.s = intersections[0].s;
//			n2.x = intersections[1].x;
//			n2.y = intersections[1].y;
//			n2.z = intersections[1].z;
//			n2.s = intersections[1].s;
//			contour.n1 = n1;
//			contour.n2 = n2;
//			isocontours_t.push_back(contour);
//		}
//	}
//}


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

//void draw_arrows(double head[2], float direct[2])
//{
//	if (g_colorPlot)
//		glColor3f(0., 0., 0.);
//	else
//		glColor3f(1., 1., 0.);
//	glPushMatrix();
//	glTranslatef(head[0], head[1], 0);
//	glRotatef(atan2(direct[1], direct[0]) * 360 / (2 * M_PI), 0, 0, 1);
//	// draw arrow head
//	glScalef(0.03, 0.03, 1);
//	glBegin(GL_TRIANGLES);
//	glVertex2f(0, 0);
//	glVertex2f(-0.35, 0.12);
//	glVertex2f(-0.35, -0.12);
//	glEnd();
//	// draw arrow body
//	glScalef(0.3, 0.3, 1);
//	//glScalef(100. / poly->nverts, 100. / poly->nverts, 1);
//	glBegin(GL_LINES);
//	glVertex2f(0, 0);
//	glVertex2f(-3, 0);
//	glEnd();
//	glPopMatrix();
//}

//void drawArrowPlot() {
//	for (int i = 0; i < poly->nverts; i++) {
//		Vertex *temp_v = poly->vlist[i];
//		double arrow_head[2] = { temp_v->x, temp_v->y };
//		float arrow_direct[2] = { temp_v->vx, temp_v->vy };
//		draw_arrows(arrow_head, arrow_direct);
//	}
//}

// TODO: user interface to scale the arrow size

void draw_3d_arrows(float head[3], float direct[3], float magnitude) {
	float rgb[3];
	colorFunction(magnitude, rgb);
	glColor3f(rgb[0], rgb[1], rgb[2]);
	glPushMatrix();
	glTranslatef(head[0], head[1], head[2]);
	glRotatef(atan2(direct[1], direct[0]) * 360 / (2 * M_PI), 0, 0, 1);
	glRotatef(atan2(direct[2], direct[1]) * 360 / (2 * M_PI), 1, 0, 0);
	glRotatef(atan2(direct[2], direct[0]) * 360 / (2 * M_PI), 0, 1, 0);
	// draw arrow head
	glScalef(0.1, 0.1, 0.1);
	glBegin(GL_TRIANGLES);
	glVertex2f(0, 0);
	glVertex2f(-0.5, 0.5);
	glVertex2f(-0.5, -0.5);
	glEnd();
	// draw arrow body
	glScalef(0.3, 0.3, 0.3);
	glBegin(GL_LINES);
	glVertex2f(0, 0);
	glVertex2f(-3, 0);
	glEnd();
	glPopMatrix();
}


void draw_3d_arrows_field1() {
	for (int i = 0; i < NX3d; i++) {
		for (int j = 0; j < NY3d; j++) {
			for (int k = 0; k < NZ3d; k++) {
				vecNode *node = &vec_field1[i][j][k];
				// TODO normalize vectors before passing to the Arrow function
				float arrow_head[3] = { node->x, node->y, node->z };
				float arrow_direct[3] = { node->vx, node->vy, node->vz };
				// draw_3d_arrows(arrow_head, arrow_direct, node->magnitude);
				Arrow(arrow_head, arrow_direct);
			}
		}
	}
}

void draw_3d_arrows_field2() {
	for (int i = 0; i < NX3d; i++) {
		for (int j = 0; j < NY3d; j++) {
			for (int k = 0; k < NZ3d; k++) {
				vecNode *node = &vec_field2[i][j][k];
				float arrow_head[3] = { node->x, node->y, node->z };
				float arrow_direct[3] = { node->vx, node->vy, node->vz };
				// draw_3d_arrows(arrow_head, arrow_direct, node->magnitude);
				Arrow(arrow_head, arrow_direct);
			}
		}
	}
}

void draw_3d_arrows_field3() {
	for (int i = 0; i < NX3d; i++) {
		for (int j = 0; j < NY3d; j++) {
			for (int k = 0; k < NZ3d; k++) {
				vecNode *node = &vec_field3[i][j][k];
				float arrow_head[3] = { node->x, node->y, node->z };
				float arrow_direct[3] = { node->vx, node->vy, node->vz };
				// draw_3d_arrows(arrow_head, arrow_direct, node->magnitude);
				Arrow(arrow_head, arrow_direct);
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

//void setExtremePointers() {
//	switch (whichPlot)
//	{
//	case 0:
//		max_ptr = &s_max;
//		min_ptr = &s_min;
//		break;
//	case 1:
//		max_ptr = &a_max;
//		min_ptr = &a_min;
//		break;
//	case 2:
//		max_ptr = &vx_max;
//		min_ptr = &vx_min;
//		break;
//	case 3:
//		max_ptr = &vy_max;
//		min_ptr = &vy_min;
//	default:
//		break;
//	}
//}

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

//void calcLimits() {
//	if (isPoly == 0) {
//		s_max = s_min = poly->tlist[0]->verts[0]->magnitude;
//		a_max = a_min = poly->tlist[0]->verts[0]->angle;
//		vx_max = vx_min = poly->tlist[0]->verts[0]->vx;
//		vy_max = vy_min = poly->tlist[0]->verts[0]->vy;
//		for (int i = 0; i < poly->ntris; i++) {
//			Triangle *temp_t = poly->tlist[i];
//			for (int j = 0; j < 3; j++) {
//				Vertex *temp_v = temp_t->verts[j];
//				float magnitude = temp_v->magnitude;
//				float angle = temp_v->angle;
//				float vx = temp_v->vx;
//				float vy = temp_v->vy;
//				if (magnitude > s_max)
//					s_max = magnitude;
//				if (magnitude < s_min)
//					s_min = magnitude;
//				if (angle > a_max)
//					a_max = angle;
//				if (angle < a_min)
//					a_min = angle;
//				if (vx > vx_max)
//					vx_max = vx;
//				if (vx < vx_min)
//					vx_min = vx;
//				if (vy > vy_max)
//					vy_max = vy;
//				if (vy < vy_min)
//					vy_min = vy;
//			}
//		}
//	}
//	else if (isPoly == 1) {
//		s_max = s_min = grid[0][0].s;
//		for (int i = 0; i < NY; i++) {
//			for (int j = 0; j < NX; j++) {
//				if (grid[i][j].s < s_min)
//					s_min = grid[i][j].s;
//				if (grid[i][j].s > s_max)
//					s_max = grid[i][j].s;
//			}
//		}
//	}
//	else {
//		s_max = s_min = grid3d[0][0][0].T;
//		g_gradientMax = g_gradientMin = grid3d[0][0][0].grad;
//		for (int i = 0; i < NX3d; i++) {
//			for (int j = 0; j < NY3d; j++) {
//				for (int k = 0; k < NZ3d; k++) {
//					if (grid3d[i][j][k].T < s_min)
//						s_min = grid3d[i][j][k].T;
//					if (grid3d[i][j][k].T > s_max)
//						s_max = grid3d[i][j][k].T;
//					if (grid3d[i][j][k].grad < g_gradientMin)
//						g_gradientMin = grid3d[i][j][k].grad;
//					if (grid3d[i][j][k].grad > g_gradientMax)
//						g_gradientMax = grid3d[i][j][k].grad;
//				}
//			}
//		}
//		abs_s_max = s_max;
//		abs_s_min = s_min;
//		g_gradientAbsMax = g_gradientMax;
//	}
//}

//void drawSquareObject() {
//	for (int i = 0; i < NY - 1; i++) {
//		for (int j = 0; j < NX - 1; j++) {
//			quad face = faces[i][j];
//			float rgb[3];
//			node v;
//			glBegin(GL_QUADS);
//			v = face.v0;
//			colorFunction(v.s, rgb);
//			glColor3f(rgb[0], rgb[1], rgb[2]);
//			glVertex3f(v.x, v.y, v.z);
//			v = face.v1;
//			colorFunction(v.s, rgb);
//			glColor3f(rgb[0], rgb[1], rgb[2]);
//			glVertex3f(v.x, v.y, v.z);
//			v = face.v2;
//			colorFunction(v.s, rgb);
//			glColor3f(rgb[0], rgb[1], rgb[2]);
//			glVertex3f(v.x, v.y, v.z);
//			v = face.v3;
//			colorFunction(v.s, rgb);
//			glColor3f(rgb[0], rgb[1], rgb[2]);
//			glVertex3f(v.x, v.y, v.z);
//			glEnd();
//		}
//	}
//	// draw the contours
//	glColor3f(0, 0, 0);
//	for (int i = 0; i < isocontours.size(); i++) {
//		node v1 = isocontours[i].n1;
//		node v2 = isocontours[i].n2;
//		glBegin(GL_LINES);
//		glVertex3f(v1.x, v1.y, v1.z);
//		glVertex3f(v2.x, v2.y, v2.z);
//		glEnd();
//	}
//}
//
//void drawTriangularObject() {
//	// draw and color object
//	for (int i = 0; i < poly->ntris; i++) {
//		Triangle *temp_t = poly->tlist[i];
//		glBegin(GL_POLYGON);
//		for (int j = 0; j < 3; j++) {
//			Vertex *temp_v = temp_t->verts[j];
//			glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
//			float rgb[3];
//			float x = temp_v->x;
//			float y = temp_v->y;
//			float z = temp_v->z;
//			// float s = temp_v->s;
//			if (whichPlot == 0)
//				colorFunction(temp_v->magnitude, rgb);
//			else if (whichPlot == 1)
//				colorFunction(temp_v->angle, rgb);
//			else if (whichPlot == 2)
//				colorFunction(temp_v->vx, rgb);
//			else if (whichPlot == 3)
//				colorFunction(temp_v->vy, rgb);
//			glColor4f(rgb[0], rgb[1], rgb[2], 1);
//			glVertex3d(x, y, z);
//		}
//		glEnd();
//	}
//	if (g_arrows)
//		drawArrowPlot();
//	// draw contours
//	/*glColor3f(0, 0, 0);
//	for (int i = 0; i < isocontours_t.size(); i++) {
//		node v1 = isocontours_t[i].n1;
//		node v2 = isocontours_t[i].n2;
//		glBegin(GL_LINES);
//		glVertex3f(v1.x, v1.y, v1.z);
//		glVertex3f(v2.x, v2.y, v2.z);
//		glEnd();
//	}*/
//}

//void draw3dObject(void(*colorFunction)(float s, float rgb[3])) {
//	for (int i = 0; i < NX3d; i++) {
//		for (int j = 0; j < NY3d; j++) {
//			for (int k = 0; k < NZ3d; k++) {
//				float rgb[3];
//				isoSurfaceNode current = grid3d[i][j][k];
//				glBegin(GL_POINTS);
//				colorFunction(current.T, rgb);
//				glColor3d(rgb[0], rgb[1], rgb[2]);
//				glVertex3d(current.x, current.y, current.z);
//				glEnd();
//			}
//		}
//	}
//}
//
//void draw3dVis() {
//	isoSurfaceNode curr;
//	if (g_enableSlices) {
//		if (g_XYplane) {
//			for (int i = 0; i < NX3d - 1; i++) {
//				for (int j = 0; j < NY3d - 1; j++) {
//					// construct a plane and color it
//					curr = grid3d[i][j][g_Zslice];
//					if (!curr.draw)
//						continue;
//					glBegin(GL_QUADS);
//					glColor3d(curr.rgb[0], curr.rgb[1], curr.rgb[2]);
//					glVertex3f(curr.x, curr.y, curr.z);
//					curr = grid3d[i + 1][j][g_Zslice];
//					glColor3d(curr.rgb[0], curr.rgb[1], curr.rgb[2]);
//					glVertex3f(curr.x, curr.y, curr.z);
//					curr = grid3d[i + 1][j + 1][g_Zslice];
//					glColor3d(curr.rgb[0], curr.rgb[1], curr.rgb[2]);
//					glVertex3f(curr.x, curr.y, curr.z);
//					curr = grid3d[i][j + 1][g_Zslice];
//					glColor3d(curr.rgb[0], curr.rgb[1], curr.rgb[2]);
//					glVertex3f(curr.x, curr.y, curr.z);
//					glEnd();
//				}
//			}
//		}
//		if (g_YZplane) {
//			for (int j = 0; j < NY3d - 1; j++) {
//				for (int k = 0; k < NZ3d - 1; k++) {
//					// construct a plane and color it
//					curr = grid3d[g_Xslice][j][k];
//					if (!curr.draw)
//						continue;
//					glBegin(GL_QUADS);
//					glColor3d(curr.rgb[0], curr.rgb[1], curr.rgb[2]);
//					glVertex3f(curr.x, curr.y, curr.z);
//					curr = grid3d[g_Xslice][j + 1][k];
//					glColor3d(curr.rgb[0], curr.rgb[1], curr.rgb[2]);
//					glVertex3f(curr.x, curr.y, curr.z);
//					curr = grid3d[g_Xslice][j + 1][k + 1];
//					glColor3d(curr.rgb[0], curr.rgb[1], curr.rgb[2]);
//					glVertex3f(curr.x, curr.y, curr.z);
//					curr = grid3d[g_Xslice][j][k + 1];
//					glColor3d(curr.rgb[0], curr.rgb[1], curr.rgb[2]);
//					glVertex3f(curr.x, curr.y, curr.z);
//					glEnd();
//				}
//			}
//		}
//		if (g_XZplane) {
//			for (int i = 0; i < NX3d - 1; i++) {
//				for (int k = 0; k < NZ3d - 1; k++) {
//					// construct a plane and color it
//					curr = grid3d[i][g_Yslice][k];
//					if (!curr.draw)
//						continue;
//					glBegin(GL_QUADS);
//					glColor3d(curr.rgb[0], curr.rgb[1], curr.rgb[2]);
//					glVertex3f(curr.x, curr.y, curr.z);
//					curr = grid3d[i + 1][g_Yslice][k];
//					glColor3d(curr.rgb[0], curr.rgb[1], curr.rgb[2]);
//					glVertex3f(curr.x, curr.y, curr.z);
//					curr = grid3d[i + 1][g_Yslice][k + 1];
//					glColor3d(curr.rgb[0], curr.rgb[1], curr.rgb[2]);
//					glVertex3f(curr.x, curr.y, curr.z);
//					curr = grid3d[i][g_Yslice][k + 1];
//					glColor3d(curr.rgb[0], curr.rgb[1], curr.rgb[2]);
//					glVertex3f(curr.x, curr.y, curr.z);
//					glEnd();
//				}
//			}
//		}
//	}
//
//	// Draw contours if flag is set
//	if (g_isoSurfaces) {
//		glColor3f(0, 1, 1);
//		for (int i = 0; i < isosurfacecontours.size(); i++) {
//			node* v1 = &isosurfacecontours[i].n1;
//			node* v2 = &isosurfacecontours[i].n2;
//			glBegin(GL_LINES);
//			glVertex3f(v1->x, v1->y, v1->z);
//			glVertex3f(v2->x, v2->y, v2->z);
//			glEnd();
//		}
//	}
//}

//void DisplayNew(void) {
//	glViewport(0, 0, (GLsizei)IMG_RES, (GLsizei)IMG_RES);
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//	gluOrtho2D(0, 1, 0, 1);
//	glEnable(GL_TEXTURE_2D);
//	glEnable(GL_BLEND);
//	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//	glClear(GL_COLOR_BUFFER_BIT);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
//	glShadeModel(GL_FLAT);
//	////Test noise texture (for debugging purpose)
//	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, IMG_RES, IMG_RES, 0,	GL_RGB, GL_UNSIGNED_BYTE, noise_tex);
//	//// Test vector field image (for debugging purpose)
//	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, IMG_RES, IMG_RES, 0,	GL_RGB, GL_UNSIGNED_BYTE, vec_img);
//	// Display LIC image using texture mapping
//	if (!g_enhanceLIC)
//		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, IMG_RES, IMG_RES, 0, GL_RGB, GL_UNSIGNED_BYTE, lic_tex);
//	else if (g_enhanceLIC)
//		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, IMG_RES, IMG_RES, 0, GL_RGB, GL_UNSIGNED_BYTE, lic_tex_enhanced);
//	glBegin(GL_QUAD_STRIP);
//	glTexCoord2f(0.0, 0.0); glVertex2f(0.0, 0.0);
//	glTexCoord2f(0.0, 1.0); glVertex2f(0.0, 1.0);
//	glTexCoord2f(1.0, 0.0); glVertex2f(1.0, 0.0);
//	glTexCoord2f(1.0, 1.0); glVertex2f(1.0, 1.0);
//	glEnd();
//	glDisable(GL_TEXTURE_2D);
//	if (g_arrows)
//		drawArrowPlot();
//	if (g_coloredLIC) {
//		glShadeModel(GL_SMOOTH);
//		setColorFunction();
//		/*setExtremePointers();*/
//		for (int i = 0; i < poly->ntris; i++) {
//			Triangle *t = poly->tlist[i];
//			glBegin(GL_POLYGON);
//			for (int j = 0; j < 3; j++) {
//				Vertex *v = t->verts[j];
//				glNormal3d(v->normal.entry[0], v->normal.entry[1], v->normal.entry[2]);
//				float rgb[3];
//				if (whichPlot == 0)
//					colorFunction(v->magnitude, rgb);
//				else if (whichPlot == 1)
//					colorFunction(v->angle, rgb);
//				else if (whichPlot == 2)
//					colorFunction(v->vx, rgb);
//				else if (whichPlot == 3)
//					colorFunction(v->vy, rgb);
//				glColor4f(rgb[0], rgb[1], rgb[2], 0.45);
//				glVertex3d(v->x, v->y, v->z);
//			}
//			glEnd();
//		}
//	}
//	if (g_colorPlot) {
//		g_coloredLIC = false;
//		glDisable(GL_BLEND);
//		glShadeModel(GL_SMOOTH);
//		setColorFunction();
//		/*setExtremePointers();*/
//		drawTriangularObject();
//	}
//	TwDraw();
//	glutSwapBuffers();
//	glFlush();
//	glutPostRedisplay();
//}


void draw_streamlines() {
	for (int i = 0; i < (int) streamlines.size(); i++) {
		std::vector<streampoint> streamline = streamlines[i];
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j < streamline.size(); j++) {
			streampoint *point = &streamline[j];
			glVertex3f(point->nextX, point->nextY, point->nextZ);
		}
		glEnd();
	}
}

// TODO: color streamline

void draw_streamribbon() {
	for (int i = 0; i < (int) streamr1.size() - 1; i++) {
		glBegin(GL_QUADS);
		streampoint *p1 = &streamr1[i];
		streampoint *p2 = &streamr2[i];
		glVertex3f(p1->nextX, p1->nextY, p1->nextZ);
		glVertex3f(p2->nextX, p2->nextY, p2->nextZ);
		++p1;
		++p2;
		glVertex3f(p2->nextX, p2->nextY, p2->nextZ);
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

	//setColorFunction();
	//setExtremePointers();
	//// Draw the 3D object
	//if (isPoly == 1) drawSquareObject();
	//else if (isPoly == 0) drawTriangularObject();
	////else draw3dObject(colorFunction);
	//else {
	//	draw3dVis();
	//	if (g_enableDVR) {
	//		determineVisibility(mat);
	//		CompositeXY();
	//		CompositeXZ();
	//		CompositeYZ();
	//		drawTexture();
	//	}
	//}

	draw_cube();

	if (g_arrows) {
		if (currentField == 1)
			draw_3d_arrows_field1();
		else if (currentField == 2)
			draw_3d_arrows_field2();
		else
			draw_3d_arrows_field3();
	}

	if (g_streamlines)
		draw_streamlines();

	if (g_streamribbon)
		draw_streamribbon();

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

void ReshapeNew(int width, int height) {
	// Set OpenGL viewport and camera
	glViewport(0, 0, (GLsizei)width, (GLsizei)height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, 1, 0, 1);
	glClearColor(0, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Send the new window size to AntTweakBar
	TwWindowSize(width, height);
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

// TODO include note in report about the benefit of RK2 but visualizing the stream ribbon for step size = 0.5 and obesrving the error for euler vs rk2
// TODO include note in report about the origin being a fixed point for field 3 which is why the vector field won't show
// TODO include note in report about streamlines in bunch being cut off because of the separation condition  - Lec10.pdf slide 48 condition 1

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

bool checkConditions(float next_x, float next_y, float next_z, int n_steps) {
	bool conditions, condition1, condition2, condition3, condition4, condition5;	
	// next value should be within bounds
	condition1 = (next_x >= -1. && next_x <= 1.) && (next_y >= -1. && next_y <= 1.) && (next_z >= -1. && next_z <= 1.);
	// next value is not a fixed point
	condition2 = getMagnitude(next_x, next_y, next_z) > 0.000001;
	// streamline isn't curving back on itself - check if the euclidean distance is greater than threshold
	condition3 = getDistance(next_x, next_y, next_z) > 0.001;
	// streamline length in steps has been reached
	condition4 = n_steps <= g_streamLength;
	// streamline is too close to other streamlines
	condition5 = getSeparation(next_x, next_y, next_z) > 0.05;
	// check if all conditions are satisfied
	conditions = condition1 && condition2 && condition3 && condition4 && condition5;
	
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
		while (checkConditions(next_x, next_y, next_z, n_steps)) {
			point.nextX = next_x;
			point.nextY = next_y;
			point.nextZ = next_z;
			streamline.push_back(point);
			n_steps++;
			vectorFunction(next_x, next_y, next_z, vx, vy, vz, magnitude);
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
		while (checkConditions(next_x, next_y, next_z, n_steps)) {
			point.nextX = next_x;
			point.nextY = next_y;
			point.nextZ = next_z;
			streamline.push_back(point);
			n_steps++;
			vectorFunction(next_x, next_y, next_z, vx, vy, vz, magnitude);
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
		while (checkConditions(next_x1, next_y1, next_z1, n_steps) && checkConditions(next_x2, next_y2, next_z2, n_steps)) {
			point1.nextX = next_x1;
			point1.nextY = next_y1;
			point1.nextZ = next_z1;
			streamr1.push_back(point1);
			point2.nextX = next_x2;
			point2.nextY = next_y2;
			point2.nextZ = next_z2;
			streamr2.push_back(point2);
			n_steps++;
			vectorFunction(next_x1, next_y1, next_z1, vx, vy, vz, magnitude);
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
		while (checkConditions(next_x1, next_y1, next_z1, n_steps) && checkConditions(next_x2, next_y2, next_z2, n_steps)) {
			point1.nextX = next_x1;
			point1.nextY = next_y1;
			point1.nextZ = next_z1;
			streamr1.push_back(point1);
			point2.nextX = next_x2;
			point2.nextY = next_y2;
			point2.nextZ = next_z2;
			streamr2.push_back(point2);
			n_steps++;
			vectorFunction(next_x1, next_y1, next_z1, vx, vy, vz, magnitude);
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


//void TW_CALL nContours(void *ClientData) { // first draws one contour based on the scalar value chosen and then draws the n contours equally spaced across the scalar range
//	if (isPoly == 1) {
//		isocontours.clear();
//		computeContours(g_sprime); // draw the first contour based on the iso value chosen
//		for (int i = 1; i < g_ncontours; i++) { // draw the rest equally spaced out
//			float buffer = (s_max - s_min) / 1000;
//			float sprime_i = s_min + (i * (s_max - s_min) / (g_ncontours - 1));
//			if (sprime_i == s_min) sprime_i += buffer; // buffer to protect against extreme values where only 1 point may be present
//			else if (sprime_i == s_max) sprime_i -= buffer;
//			computeContours(sprime_i);
//		}
//	}
//	else if (isPoly == 0) {
//		isocontours_t.clear();
//		computeContoursTriangles(g_sprime);
//		for (int i = 1; i < g_ncontours; i++) {
//			float buffer = (s_max - s_min) / 1000;
//			float sprime_i = s_min + (i * (s_max - s_min) / (g_ncontours - 1));
//			if (sprime_i == s_min) sprime_i += buffer; // buffer to protect against extreme values where only 1 point may be present
//			else if (sprime_i == s_max) sprime_i -= buffer;
//			computeContoursTriangles(sprime_i);
//		}
//	}
//}

//void computeLIC() {
//	/*For each pixel
//		Compute a streamline using the vector field image in forward and backward direction (the streamline computation is terminated when the desired number of pixels is reached).
//		Accumulate the color values from the pixels obtained in the previous step*/
//	streamlines.clear();
//	for (int i = 0; i < IMG_RES; i++) {
//		for (int j = 0; j < IMG_RES; j++) {
//			int next_i, next_j;
//			float x, y, vx, vy, mag;
//			bool conditions;
//			streamline.clear();
//			// compute streamline in forward direction
//			x = j + 0.5;
//			y = i + 0.5;
//			next_i = i;
//			next_j = j;
//			conditions = !(next_i < 0 || next_j < 0 || next_i >= IMG_RES || next_j >= IMG_RES);
//			for (int l = 0; l < g_streamLength / 2 && conditions; l++) {
//				vx = vx_min + ((vx_max - vx_min) * vec_img[next_i][next_j][0] / 255.);
//				vy = vy_min + ((vy_max - vy_min) * vec_img[next_i][next_j][1] / 255.);
//				mag = sqrt(pow(vx, 2) + pow(vy, 2));
//				if (mag < 0.000001)
//					break;
//				vx /= mag;
//				vy /= mag;
//				x += vx;
//				y += vy;
//				next_j = (int)x;
//				next_i = (int)y;
//				conditions = !(next_i < 0 || next_j < 0 || next_i >= IMG_RES || next_j >= IMG_RES);
//				if (!conditions)
//					break;
//				streampoint point;
//				point.nextX = next_i;
//				point.nextY = next_j;
//				streamline.push_back(point);
//			}
//			// compute streamlines in backward direction
//			x = j - 0.5;
//			y = i - 0.5;
//			next_i = i;
//			next_j = j;
//			conditions = !(next_i < 0 || next_j < 0 || next_i >= IMG_RES || next_j >= IMG_RES);
//			for (int l = 0; l < g_streamLength / 2 && conditions; l++) {
//				vx = vx_min + ((vx_max - vx_min) * vec_img[next_i][next_j][0] / 255.);
//				vy = vy_min + ((vy_max - vy_min) * vec_img[next_i][next_j][1] / 255.);
//				mag = sqrt(pow(vx, 2) + pow(vy, 2));
//				if (mag < 0.000001)
//					break;
//				vx /= mag;
//				vy /= mag;
//				x -= vx;
//				y -= vy;
//				next_j = (int)x;
//				next_i = (int)y;
//				conditions = !(next_i < 0 || next_j < 0 || next_i >= IMG_RES || next_j >= IMG_RES);
//				if (!conditions)
//					break;
//				streampoint point;
//				point.nextX = next_i;
//				point.nextY = next_j;
//				streamline.push_back(point);
//			}
//			// compute the averaged color for the pixel
//			float totalNr, totalNg, totalNb;
//			totalNr = totalNg = totalNb = 0.;
//			for (int k = 0; k < streamline.size(); k++) {
//				int xp = streamline[k].nextX;
//				int yp = streamline[k].nextY;
//				float r, g, b;
//				r = (float)noise_tex[xp][yp][0];
//				g = (float)noise_tex[xp][yp][1];
//				b = (float)noise_tex[xp][yp][2];
//				totalNr += r;
//				totalNg += g;
//				totalNb += b;
//			}
//			lic_tex[i][j][0] = (unsigned char)(totalNr / streamline.size());
//			lic_tex[i][j][1] = (unsigned char)(totalNg / streamline.size());
//			lic_tex[i][j][2] = (unsigned char)(totalNb / streamline.size());
//			streamlines.push_back(streamline);
//		}
//	}
//}

//void computeELIC() {
//	// compute enhanced LIC
//	for (int i = 0; i < streamlines.size(); i++) {
//		float totalNr, totalNg, totalNb;
//		totalNr = totalNg = totalNb = 0.;
//		std::vector<streampoint> stream = streamlines[i];
//		for (int k = 0; k < stream.size(); k++) {
//			int xp = stream[k].nextX;
//			int yp = stream[k].nextY;
//			float r, g, b;
//			r = (float)lic_tex[xp][yp][0];
//			g = (float)lic_tex[xp][yp][1];
//			b = (float)lic_tex[xp][yp][2];
//			totalNr += r;
//			totalNg += g;
//			totalNb += b;
//		}
//		lic_tex_enhanced[i / IMG_RES][i % IMG_RES][0] = (unsigned char)(totalNr / stream.size());
//		lic_tex_enhanced[i / IMG_RES][i % IMG_RES][1] = (unsigned char)(totalNg / stream.size());
//		lic_tex_enhanced[i / IMG_RES][i % IMG_RES][2] = (unsigned char)(totalNb / stream.size());
//	}
//}

//void updateGradient(void* clientData) {
//	setupTwBar();
//	color3dStruct();
//}

void TW_CALL setDVRCB(const void* value, void* clientData) {
	g_enableDVR = *(const int *)value;
}

void TW_CALL getDVRCB(void* value, void* clientData) {
	*(int *)value = g_enableDVR;
}

void TW_CALL setSlicesCB(const void* value, void* clientData) {
	g_enableSlices = *(const int *)value;
}

void TW_CALL getSlicesCB(void* value, void* clientData) {
	*(int *)value = g_enableSlices;
}

void TW_CALL setSurfaceCB(const void* value, void* clientData) {
	g_isoSurfaces = *(const int *)value;
}

void TW_CALL getSurfaceCB(void* value, void* clientData) {
	*(int *)value = g_isoSurfaces;
}

void TW_CALL setTextureCB(const void* value, void* clientData) {
	g_bilinear = *(const int *)value;
}

void TW_CALL getTextureCB(void* value, void* clientData) {
	*(int *)value = g_bilinear;
}

//void TW_CALL recomputeIsoSurface(void* clientData) {
//	computeIsoSurfaces(g_sprime);
//}

//void setupTwBar() {
//	if (isPoly == 2) {
//		TwRemoveVar(bar, "Iso value");
//		TwRemoveVar(bar, "Update Iso value / no. of contours");
//		TwRemoveVar(bar, "controlSmin");
//		TwRemoveVar(bar, "controlSmax");
//		TwRemoveVar(bar, "updateMinMax");
//		TwRemoveVar(bar, "toggleXY");
//		TwRemoveVar(bar, "toggleYZ");
//		TwRemoveVar(bar, "toggleXZ");
//		TwRemoveVar(bar, "controlX");
//		TwRemoveVar(bar, "controlY");
//		TwRemoveVar(bar, "controlZ");
//		TwRemoveVar(bar, "controlGradientMin");
//		TwRemoveVar(bar, "controlGradientMax");
//		TwRemoveVar(bar, "updateGradient");
//		TwRemoveVar(bar, "toggleSurfaces");
//		TwRemoveVar(bar, "No. of Iso contours");
//		TwRemoveVar(bar, "updateTexture");
//		TwRemoveVar(bar, "updateOpacity");
//		TwRemoveVar(bar, "toggleDVR");
//		TwRemoveVar(bar, "toggleSlices");
//		// Control the range of values
//		std::string definition = "label='Increase minimum threshold' min=" + std::to_string(0.0) + " max= " + std::to_string(s_max) + " step=" + std::to_string((s_max - s_min) / 1000.);
//		const char* def = definition.c_str();
//		TwAddVarRW(bar, "controlSmin", TW_TYPE_DOUBLE, &s_min, def);
//		definition = "label='Decrease maximum threshold' min=" + std::to_string(s_min) + " max= " + std::to_string(100.0) + " step=" + std::to_string((s_max - s_min) / 1000.);
//		def = definition.c_str();
//		TwAddVarRW(bar, "controlSmax", TW_TYPE_DOUBLE, &s_max, def);
//		TwAddButton(bar, "updateMinMax", updateDataRange, NULL, "label='Update Temp. min/max'");
//		TwAddVarCB(bar, "toggleXY", TW_TYPE_BOOL32, setXYCB, getXYCB, NULL, "label='Toggle XY cutting plane'");
//		TwAddVarCB(bar, "toggleYZ", TW_TYPE_BOOL32, setYZCB, getYZCB, NULL, "label='Toggle YZ cutting plane'");
//		TwAddVarCB(bar, "toggleXZ", TW_TYPE_BOOL32, setXZCB, getXZCB, NULL, "label='Toggle XZ cutting plane'");
//		TwAddVarRW(bar, "controlX", TW_TYPE_INT32, &g_Xslice, "label='Control slice along X axis' min=0 max=49 step=1");
//		TwAddVarRW(bar, "controlY", TW_TYPE_INT32, &g_Yslice, "label='Control slice along Y axis' min=0 max=49 step=1");
//		TwAddVarRW(bar, "controlZ", TW_TYPE_INT32, &g_Zslice, "label='Control slice along Z axis' min=0 max=49 step=1");
//		definition = "label='Increase Min Gradient' min=" + std::to_string(0.0) + " max=" + std::to_string(g_gradientMax) + " step=" +
//			std::to_string((g_gradientMax - g_gradientMin) / 100.);
//		def = definition.c_str();
//		TwAddVarRW(bar, "controlGradientMin", TW_TYPE_FLOAT, &g_gradientMin, def);
//		definition = "label='Decrease Max Gradient' min=" + std::to_string(g_gradientMin) + " max=" + std::to_string(g_gradientAbsMax) + " step=" +
//			std::to_string((g_gradientMax - g_gradientMin) / 100.);
//		def = definition.c_str();
//		TwAddVarRW(bar, "controlGradientMax", TW_TYPE_FLOAT, &g_gradientMax, def);
//		TwAddButton(bar, "updateGradient", updateGradient, NULL, "label='Update Gradient limits'");
//		TwAddVarCB(bar, "toggleSlices", TW_TYPE_BOOL32, setSlicesCB, getSlicesCB, NULL, "label='Display slices'");
//		TwAddVarCB(bar, "toggleSurfaces", TW_TYPE_BOOL32, setSurfaceCB, getSurfaceCB, NULL, "label='Display IsoSurfaces'");
//		float difference = (abs_s_max - abs_s_min) / 1000; // buffer to protect against minimas and maximas where there may be only a single scalar value or a plateau
//		definition = "label='Adjust iso scalar value' min=" + std::to_string(abs_s_min) + " max=" + std::to_string(abs_s_max - difference) + " step=" + std::to_string(difference) +
//			" help='Increase/decrease iso scalar value'";
//		def = definition.c_str();
//		g_sprime = (abs_s_max + abs_s_min) / 2;
//		TwAddVarRW(bar, "Iso value", TW_TYPE_FLOAT, &g_sprime, def);
//		TwAddButton(bar, "Update Iso value / no. of contours", recomputeIsoSurface, NULL, " label = 'Load new iso surface after changing value'");
//		TwAddVarCB(bar, "toggleDVR", TW_TYPE_BOOL32, setDVRCB, getDVRCB, NULL, "label='Display 3D Volume'");
//		TwAddVarRW(bar, "updateOpacity", TW_TYPE_FLOAT, &g_opacity, "label='Update opacity' min=0 max=1 step=0.01");
//		TwAddVarCB(bar, "updateTexture", TW_TYPE_BOOL32, setTextureCB, getTextureCB, NULL, "label='Bilinear'");
//	}
//	else {
//		TwRemoveVar(bar, "toggleXY");
//		TwRemoveVar(bar, "toggleYZ");
//		TwRemoveVar(bar, "toggleXZ");
//		TwRemoveVar(bar, "controlX");
//		TwRemoveVar(bar, "controlY");
//		TwRemoveVar(bar, "controlZ");
//		TwRemoveVar(bar, "controlGradientMin");
//		TwRemoveVar(bar, "controlGradientMax");
//		TwRemoveVar(bar, "updateGradient");
//		TwRemoveVar(bar, "Iso value");
//		TwRemoveVar(bar, "Update Iso value / no. of contours");
//		TwRemoveVar(bar, "updateTexture");
//		TwRemoveVar(bar, "updateOpacity");
//		TwRemoveVar(bar, "toggleSurfaces");
//		TwRemoveVar(bar, "toggleDVR");
//		TwRemoveVar(bar, "toggleSlices");
//		//TwAddVarRW(bar, "No. of Iso contours", TW_TYPE_UINT16, &g_ncontours, " label = 'Adjust no. of contours shown' min=1 max=256 step=1 help='Increase/decrease the no. of contours'");
//		//float difference = (s_max - s_min) / 1000; // buffer to protect against minimas and maximas where there may be only a single scalar value or a plateau
//		//std::string definition = "label='Adjust iso scalar value' min=" + std::to_string(s_min + difference) + " max=" + std::to_string(s_max - difference) + " step=" + std::to_string(difference) +
//		//	" help='Increase/decrease iso scalar value'";
//		//const char* def = definition.c_str();
//		//g_sprime = (s_max + s_min) / 2;
//		//g_ncontours = 1;
//		//TwAddVarRW(bar, "Iso value", TW_TYPE_FLOAT, &g_sprime, def);
//		//TwAddButton(bar, "Update Iso value / no. of contours", nContours, NULL, " label = 'Load new iso contour after changing value or update no. of contours' ");
//
//		//// draw contours
//		//nContours(nullptr);
//	}
//}

//void updateDataRange(void* clientData) {
//	// Currently only works for temperature in 3DVis. Support for other objects to be added.
//	setupTwBar();
//	color3dStruct();
//}

void TW_CALL setXYCB(const void* value, void* clientData) {
	g_XYplane = *(const int *)value;
}

void TW_CALL getXYCB(void* value, void* clientData) {
	*(int *)value = g_XYplane;
}

void TW_CALL setYZCB(const void* value, void* clientData) {
	g_YZplane = *(const int *)value;
}

void TW_CALL getYZCB(void* value, void* clientData) {
	*(int *)value = g_YZplane;
}

void TW_CALL setXZCB(const void* value, void* clientData) {
	g_XZplane = *(const int *)value;
}

void TW_CALL getXZCB(void* value, void* clientData) {
	*(int *)value = g_XZplane;
}

void TW_CALL setArrowCB(const void* value, void* clientData) {
	g_arrows = *(const int *)value;
}

void TW_CALL getArrowCB(void* value, void* clientData) {
	*(int *)value = g_arrows;
}

void TW_CALL setStreamlinesCB(const void* value, void* clientData) {
	g_streamlines = *(const int *)value;
	g_streamribbon = 0;
}

void TW_CALL getStreamlinesCB(void* value, void* clientData) {
	*(int *)value = g_streamlines;
}

void TW_CALL setStreamribbonCB(const void* value, void* clientData) {
	g_streamribbon = *(const int *)value;
	g_streamlines = 0;
}

void TW_CALL getStreamribbonCB(void* value, void* clientData) {
	*(int *)value = g_streamribbon;
}

void TW_CALL setEnhancedCB(const void* value, void* clientData) {
	g_enhanceLIC = *(const int *)value;
}

void TW_CALL getEnhancedCB(void* value, void* clientData) {
	*(int *)value = g_enhanceLIC;
}

void TW_CALL getColorLICCB(void* value, void* clientData) {
	*(int *)value = g_coloredLIC;
}

void TW_CALL setColorLICCB(const void* value, void* clientData) {
	g_coloredLIC = *(const int *)value;
}

void TW_CALL getColorCB(void* value, void* clientData) {
	*(int *)value = g_colorPlot;
}

void TW_CALL setColorCB(const void* value, void* clientData) {
	g_colorPlot = *(const int *)value;
}

//void renderVecImg() {
//	glViewport(0, 0, (GLsizei)IMG_RES, (GLsizei)IMG_RES);
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//	gluOrtho2D(0, 1, 0, 1);
//	glClear(GL_COLOR_BUFFER_BIT);
//	glDrawBuffer(GL_BACK);
//	int i, j;
//	// render the mesh
//	for (i = 0; i < poly->ntris; i++) {
//		Triangle *temp_t = poly->tlist[i];
//		float rgb[3];
//		glBegin(GL_TRIANGLES);
//		for (j = 0; j < 3; j++)
//		{
//			Vertex *v = temp_t->verts[j];
//			//determine the color for this vertex based on its vector value
//			rgb[0] = (v->vx - vx_min) / (vx_max - vx_min);
//			rgb[1] = (v->vy - vy_min) / (vy_max - vy_min);
//			rgb[2] = 0.;
//			glColor3f(rgb[0], rgb[1], rgb[2]);
//			glVertex2f(v->x, v->y);
//		}
//		glEnd();
//	}
//	// save the rendered image into the vec_img
//	glReadBuffer(GL_BACK);
//	glReadPixels(0, 0, IMG_RES, IMG_RES, GL_RGB, GL_UNSIGNED_BYTE, vec_img);
//}

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
	/*char object_name[128] = "Field 1";*/

	switch (g_CurrentShape) {
	case 0:
		/*strcpy(object_name, "Field 1");*/
		currentField = 1;
		max_ptr = &max_sv1;
		min_ptr = &min_sv1;
		break;

	case 1:
		/*strcpy(object_name, "Field 2");*/
		currentField = 2;
		max_ptr = &max_sv2;
		min_ptr = &min_sv2;
		break;

	case 2:
		/*strcpy(object_name, "Field 3");*/
		currentField = 3;
		max_ptr = &max_sv3;
		min_ptr = &min_sv3;
		break;/*

	case 3:
		strcpy(object_name, "bnoise");
		break;

	case 4:
		strcpy(object_name, "vnoise");
		break;*/

		/*case 5:
			strcpy(object_name, "temperature1.dat");
			break;

		case 6:
			strcpy(object_name, "temperature2.dat");
			break;
		case 7:
			strcpy(object_name, "3dVis");*/
	}
	setVecFieldPointer();
	computeStreamBunch();
	computeStreamribbon();

	//if (isPoly == 0)
		//poly->finalize();
		//Reset();
		//char tmp_str[512];
		//if (strstr(object_name, "dat")) {
		//	// load dat file and put in same data structure
		//	sprintf(tmp_str, "./models/%s", object_name);
		//	isPoly = 1;
		//	Load_data_on_uniformGrids(tmp_str);
		//	build_edge_list();
		//	build_face_list();
		//}
		//else if (strstr(object_name, "3dVis")) {
		//	isPoly = 2;
		//	populate3dStruct();
		//	computeGradient();
		//	calcLimits();
		//	color3dStruct();
		//	computeIsoSurfaces(g_sprime);
		//	setupTwBar();
		//}
		//else {
		//	sprintf(tmp_str, "./models/%s.ply", object_name);
		//	FILE *this_file = fopen(tmp_str, "r");
		//	poly = new Polyhedron(this_file);
		//	fclose(this_file);
		//	isPoly = 0;
		//	poly->initialize(); // initialize everything
		//	poly->calc_bounding_sphere();
		//	poly->calc_face_normals_and_area();
		//	poly->average_normals();
		//	calcLimits();
		//	/*setColorFunction();
		//	setExtremePointers();*/
		//	//setupTwBar();
		//	renderVecImg();
		//	computeLIC();
		//	computeELIC();
		//}

		//if (isPoly != 2) {
		//	calcLimits(); // calc s_max and s_min for the new objects
		//	setupTwBar();
		//}

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

//void updateLIC(void* clientData) {
//	computeLIC();
//	computeELIC();
//}

void recompStreamline(void* clientData) {
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
		/*ShapeEV associates Shape enum values with labels that will be displayed instead of enum values
		TwEnumVal shapeEV[NUM_SHAPES] = { { SHAPE_TEAPOT, "Teapot" },{ SHAPE_TORUS, "Torus" },{ SHAPE_CONE, "Cone" },{ BUNNY, "Bunny" } };
		TwEnumVal shapeEV[NUM_SHAPES] = { { 0, "torus_field" },{ 1, "iceland_current_field" },{ 2, "diesel_field1" },{ 3, "distance_field1" },{ 4, "distance_field2" },
		{ 5, "temperature1.dat" },{ 6, "temperature2.dat" }, {7, "3D Vis" } };
		TwEnumVal shapeEV[NUM_SHAPES] = { {0, "Dipole"}, {1, "Bruno3"}, {2, "Cnoise"}, {3, "Bnoise"}, {4, "Vnoise"} };*/

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

	// Control the range of values
	//std::string definition = "label='Increase minimum threshold' min=" + std::to_string(s_min) + " max= " + std::to_string(s_max) + " step=" + std::to_string((s_max - s_min) / 1000.);
	//const char* def = definition.c_str();
	//TwAddVarRW(bar, "controlSmin", TW_TYPE_DOUBLE, &s_min, def);
	//definition = "label='Decrease maximum threshold' min=" + std::to_string(s_min) + " max= " + std::to_string(s_max) + " step=" + std::to_string((s_max - s_min) / 1000.);
	//def = definition.c_str();
	//TwAddVarRW(bar, "controlSmax", TW_TYPE_DOUBLE, &s_max, def);
	//TwAddButton(bar, "updateMinMax", updateDataRange, NULL, "label='Update Temp. min/max'");
	////Add modifier for the no. of iso-contours computed and displayed
	//TwAddVarRW(bar, "No. of Iso contours", TW_TYPE_UINT16, &g_ncontours, " label = 'Adjust no. of contours shown' min=1 max=256 step=1 help='Increase/decrease the no. of contours'");
	////Add modifier for the iso-contour value
	//float difference = (s_max - s_min) / 1000; // buffer to protect against minimas and maximas where there may be only a single scalar value
	//definition = "label='Adjust iso scalar value' min=" + std::to_string(s_min + difference) + " max=" + std::to_string(s_max - difference) + " step=" + std::to_string(difference) +
	//	" help='Increase/decrease iso scalar value'";
	//def = definition.c_str();
	//TwAddVarRW(bar, "Iso value", TW_TYPE_FLOAT, &g_sprime, def);
	//TwAddButton(bar, "Update Iso value / no. of contours", nContours, NULL, " label = 'Load new iso contour after changing value or update no. of contours' ");	
	//TwEnumVal VectorPlotEV[4] = { {0, "Magnitude"}, {1, "Angle"}, {2, "X-component"}, {3, "Y-component"} };
	//TwType PlotType = TwDefineEnum("VectorPlotType", VectorPlotEV, 4);
	//TwAddVarRW(bar, "VectorPlot", PlotType, &whichPlot, "label='Plot Type'");
	////TwAddVarCB(bar, "colorLIC", TW_TYPE_BOOL32, setColorLICCB, getColorLICCB, NULL, "label='Toggle colored LIC'");
	//TwAddVarCB(bar, "colorPlot", TW_TYPE_BOOL32, setColorCB, getColorCB, NULL, "label='Toggle Color Plots'");
	TwAddVarCB(bar, "toggleArrows", TW_TYPE_BOOL32, setArrowCB, getArrowCB, NULL, "label='Toggle Arrows'");
	TwAddVarCB(bar, "toggleStreamlines", TW_TYPE_BOOL32, setStreamlinesCB, getStreamlinesCB, NULL, "label='Toggle Streamlines'");
	TwAddVarCB(bar, "toggleStreamribbon", TW_TYPE_BOOL32, setStreamribbonCB, getStreamribbonCB, NULL, "label='Toggle Streamribbon'");
	TwAddVarRW(bar, "modifyStreamLength", TW_TYPE_INT32, &g_streamLength, "label='Change streamline length' min=2 max=500 step=1 help='NOTE: Large values will take longer to compute'");
	TwAddVarRW(bar, "moveProbeX", TW_TYPE_FLOAT, &g_probeX, "label='Change X coordinate' min=-1.0 max=1.0 step=0.01");
	TwAddVarRW(bar, "moveProbeY", TW_TYPE_FLOAT, &g_probeY, "label='Change Y coordinate' min=-1.0 max=1.0 step=0.01");
	TwAddVarRW(bar, "moveProbeZ", TW_TYPE_FLOAT, &g_probeZ, "label='Change Z coordinate' min=-1.0 max=1.0 step=0.01");
	TwEnumVal Integrators[2] = { {0, "Euler"}, {1, "Runga-Kutta Second Order"} };
	TwType Integrator = TwDefineEnum("Integrators", Integrators, 2);
	TwAddVarRW(bar, "changeIntegrator", Integrator, &whichIntegrator, "label='Change integration method'");
	TwAddVarRW(bar, "changeDt", TW_TYPE_FLOAT, &g_step, "label='Change step size' min=0.01 max=0.5 step=0.01");
	TwAddVarRW(bar, "changeSepDist", TW_TYPE_FLOAT, &g_linedist, "label='Change dist. between lines' min=0.1 max=0.5 step=0.01");
	TwAddButton(bar, "updateProbe", recompStreamline, NULL, "label='Update Streamline/Streamribbon'");
	//TwAddVarCB(bar, "toggleEnhancedLIC", TW_TYPE_BOOL32, setEnhancedCB, getEnhancedCB, NULL, "label='Toggle Enhanced LIC'");
	//TwAddVarCB(bar, "toggleColoredLIC", TW_TYPE_BOOL32, setColorLICCB, getColorLICCB, NULL, "label='Toggle Colored LIC'");
	//TwAddButton(bar, "updateLIC", updateLIC, NULL, "label='Recompute LIC'");
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
	MainWindow = glutCreateWindow("Assignment 7 – Binoy Dalal");
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

	// setExtremePointers();
	// Load the model and data here
	//FILE *this_file = fopen("./models/dipole.ply", "r");
	//poly = new Polyhedron(this_file);
	//fclose(this_file);
	//poly->initialize(); // initialize everything
	//poly->calc_bounding_sphere();
	//poly->calc_face_normals_and_area();
	//poly->average_normals();
	//calcLimits(); // calculate limits for the default figure
	//setColorFunction();
	//setExtremePointers();
	//gen_noise_tex();
	//renderVecImg();
	//computeLIC();
	//computeELIC();
	// nContours(nullptr);

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