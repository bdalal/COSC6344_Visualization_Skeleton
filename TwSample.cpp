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
//  @date       2018/09/07
//  ---------------------------------------------------------------------------


#include <AntTweakBar.h>

#include <stdlib.h>
#include <stdio.h>

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

#define NUM_SHAPES 4
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
float g_MatAmbient[] = { 0.5f, 0.5f, 0.5f, 0.5f };
float g_MatDiffuse[] = { 1.0f, 1.0f, 1.0f, 0.5f };
// Light parameter
float g_LightMultiplier = 1.0f;
float g_LightDirection[] = { -0.57735f, -0.57735f, -0.57735f };
float s_min, s_max, s_mid;
// White threshold
float g_WhiteThreshold = 0.5;


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

void(*colorFunction)(float, float[]); // pointer to color function of choice

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

int g_Axes = 0;   // Toggle Axes
int g_Box = 0;    // Toggle Box

const int IMG_RES_NOISE = 512; // resolution of the noise images
unsigned char noise_tex[IMG_RES_NOISE][IMG_RES_NOISE][4];

const int IMG_RES = 512; // resolution of image for IBFV
unsigned char ft[IMG_RES][IMG_RES][4]; // texture for Ft image in first step of algo.
unsigned char f[IMG_RES][IMG_RES][4]; // texture for F image in second step of algo.

const int N_Noise = 32; // no. of noise images

float alpha = (0.12 * 255); // Blending factor for noise images

int g_color = 0;

int frame_counter = 0; // counter for noise images

float SCALE = 2.0;

#define dmax  SCALE / IMG_RES // limit texture warping

#define tmax  IMG_RES / (SCALE * IMG_RES_NOISE)

float objXmat[16]; // to store the model view matrix for transforming the texture

#include "Skeleton.h"
Polyhedron *poly = NULL;


// TODO - add arrows to check correctness
// TODO - make dmax, tmax, scale configurable

// PROBLEMS - Aspect ratio gets kind of messed up for bigger screen sizes

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
void HsvRgb(float hsv[3], float rgb[3])
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

void
Axes(float length)
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

void
InitAxesLists(void)
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
	float t = (s - s_min) / (s_max - s_min);
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
	float t = (s - s_min) / (s_max - s_min);
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

// Q 3.1
void HeatMap(float s, float rgb[3]) {
	float t = (s - s_min) / (s_max - s_min);
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
	int t = floor((s - s_min) / (s_max - s_min) * 10);
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
	float t = (sqrt(s) - sqrt(s_min)) / (sqrt(s_max) - sqrt(s_min));
	float hsv[4];
	hsv[1] = hsv[2] = hsv[3] = 1.;
	hsv[0] = (1 - t) * 240;
	HsvRgb(hsv, rgb);
}

void calcLimits() {
	s_max = s_min = poly->tlist[0]->verts[0]->magnitude;
	for (int i = 0; i < poly->ntris; i++) {
		Triangle *temp_t = poly->tlist[i];
		for (int j = 0; j < 3; j++) {
			Vertex *temp_v = temp_t->verts[j];
			float s = temp_v->magnitude;
			if (s > s_max)
				s_max = s;
			if (s < s_min)
				s_min = s;
		}
	}
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

// TODO provide option to use either modelview advection or plain old advection
void getDistortedVertices(Vertex* temp_v, double &px, double &py) {
	float p[4], pr[4];
	float x = temp_v->x;
	float y = temp_v->y;
	float z = temp_v->z;
	float vx = temp_v->nx / temp_v->magnitude;
	float vy = temp_v->ny / temp_v->magnitude;
	float vz = temp_v->nz / temp_v->magnitude;
	pr[0] = pr[1] = pr[2] = pr[3] = 0;
	vx *= dmax * 2;
	vy *= dmax * 2;
	vz *= dmax * 2;
	p[0] = x - vx;
	p[1] = y - vy;
	p[2] = z - vz;
	p[3] = 1;
	glGetFloatv(GL_MODELVIEW_MATRIX, objXmat);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			pr[i] += objXmat[i + 4 * j] * p[j];
	px = pr[0];
	py = pr[1];
	/*px = x + vx;
	py = x + vy;*/
}

// Algo step 1 - Initialize Ft with background color (gray)
void initFt() {
	for (int i = 0; i < IMG_RES; i++) {
		for (int j = 0; j < IMG_RES; j++) {
			ft[i][j][0] = ft[i][j][1] = ft[i][j][2] = (unsigned char)128;
			ft[i][j][3] = (unsigned char)(0.4 * 255);
		}
	}
}

// Algo step 2 - Initialize F with background color (gray)
void initF() {
	for (int i = 0; i < IMG_RES; i++) {
		for (int j = 0; j < IMG_RES; j++) {
			f[i][j][0] = f[i][j][1] = f[i][j][2] = (unsigned char)128;
			f[i][j][3] = (unsigned char)(0.5 * 255);
		}
	}
}


// Generate noise textures
void gen_noise_tex()
{
	for (int n = 0; n < N_Noise; n++) {
		for (int x = 0; x < IMG_RES_NOISE; x++)
			for (int y = 0; y < IMG_RES_NOISE; y++)
			{
				noise_tex[x][y][0] =
					noise_tex[x][y][1] =
					noise_tex[x][y][2] = (unsigned char)255 * (rand() % 32768) / 32768.0;
				noise_tex[x][y][3] = (unsigned char)(0.5 * 255);
			}
		glNewList(n + 1, GL_COMPILE);
		glTexImage2D(GL_TEXTURE_2D, 0, 4, IMG_RES_NOISE, IMG_RES_NOISE, 0, GL_RGBA, GL_UNSIGNED_BYTE, noise_tex);
		glEndList();
	}
}

void gen_noise_tex_st() {
	int lut[256]; // lookup table to store dynamic profiles
	int phase[IMG_RES_NOISE][IMG_RES_NOISE]; // stores the phase for each pixel in the noise pattern
	int t;

	for (int i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255; // fill lookup table with block pulse

	for (int i = 0; i < IMG_RES_NOISE; i++)
		for (int j = 0; j < IMG_RES_NOISE; j++) phase[i][j] = rand() % 256;

	for (int n = 0; n < N_Noise; n++) {
		t = n * 256 / N_Noise;
		for (int x = 0; x < IMG_RES_NOISE; x++)
			for (int y = 0; y < IMG_RES_NOISE; y++)
			{
				noise_tex[x][y][0] =
					noise_tex[x][y][1] =
					noise_tex[x][y][2] = (unsigned char)lut[(t + phase[x][y]) % 255];
				noise_tex[x][y][3] = (unsigned char)(0.5 * 255);
			}
		glNewList(n + 1, GL_COMPILE);
		glTexImage2D(GL_TEXTURE_2D, 0, 4, IMG_RES_NOISE, IMG_RES_NOISE, 0, GL_RGBA, GL_UNSIGNED_BYTE, noise_tex);
		glEndList();
	}
}

void set_world_space_view(float mat[]) {
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
}

void reset_matrices() {
	// reset the projection matrix
	// float pr_mat[16];
	float mv_mat[16];
	// GLint m_mode;
	// glGetFloatv(GL_PROJECTION_MATRIX, pr_mat);
	glGetFloatv(GL_MODELVIEW_MATRIX, mv_mat);
	// glGetIntegerv(GL_MATRIX_MODE, &m_mode);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// glMultMatrixf(pr_mat);
	// glGetIntegerv(GL_MATRIX_MODE, &m_mode);


	//glGetFloatv(GL_PROJECTION_MATRIX, pr_mat);

	// set the orthographic view and clipping planes
	glOrtho(-2.5, 2.5, -2.5, 2.5, -1000.0, 4000.0);
	// gluPerspective(45.0, 1.0, 0.1, 40.0);

	//glGetFloatv(GL_PROJECTION_MATRIX, pr_mat);

	// reset the modelview matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixf(mv_mat);
}

void noise_blend_test() {
	glClearColor(0, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, IMG_RES, IMG_RES, 0, GL_RGBA, GL_UNSIGNED_BYTE, ft);
	glBegin(GL_QUAD_STRIP);
	glTexCoord2f(0.0, 0.0); glVertex2f(0, 0);
	glTexCoord2f(0.0, 1.0); glVertex2f(1, 0);
	glTexCoord2f(1.0, 0.0); glVertex2f(0, 1);
	glTexCoord2f(1.0, 1.0); glVertex2f(1, 1);
	glEnd();
}

// Callback function called by GLUT to render screen
void Display(void)
{
	float v[4]; // will be used to set light parameters
	float mat[4 * 4]; // rotation matrix
	Sleep(50);
	// Clear frame buffer
	glClearColor(0, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glEnable(GL_NORMALIZE);

	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glShadeModel(GL_FLAT);

	// Set light
	//glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);
	//v[0] = v[1] = v[2] = g_LightMultiplier * 0.4f; v[3] = 1.0f;
	//glLightfv(GL_LIGHT0, GL_AMBIENT, v);
	//v[0] = v[1] = v[2] = g_LightMultiplier * 0.8f; v[3] = 1.0f;
	//glLightfv(GL_LIGHT0, GL_DIFFUSE, v);
	//v[0] = -g_LightDirection[0]; v[1] = -g_LightDirection[1]; v[2] = -g_LightDirection[2]; v[3] = 0.0f;
	//glLightfv(GL_LIGHT0, GL_POSITION, v);

	//// Set material
	//glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, g_MatAmbient);
	//glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, g_MatDiffuse);

	// Setup texutring params and enable texturing
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glEnable(GL_TEXTURE_2D);
		
	// get base texture image
	// glDrawBuffer(GL_BACK);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, IMG_RES, IMG_RES, 0, GL_RGBA, GL_UNSIGNED_BYTE, ft);
	glGetFloatv(GL_MODELVIEW_MATRIX, objXmat);
	// reset transformation matrices from previous executions
	reset_matrices();
	//glGetFloatv(GL_MODELVIEW_MATRIX, objXmat);
	/*float pr_mat[16];
	float mv_mat[16];
	glGetFloatv(GL_PROJECTION_MATRIX, pr_mat);
	glGetFloatv(GL_MODELVIEW_MATRIX, mv_mat);*/
	// go to world space
	glPushMatrix();
	//glGetFloatv(GL_MODELVIEW_MATRIX, objXmat);


	// set the world space view to enable rotation etc.
	set_world_space_view(mat);
	/*glGetFloatv(GL_PROJECTION_MATRIX, pr_mat);
	glGetFloatv(GL_MODELVIEW_MATRIX, mv_mat);*/
	// choose color mapping for object
	// setColorFunction();

	// Currently in world space

	// advect base texture and draw object

	// advect texture
	double px, py;
	for (int i = 0; i < poly->ntris; i++) {
		Triangle *t = poly->tlist[i];
		glBegin(GL_POLYGON);
		for (int j = 0; j < 3; j++) {
			Vertex *v = t->verts[j];
			// get advected coordinates for texture
			getDistortedVertices(v, px, py);
			v->tx[0] = px;
			v->tx[1] = py;
			// map advected coordinates to the object
			glTexCoord2dv(v->tx);
			glVertex3f(v->x, v->y, v->z);
		}
		glEnd();
	}
	// world space ends here - pop matrix to go to view space and blend the noise texture
	glPopMatrix();

	/*glGetFloatv(GL_PROJECTION_MATRIX, pr_mat);
	glGetFloatv(GL_MODELVIEW_MATRIX, mv_mat);*/

	//increment noise frame counter
	frame_counter++;

	// enable noise blending
	glEnable(GL_BLEND);
	// specify blend function
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	// specify depth function for depth testing
	glDepthFunc(GL_GREATER);
	// get noise image based on frame counter
	glCallList(frame_counter % N_Noise + 1);
	// draw the noise texture in the view space
	glBegin(GL_QUAD_STRIP);
	glTexCoord2f(0.0, 0.0);   glVertex3f(-5.0, -5.0, -39.9);
	glTexCoord2f(0.0, tmax);  glVertex3f(-5.0, 5.0, -39.9);
	glTexCoord2f(tmax, 0.0);  glVertex3f(5.0, -5.0, -39.9);
	glTexCoord2f(tmax, tmax); glVertex3f(5.0, 5.0, -39.9);
	glEnd();
	glDepthFunc(GL_LESS);
	// blending done - disable it
	glDisable(GL_BLEND);

	// read the advected and blended texture into Ft
	// glReadBuffer(GL_BACK);
	// glReadPixels(0, 0, IMG_RES, IMG_RES, GL_RGBA, GL_UNSIGNED_BYTE, ft);
	glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, ft);
	// glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, IMG_RES, IMG_RES, 0);

	// all drawing finished - go back to world space and draw the shading, axes and the tweak bar
	glPushMatrix();

	// set the world space view to enable rotation etc.
	set_world_space_view(mat);

	// set texture mode to modulate to be able to draw the shaded texture on top of the original
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	///*glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glShadeModel(GL_SMOOTH);*/
	glShadeModel(GL_SMOOTH);
	////// create new texture for blending with the Ft saved earlier
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, IMG_RES, IMG_RES, 0, GL_RGBA, GL_UNSIGNED_BYTE, ft);
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
	
	// draw object with normals and texture mapping
	for (int i = 0; i < poly->ntris; i++) {
		Triangle *t = poly->tlist[i];
		glBegin(GL_POLYGON);
		for (int j = 0; j < 3; j++) {
			Vertex *v = t->verts[j];
			glNormal3dv(v->normal.entry);
			glTexCoord2dv(v->tx);
			// glColor3f(0.5, 0.5, 0.5);
			glVertex3f(v->x, v->y, v->z);
		}
		glEnd();
	}
	// noise_blend_test();
	glPopMatrix();

	/*glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, IMG_RES, IMG_RES, 0, GL_RGB, GL_UNSIGNED_BYTE, f);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glBegin(GL_QUAD_STRIP);
	glTexCoord2f(0.0, 0.0);  glVertex3f(-5.0, -5.0, -5.0);
	glTexCoord2f(0.0, 1.0); glVertex3f(-5.0, 5.0, -5.0);
	glTexCoord2f(1.0, 0.0);  glVertex3f(5.0, -5.0, -5.0);
	glTexCoord2f(1.0, 1.0); glVertex3f(5.0, 5.0, -5.0);
	glEnd();
	glDisable(GL_BLEND);*/

	// done texturing - disable it.
	glDisable(GL_TEXTURE);

	glPushMatrix();
	set_world_space_view(mat);
	// Draw axes
	if (g_Axes)
		glCallList(AxesList);

	//	// Draw tweak bars
	TwDraw();
	// done drawing axes and tweak bar - pop matrix to go to view space
	glPopMatrix();
	// Present frame buffer
	// glutSwapBuffers();
	glFlush();

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


void TW_CALL setColorCB(const void *value, void *clientData)
{
	g_color = *(const int *)value; // copy value to g_Axes
}

void TW_CALL getColorCB(void *value, void *clientData)
{
	*(int *)value = g_color; // copy g_AutoRotate to value
}


void TW_CALL loadNewObjCB(void *clientData)
{
	char object_name[128] = "Bunny";

	switch (g_CurrentShape) {
	case 0:
		strcpy(object_name, "bunny1");
		break;

	case 1:
		strcpy(object_name, "sphere1");
		break;

	case 2:
		strcpy(object_name, "torus1");
		break;

	case 3:
		strcpy(object_name, "torus2");
		break;
	}

	poly->finalize();

	//Reset();

	char tmp_str[512];

	sprintf(tmp_str, "./models/%s.ply", object_name);

	FILE *this_file = fopen(tmp_str, "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);

	poly->initialize(); // initialize everything

	// Q 3.3a	
	calcLimits(); // calc s_max and s_min for the new objects
	g_WhiteThreshold = 0.5; // reset g_WhiteThreshold

	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();

	glutSetWindow(MainWindow);
	glutPostRedisplay();
}


void InitTwBar(TwBar *bar)
{
	// Create a tweak bar
	bar = TwNewBar("TweakBar");
	TwDefine(" GLOBAL help='This example shows how to integrate AntTweakBar with GLUT and OpenGL.' "); // Message added to the help bar.
	TwDefine(" TweakBar size='200 800' color='0 128 255' alpha=128  "); // change default tweak bar size and color
	TwDefine(" TweakBar  label='Visual Parameters'");        // change the title of the Tweakbar

	// Add callback to toggle reference axes (callback functions are defined above).
	TwAddVarCB(bar, "Axes", TW_TYPE_BOOL32, SetAxesCB, GetAxesCB, NULL,
		" label='Axes' key=a help='Toggle reference axes.' ");

	TwAddSeparator(bar, NULL, NULL);

	// Add 'g_Zoom' to 'bar': this is a modifable (RW) variable of type TW_TYPE_FLOAT. Its key shortcuts are [z] and [Z].
	TwAddVarRW(bar, "Zoom", TW_TYPE_FLOAT, &g_Zoom,
		" min=0.01 max=2.5 step=0.01 keyIncr=z keyDecr=Z help='Scale the object (1=original size).' ");

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
	// TwAddVarRW(bar, "Ambient", TW_TYPE_COLOR3F, &g_MatAmbient, " group='Material' ");

	// Add 'g_MatDiffuse' to 'bar': this is a variable of type TW_TYPE_COLOR3F (3 floats color, alpha is ignored)
	// and is inserted into group 'Material'.
	// TwAddVarRW(bar, "Diffuse", TW_TYPE_COLOR3F, &g_MatDiffuse, " group='Material' ");

	// TwAddSeparator(bar, " objects ", NULL);

	// Add the enum variable 'g_CurrentShape' to 'bar'
	// (before adding an enum variable, its enum type must be declared to AntTweakBar as follow)
	{
		// ShapeEV associates Shape enum values with labels that will be displayed instead of enum values
		//TwEnumVal shapeEV[NUM_SHAPES] = { { SHAPE_TEAPOT, "Teapot" },{ SHAPE_TORUS, "Torus" },{ SHAPE_CONE, "Cone" },{ BUNNY, "Bunny" } };

		TwEnumVal shapeEV[NUM_SHAPES] = { { 0, "Bunny" },{ 1, "Sphere" },{ 2, "Torus - 1" },{ 3, "Torus - 2" } };

		// Create a type for the enum shapeEV
		TwType shapeType = TwDefineEnum("ShapeType", shapeEV, NUM_SHAPES);

		// add 'g_CurrentShape' to 'bar': this is a variable of type ShapeType. Its key shortcuts are [<] and [>].
		TwAddVarRW(bar, "Shape", shapeType, &g_CurrentShape, " keyIncr='<' keyDecr='>' help='Change object shape.' ");

		// add a button to reload the selected object
		TwAddButton(bar, "Update (Re-load)", loadNewObjCB, NULL, " label='Load new object after selection' ");
	}

	TwAddSeparator(bar, " others ", NULL);

	// Enable/disable coloring
	TwAddVarCB(bar, "enableColor", TW_TYPE_BOOL32, setColorCB, getColorCB, NULL, "label='Enable/Disable Coloring of object'");

	// Add the enum variable 'whichColor' to 'bar' 
	{
		TwEnumVal ColorEV[NUM_COLORS] = { {0, "Rainbow"}, {1, "Blue-White-Red"}, {2, "Heat map"}, {3, "Discrete"}, {4, "NonLinear - Extremes"} };
		// Create a type for the enum ColorEV
		TwType ColorType = TwDefineEnum("ColoType", ColorEV, NUM_COLORS);

		// add 'whichColor' to 'bar': this is a variable of type ShapeType. Its key shortcuts are [+] and [-].
		TwAddVarRW(bar, "Object colors", ColorType, &whichColor, " help='Change object color.' ");
	}

	TwAddVarRW(bar, "WhiteThreshold", TW_TYPE_FLOAT, &g_WhiteThreshold,
		" label = 'Adjust white threshold' min=0 max=1 step=0.01 keyIncr = 'w' keyDecr = 's' help='Increase/decrease white threshold' ");

	//TwAddVarRW(bar, "modifyTexWarp", TW_TYPE_FLOAT, &dmax, "label='Modify degree of texture warping' min=0.01 max=10 step=0.05");

	//TwAddVarRW(bar, "modifyTexMap", TW_TYPE_FLOAT, &tmax, "label='Modify texture UV' min=1. max=10. step=0.5");
	
	TwAddVarRW(bar, "modifyScale", TW_TYPE_FLOAT, &SCALE, "label='Modify Scale' min=0.0001 max=512. step=0.1");
}

// TODO enable colors
// TODO light object from all sides

//// Main
int main(int argc, char *argv[])
{
	TwBar *bar = NULL; // Pointer to the tweak bar

	float axis[] = { 0.7f, 0.7f, 0.0f }; // initial model rotation
	float angle = 0.8f;

	// Initialize GLUT
	glutInit(&argc, argv);
	// First parameter is the buffer - single/double
	// probably no noticeable difference between single and double
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(1024, 1024);
	MainWindow = glutCreateWindow("Final Project – IBFVS – Binoy Dalal (1794070)");
	glutCreateMenu(NULL);

	// Set GLUT callbacks
	glutDisplayFunc(Display); // Display is a function pointer - function contains info. about objects to be drawn on the screen
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


	// Load the model and data here
	FILE *this_file = fopen("./models/bunny1.ply", "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);
	poly->initialize(); // initialize everything

	calcLimits(); // calculate s_max and s_min for the default figure
	setColorFunction();
	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();

	// generate noise textures
	gen_noise_tex();
	initFt();

	// Build a display list for the axes

	InitAxesLists();

	// Initialize the AntTweakBar interface

	InitTwBar(bar);

	// Store time
	g_RotateTime = GetTimeMs();
	// Init rotation
	SetQuaternionFromAxisAngle(axis, angle, g_Rotation);
	SetQuaternionFromAxisAngle(axis, angle, g_RotateStart);

	// Call the GLUT main loop
	glutMainLoop(); // Constantly renders the object on screen

	return 0;
}
