#include <Windows.h>
#include <AntTweakBar.h>
#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "Skeleton.h"

#define NUM_SHAPES 7

int MainWindow;
const int IMG_RES = 512;
const int NOISE_RES = 64;
const int NNoise = 32; // No. of noise images
const float BOXSIZE = 2.f;
const float AXES_WIDTH = 3.;
const float LENFRAC = 0.10f; // fraction of the length to use as height of the characters
const float BASEFRAC = 1.10f; // fraction of length to use as start location of the characters
static float alpha = 0.4; // alpha channel for noise images
static float xx[] = {0.f, 1.f, 0.f, 1.f};
static float xy[] = {-.5f, .5f, .5f, -.5f};
static int xorder[] = {1, 2, -3, 4};
static float yx[] = {0.f, 0.f, -.5f, .5f};
static float yy[] = {0.f, .6f, 1.f, 1.f};
static int yorder[] = {1, 2, 3, -2, 4};
static float zx[] = {1.f, 0.f, 1.f, 0.f, .25f, .75f};
static float zy[] = {.5f, .5f, -.5f, -.5f, 0.f, 0.f};
static int zorder[] = {1, 2, 3, 4, -5, 6};
float g_Zoom = 1.0f;
float g_WhiteThreshold = 0.5;
float g_Rotation[] = { 0.0f, 0.0f, 0.0f, 1.0f };
float g_RotateStart[] = { 0.0f, 0.0f, 0.0f, 1.0f };
float g_LightMultiplier = 1.0f;
float g_LightDirection[] = { -0.57735f, -0.57735f, -0.57735f };
float SCALE = 0.00225;
float dmax = 0.5;
static float tmax = 2 * IMG_RES / (1. * NOISE_RES); // texture repeat parameter
int g_AutoRotate = 0;
int g_RotateTime = 0;
int g_CurrentShape = 0;
int g_Axes = 0;
int g_arrows = 0;
static int fctr = 0; // noise frame counter
GLuint BoxList = 100;
GLuint AxesList = 101;
static GLuint noiseTexImg[NNoise]; // 32 noise images
static GLuint FtTex, FsTex, FtTexA; // Ft image
static GLubyte noiseTex[NOISE_RES][NOISE_RES][4]; // to store a noise pattern
static GLubyte baseTex[IMG_RES][IMG_RES][3]; // to store the base image pattern
static GLubyte baseTexA[IMG_RES][IMG_RES][4]; // to store the base image pattern
Polyhedron *poly;

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

void MultiplyQuaternions(const float *q1, const float *q2, float *qout)
{
	float qr[4];
	qr[0] = q1[3] * q2[0] + q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1];
	qr[1] = q1[3] * q2[1] + q1[1] * q2[3] + q1[2] * q2[0] - q1[0] * q2[2];
	qr[2] = q1[3] * q2[2] + q1[2] * q2[3] + q1[0] * q2[1] - q1[1] * q2[0];
	qr[3] = q1[3] * q2[3] - (q1[0] * q2[0] + q1[1] * q2[1] + q1[2] * q2[2]);
	qout[0] = qr[0]; qout[1] = qr[1]; qout[2] = qr[2]; qout[3] = qr[3];
}

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

int GetTimeMs()
{
#if !defined(_WIN32)
	return glutGet(GLUT_ELAPSED_TIME);
#else
	// glutGet(GLUT_ELAPSED_TIME) seems buggy on Windows
	return (int)GetTickCount();
#endif
}

void clearFrame() {
	glClearColor(0, 0, 0, 1);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void setupDisplayParams() {
	glShadeModel(GL_FLAT);
	glEnable(GL_DEPTH_TEST);
}

void enableTexturing() {
	glEnable(GL_TEXTURE_2D);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}

void setupView(float mat[]) {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-2, 2, -2, 2, -5, 4000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void disableTexturing() {
	glDisable(GL_TEXTURE_2D);
}

void getAdvectedCoords(Vertex *v) {
	float mat[16], p[4], pr[4];
	glGetFloatv(GL_MODELVIEW_MATRIX, mat);
	pr[0] = pr[1] = pr[2] = pr[3] = 0;
	p[0] = v->x - ((v->t1 / v->t_magnitude) * (SCALE * 2));
	p[1] = v->y - ((v->t2 / v->t_magnitude) * (SCALE * 2));
	p[2] = v->z - ((v->t3 / v->t_magnitude) * (SCALE * 2));
	p[3] = 1;
	for (int k = 0; k < 4; k++)
		for (int j = 0; j < 4; j++)
			pr[k] += mat[k + 4 * j] * p[j];
	v->tx[0] = pr[0];
	v->tx[1] = pr[1];
}

void calcTextureCoords() {
	for (int i = 0; i < poly->nverts; i++) {
		Vertex *v = poly->vlist[i];
		getAdvectedCoords(v);
	}
}

void drawObjectWithFt() {
	glDrawBuffer(GL_BACK);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, IMG_RES, IMG_RES, 0, GL_RGB, GL_UNSIGNED_BYTE, baseTex);
	for (int i = 0; i < poly->ntris; i++) {
		Triangle* t = poly->tlist[i];
		glBegin(GL_TRIANGLES);
		for (int j = 0; j < 3; j++) {
			Vertex *v = t->verts[j];
			glTexCoord2dv(v->tx); glVertex3f(v->x, v->y, v->z);
		}
		glEnd();
	}
}

void blendNoise() {
	glBindTexture(GL_TEXTURE_2D, noiseTexImg[fctr % NNoise]);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDepthFunc(GL_GREATER);
	glBegin(GL_QUAD_STRIP);
	glTexCoord2f(0, 0);       glVertex3f(-2, -2, -40);
	glTexCoord2f(0, tmax);    glVertex3f(-2,  2, -40);
	glTexCoord2f(tmax, 0);    glVertex3f( 2, -2, -40);
	glTexCoord2f(tmax, tmax); glVertex3f( 2,  2, -40);
	glEnd();
	glDisable(GL_BLEND);
	glDepthFunc(GL_LESS);
	fctr++;
}

void readImage() {
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, IMG_RES, IMG_RES, GL_RGB, GL_UNSIGNED_BYTE, baseTex);
	glBindTexture(GL_TEXTURE_2D, FtTex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, IMG_RES, IMG_RES, 0, GL_RGB, GL_UNSIGNED_BYTE, baseTex);
	glReadPixels(0, 0, IMG_RES, IMG_RES, GL_RGBA, GL_UNSIGNED_BYTE, baseTexA);
	glBindTexture(GL_TEXTURE_2D, FtTexA);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, IMG_RES, IMG_RES, 0, GL_RGBA, GL_UNSIGNED_BYTE, baseTexA);
}

void blendWithShading() {
	clearFrame();
	glShadeModel(GL_SMOOTH);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, IMG_RES, IMG_RES, 0, GL_RGB, GL_UNSIGNED_BYTE, baseTex);
	for (int i = 0; i < poly->ntris; i++) {
		Triangle *t = poly->tlist[i];
		glBegin(GL_TRIANGLES);
		for (int j = 0; j < 3; j++) {
			Vertex *v = t->verts[j];
			glNormal3dv(v->normal.entry);			glTexCoord2dv(v->tx);
			glVertex3f(v->x, v->y, v->z);
		}
		glEnd();
	}
}

void setScene(float mat[]){
	float matm[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, matm);
	glTranslatef(0.0, 0.0, -0.50);
	glGetFloatv(GL_MODELVIEW_MATRIX, matm);
	glTranslatef(poly->center.entry[0], poly->center.entry[1], poly->center.entry[2]);
	glGetFloatv(GL_MODELVIEW_MATRIX, matm);
	glGetFloatv(GL_MODELVIEW_MATRIX, matm);
	glScalef(0.7, 0.7, 0.7);
	glGetFloatv(GL_MODELVIEW_MATRIX, matm);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
	glGetFloatv(GL_MODELVIEW_MATRIX, matm);
	glGetFloatv(GL_MODELVIEW_MATRIX, matm);
}

void Display(void) {
	float mat[16];

	clearFrame();

	setupDisplayParams();

	enableTexturing();

	calcTextureCoords();
	
	setupView(mat);

	glPushMatrix();
    setScene(mat);
	drawObjectWithFt();
	glPopMatrix();

	blendNoise();

	readImage();

	blendWithShading();

	disableTexturing();

	if (g_Axes)
		glCallList(AxesList);
	TwDraw();

	glFlush();
	glutPostRedisplay();

}

void Reshape(int width, int height) {
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-2, 2, -2, 2, -5, 4000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);
	glTranslatef(0, 0.6f, -1);
	TwWindowSize(width, height);
}

void Terminate(void) {
	TwTerminate();
}

void genNoise() {
	glGenTextures(32, noiseTexImg);
	int lut[256]; // lookup table to store dynamic profiles
	int phase[NOISE_RES][NOISE_RES]; // stores the phase for each pixel in the noise pattern
	int t;
	for (int i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255; // fill lookup table with block pulse
	for (int i = 0; i < NOISE_RES; i++)
		for (int j = 0; j < NOISE_RES; j++) phase[i][j] = rand() % 256;
	for (int i = 0; i < NNoise; i++) {
		t = i * 256 / NNoise;
		for (int j = 0; j < NOISE_RES; j++) {
			for (int k = 0; k < NOISE_RES; k++) {
				noiseTex[j][k][0] =
					noiseTex[j][k][1] =
					noiseTex[j][k][2] = (GLubyte)lut[ (phase[j][k]) % 255];
				noiseTex[j][k][3] = (GLubyte)(alpha * 255);
			}
		}
		glBindTexture(GL_TEXTURE_2D, noiseTexImg[i]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NOISE_RES, NOISE_RES, 0, GL_RGBA, GL_UNSIGNED_BYTE, noiseTex);
	}
}

void initFt() {
	glGenTextures(1, &FtTex);
	for (int i = 0; i < IMG_RES; i++) {
		for (int j = 0; j < IMG_RES; j++) {
			baseTex[i][j][0] =
				baseTex[i][j][1] =
				baseTex[i][j][2] = (GLubyte)128;
		}
	}
	glBindTexture(GL_TEXTURE_2D, FtTex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, IMG_RES, IMG_RES, 0, GL_RGB, GL_UNSIGNED_BYTE, baseTex);

	glBindTexture(GL_TEXTURE_2D, FsTex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, IMG_RES, IMG_RES, 0, GL_RGB, GL_UNSIGNED_BYTE, baseTex);
}

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
	glNewList(AxesList, GL_COMPILE);
	glLineWidth(AXES_WIDTH);
	Axes(1.5);
	glLineWidth(1.);
	glEndList();
}

void TW_CALL SetAxesCB(const void *value, void *clientData)
{
	(void)clientData; // unused
	g_Axes = *(const int *)value; // copy value to g_Axes
}

void TW_CALL GetAxesCB(void *value, void *clientData)
{
	(void)clientData; // unused
	*(int *)value = g_Axes; // copy g_AutoRotate to value
}

void TW_CALL SetAutoRotateCB(const void *value, void *clientData)
{
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

void TW_CALL GetAutoRotateCB(void *value, void *clientData)
{
	*(int *)value = g_AutoRotate; // copy g_AutoRotate to value
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
	case 4:
		strcpy(object_name, "bnoise");
		break;
	case 5:
		strcpy(object_name, "cnoise");
		break;
	case 6:
		strcpy(object_name, "vnoise");
		break;
	}
	poly->finalize();
	char tmp_str[512];
	sprintf(tmp_str, "./models/%s.ply", object_name);
	FILE *this_file = fopen(tmp_str, "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);
	poly->initialize(); // initialize everything
	//calcLimits(); // calc s_max and s_min for the new objects
	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();
	poly->calc_vector_projection();
	glutSetWindow(MainWindow);
	glutPostRedisplay();
}

void TW_CALL setArrowCB(const void* value, void* clientData) {
	g_arrows = *(const int *)value;
}

void TW_CALL getArrowCB(void* value, void* clientData) {
	*(int *)value = g_arrows;
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
		" min=0.01 max=52.5 step=0.1 keyIncr=z keyDecr=Z help='Scale the object (1=original size).' ");
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
	TwAddSeparator(bar, " objects ", NULL);
	{
		TwEnumVal shapeEV[NUM_SHAPES] = { { 0, "Bunny" },{ 1, "Sphere" },{ 2, "Torus - 1" },{ 3, "Torus - 2" },{4, "BNoise"},{5,"CNoise"},{6,"VNoise"} };
		// Create a type for the enum shapeEV
		TwType shapeType = TwDefineEnum("ShapeType", shapeEV, NUM_SHAPES);
		// add 'g_CurrentShape' to 'bar': this is a variable of type ShapeType. Its key shortcuts are [<] and [>].
		TwAddVarRW(bar, "Shape", shapeType, &g_CurrentShape, " keyIncr='<' keyDecr='>' help='Change object shape.' ");
		// add a button to reload the selected object
		TwAddButton(bar, "Update (Re-load)", loadNewObjCB, NULL, " label='Load new object after selection' ");
	}
	TwAddSeparator(bar, " others ", NULL);
	TwAddVarRW(bar, "modifyScale", TW_TYPE_FLOAT, &SCALE, "label='Modify Scale' min=0.0001 max=5. step=0.001");
	TwAddVarCB(bar, "toggleArrows", TW_TYPE_BOOL32, setArrowCB, getArrowCB, NULL, "label='Toggle Arrows'");
}

int main(int argc, char *argv[]) {
	TwBar *bar = NULL; // Pointer to the tweak bar
	float axis[] = { 0.7f, 0.7f, 0.0f }; // initial model rotation
	float angle = 0.8f;
	// Initialize GLUT
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(IMG_RES, IMG_RES);
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
	FILE *this_file = fopen("./models/sphere1.ply", "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);
	poly->initialize(); // initialize everything
	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();
	poly->calc_vector_projection();
	// generate noise textures
	genNoise();
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