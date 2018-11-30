#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <windows.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "Skeleton.h"

const int BRES = 640; // Resolution of the base image
const int NRES = 64; // Resolution of the noise image
const int NNoise = 32; // No. of noise images
static GLuint noiseTexImg[NNoise]; // 32 noise images
static GLuint FtTex; // Ft image
static GLuint FsTex; // Fs image
static GLubyte noiseTex[NRES][NRES][4]; // to store a noise pattern
static GLubyte baseTex[BRES][BRES][3]; // to store the base image pattern
static int fctr = 0; // noise frame counter
static float alpha = 0.4; // alpha channel for noise images
static float tmax = BRES / NRES; // texture repeat parameter
static float dmax = 0.5; // paramter to limit texture advection
Polyhedron *poly = NULL; // 3D object loaded from ply file

void genNoiseTexSt() {
	glGenTextures(32, noiseTexImg);

	int lut[256]; // lookup table to store dynamic profiles
	int phase[NRES][NRES]; // stores the phase for each pixel in the noise pattern
	int t;

	for (int i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255; // fill lookup table with block pulse

	for (int i = 0; i < NRES; i++)
		for (int j = 0; j < NRES; j++) phase[i][j] = rand() % 256;

	for (int i = 0; i < NNoise; i++) {
		t = i * 256 / NNoise;
		for (int j = 0; j < NRES; j++) {
			for (int k = 0; k < NRES; k++) {
				noiseTex[j][k][0] =
					noiseTex[j][k][1] =
					noiseTex[j][k][2] = (GLubyte)lut[(t + phase[j][k]) % 255];
				noiseTex[j][k][3] = (GLubyte)(alpha * 255);
			}
		}
		glBindTexture(GL_TEXTURE_2D, noiseTexImg[i]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NRES, NRES, 0, GL_RGBA, GL_UNSIGNED_BYTE, noiseTex);
	}
}

void genBaseImgs() {
	// generate Ft
	glGenTextures(1, &FtTex);
	for (int i = 0; i < BRES; i++) {
		for (int j = 0; j < BRES; j++) {
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
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, BRES, BRES, 0, GL_RGB, GL_UNSIGNED_BYTE, baseTex);

	// generate Fs
	glGenTextures(1, &FsTex);
	glBindTexture(GL_TEXTURE_2D, FsTex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, BRES, BRES, 0, GL_RGB, GL_UNSIGNED_BYTE, baseTex);
}

void getAdvectedCoords(Vertex *v, float &px, float &py) {
	px = v->x - ((v->t1 / v->magnitude) * (dmax * 2));
	py = v->y - ((v->t2 / v->magnitude) * (dmax * 2));
}

void calcAdvectedTexCoords() {
	float px, py;
	for (int i = 0; i < poly->nverts; i++) {
		Vertex *v = poly->vlist[i];
		getAdvectedCoords(v, px, py);
		v->tx[0] = px;
		v->tx[1] = py;
	}
}

void drawObjectWithFt() {
	glPushMatrix();
	glBindTexture(GL_TEXTURE_2D, FtTex);
	for (int i = 0; i < poly->ntris; i++) {
		Triangle* t = poly->tlist[i];
		glBegin(GL_TRIANGLES);
		for (int j = 0; j < 3; j++) {
			Vertex *v = t->verts[j];
			glTexCoord2dv(v->tx); glVertex3f(v->x, v->y, v->z);
		}
		glEnd();
	}
	glPopMatrix();
}

void blendNoise() {
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(-1.5, -1.5, 0);
	glScalef(3, 3, 3);
	glBindTexture(GL_TEXTURE_2D, noiseTexImg[fctr % NNoise]);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDepthFunc(GL_GREATER);
	glBegin(GL_QUAD_STRIP);
	glTexCoord2f(0, 0);       glVertex3f(0, 0, -40);
	glTexCoord2f(0, tmax);    glVertex3f(0, 1, -40);
	glTexCoord2f(tmax, 0);    glVertex3f(1, 0, -40);
	glTexCoord2f(tmax, tmax); glVertex3f(1, 1, -40);
	glEnd();
	glDisable(GL_BLEND);
	glDepthFunc(GL_LESS);
	fctr++;
	glPopMatrix();
}

void drawObjectWithNoise() {
	glPushMatrix();
	glBindTexture(GL_TEXTURE_2D, noiseTexImg[fctr % NNoise]);
	for (int i = 0; i < poly->ntris; i++) {
		Triangle* t = poly->tlist[i];
		glBegin(GL_TRIANGLES);
		for (int j = 0; j < 3; j++) {
			Vertex *v = t->verts[j];
			glTexCoord2dv(v->tx); glVertex3f(v->x, v->y, v->z);
		}
		glEnd();
	}
	fctr++;
	glPopMatrix();
}

void blendWithFt() {
	glPushMatrix();
	glEnable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_GREATER);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glBindTexture(GL_TEXTURE_2D, FtTex);
	glBegin(GL_QUAD_STRIP);
	// draw the image across the entire screen
	glTexCoord2f(0, 0); glVertex3f(-2, -1.5, -39);
	glTexCoord2f(0, 1); glVertex3f(-2, +1.5, -39);
	glTexCoord2f(1, 0); glVertex3f(+2, -1.5, -39);
	glTexCoord2f(1, 1); glVertex3f(+2, +1.5, -39);
	glEnd();
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glDepthFunc(GL_LESS);
	glDisable(GL_BLEND);
	glPopMatrix();
}

void blendWithShading() {
	glPushMatrix();
	GLfloat ambient[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat diffuse[] = { 0.8, .8, 1., 1.0 };
	GLfloat specular[] = { 0.8, 0.8, 1.0, 1.0 };
	glShadeModel(GL_SMOOTH);
	/*glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_NORMALIZE);*/
	/*glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);*/
	glBindTexture(GL_TEXTURE_2D, FtTex);
	/*glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	glMateriali(GL_FRONT_AND_BACK, GL_SHININESS, 80);*/

	for (int i = 0; i < poly->ntris; i++) {
		Triangle *t = poly->tlist[i];
		glBegin(GL_TRIANGLES);
		for (int j = 0; j < 3; j++) {
			Vertex *v = t->verts[j];
			glNormal3dv(v->normal.entry);
			glTexCoord2f(v->tx[0], v->tx[1]);
			glVertex3f(v->x, v->y, v->z);
		}
		glEnd();
	}

	glDisable(GL_TEXTURE_2D);
	glPopMatrix();
}

void setView() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1., 1., -1., 1., -5., 1000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void display(void) {
	GLfloat ambient[] = { 0.3, 0.3, 0.3, 1.0 };
	//GLfloat diffuse[] = { 1., 0.6, 0, 1.0};
	GLfloat diffuse[] = { 0.8, .8, 1., 1.0 };
	GLfloat specular[] = { 0.8, 0.8, 1.0, 1.0 };
	int shiny = 100;

	glClearColor(0, 0, 0, 1);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glShadeModel(GL_FLAT);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_TEXTURE_2D);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	setView();

	calcAdvectedTexCoords(); // calculate texture advection based on flow

	/* draw object with base image Ft */
	drawObjectWithFt();

	//blendWithFt();

	///* blend noise */
	blendNoise();

	// is blank for some reason - seems to be working now though there's strange behavior for 800px base image
	glReadPixels(0, 0, BRES, BRES, GL_RGB, GL_UNSIGNED_BYTE, baseTex);
	glBindTexture(GL_TEXTURE_2D, FtTex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, BRES, BRES, 0, GL_RGB, GL_UNSIGNED_BYTE, baseTex);

	// blend with shading
	blendWithShading();

	//drawObjectWithNoise();

	glFlush(); // flush frame buffer contents to screen
	glutPostRedisplay(); // call display again

	Sleep(100000);
}

void reshape(int width, int height) {
	glViewport(0, 0, (GLsizei)width, (GLsizei)height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-2., 2., -1.5, 1.5, -5., 1000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

int main(int argc, char* argv[]) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(BRES, BRES);
	glutCreateWindow(argv[0]);

	genNoiseTexSt(); //  generate noise textures
	genBaseImgs(); // generate base images

	glutDisplayFunc(display);
	glutReshapeFunc(reshape);

	FILE *this_file = fopen("./models/sphere1.ply", "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);
	poly->initialize(); // initialize everything
	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();
	poly->calc_vector_projection();

	glEnable(GL_TEXTURE_2D);
	glEnable(GL_DEPTH_TEST);

	glutMainLoop();
}