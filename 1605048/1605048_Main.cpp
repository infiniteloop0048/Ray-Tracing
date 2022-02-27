#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<bits/stdc++.h>
#include <windows.h>
#include <GL/glut.h>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <stack>
#include <algorithm>
#include "1605048_Header.h"

using namespace std;
#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle, angleForQW, angleForAS, angleForDF;
double p1[16], p2[16], p0[16];
double twod[4][4];




/* class vectorr
{
	public:
    double x, y, z;
		vectorr operator*(double val) {
        vectorr temp;
        temp.x = x * val;
        temp.y = y * val;
        temp.z = z * val;
        return temp;
    }
    Vector3D operator+(vectorr vect) {
        Vector3D temp;
        temp.x = vect.x + x;
        temp.y = vect.y + y;
        temp.z = vect.z + z;
        return temp;
    }
    Vector3D operator-(vectorr vect) {
        Vector3D temp;
        temp.x = vect.x - x;
        temp.y = vect.y - y;
        temp.z = vect.z - z;
        return temp;
    }
}; */

/* struct vectorr
{
    double x, y, z;
}; */

point p1prime, p2prime;
vector<point> pointlist;

void forwardAlongLaxis(double m){
    pos.x += m * l.x;
    pos.y += m * l.y;
    pos.z += m * l.z;
}

void forwardAlongRaxis(double m){
    pos.x += m * r.x;
    pos.y += m * r.y;
    pos.z += m * r.z;
}

void forwardAlongUaxis(double m){
    pos.x += m * u.x;
    pos.y += m * u.y;
    pos.z += m * u.z;
}

Vector3D crossProduct(Vector3D a, Vector3D b){
    Vector3D c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}


Vector3D addTwoVector(Vector3D a, Vector3D b){
    Vector3D sum;
    sum.x = a.x + b.x;
    sum.y = a.y + b.y;
    sum.z = a.z + b.z;
    return sum;
}

Vector3D multiplyByScalar(Vector3D v,double scalar){
    v.x *= scalar;
    v.y *= scalar;
    v.z *= scalar;
    return v;
}

void setL(Vector3D temp){
    l.x = temp.x;
    l.y = temp.y;
    l.z = temp.z;
}

void setR(Vector3D temp){
    r.x = temp.x;
    r.y = temp.y;
    r.z = temp.z;
}

void setU(Vector3D temp){
    u.x = temp.x;
    u.y = temp.y;
    u.z = temp.z;
}

void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}




void drawSquare(double a)
{
  glColor3f(1, 1, 1);
	glBegin(GL_QUADS);{
		glVertex3f( 0, a/2,a/2);
		glVertex3f( 0,-a/2,a/2);
		glVertex3f(0,-a/2,-a/2);
		glVertex3f(0, a/2,-a/2);
	}glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawUpperSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=h;
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=r*cos(((double)j/(double)slices)*2*pi);
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
		    if(j % 2 ==0){
                glColor3f(1, 1, 1);
		    }
		    else glColor3f(0, 0, 0);
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

			}glEnd();
		}
	}
}

void drawLowerSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=h;
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=r*cos(((double)j/(double)slices)*2*pi);
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
		    if(j % 2 ==0){
                glColor3f(1, 1, 1);
		    }
            else glColor3f(0, 0, 0);
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(-points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(-points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(-points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(-points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

			}glEnd();
		}
	}
}


void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawCylinder(double r, double H, int slices, int stacks){

    struct point points[100][100];
	int i,j;
	double h;
	//generate points
	for(i=0;i<=stacks;i++)
	{
	    h = ((double)H / (double)stacks) * i;
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=h;
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=r*cos(((double)j/(double)slices)*2*pi);
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
		    if(j % 2 ==0){
                glColor3f(1, 1, 1);
		    }
            else glColor3f(0, 0, 0);
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

			}glEnd();
		}
	}

}

void drawLastPart(double radius, int slices, int stacks){
    struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=2 * radius - radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=h;
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=r*cos(((double)j/(double)slices)*2*pi);
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
		    if(j % 2 ==0){
                glColor3f(1, 1, 1);
		    }
            else glColor3f(0, 0, 0);
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

			}glEnd();
		}
	}
}

void drawPoints(double x, double y, double z){
    glColor3f(1, 0, 0);
    glBegin(GL_POINTS);
    glVertex3f(x, y, z);
    glEnd();
}


void drawSS()
{
    /* glColor3f(1,1,1);
    glPushMatrix();
    glRotated(angleForQW, 0, 1, 0);
    drawLowerSphere(30, 30, 20);
    glRotated(angle, 0, 0, 1);
    drawUpperSphere(30, 30, 20);
    glTranslated(30, 0, 0);
    glRotated(angleForAS, 0, 0, 1);
    glRotated(angleForDF, 1, 0, 0);
    glTranslated(15, 0, 0);
    drawLowerSphere(15, 20, 15);
    drawCylinder(15, 60, 30, 30);
    glTranslated(60, 0, 0);
    drawLastPart(15, 30, 30);
    glTranslated(80, 0, 0);
    glGetDoublev(GL_MODELVIEW_MATRIX, p2);
    glPopMatrix();
    glPushMatrix();
    glTranslated(250, 0, 0);
    drawSquare(150);
    glPopMatrix();
    //glTranslated(249, 0, 0);
    //smallDrawSquare(5);
    //smallDrawSquare(20);
    for(int i = 0; i < pointlist.size(); i++){
        glPushMatrix();
        glTranslated(pointlist[i].x, pointlist[i].y, pointlist[i].z);
        smallDrawSquare(5);
        glPopMatrix();
    } */
    for(Object *obj: objects){
      obj->draw();
    }

}


void RotateAfterKeyListening(unsigned char key, double angle){
    Vector3D crossResult1, crossResult2, temp1, temp2, sum;
    if(key == '1' || key == '2'){
        crossResult1 = crossProduct(u, l);
        temp1.x = l.x * cos(angle) + crossResult1.x * sin(angle);
        temp1.y = l.y * cos(angle) + crossResult1.y * sin(angle);
        temp1.z = l.z * cos(angle) + crossResult1.z * sin(angle);
        setL(temp1);

        crossResult2 = crossProduct(u, r);
        temp1.x = r.x * cos(angle) + crossResult2.x * sin(angle);
        temp1.y = r.y * cos(angle) + crossResult2.y * sin(angle);
        temp1.z = r.z * cos(angle) + crossResult2.z * sin(angle);
        setR(temp1);
    }
    else if(key == '3' || key == '4'){
        crossResult1 = crossProduct(r, l);
        temp1.x = l.x * cos(angle) + crossResult1.x * sin(angle);
        temp1.y = l.y * cos(angle) + crossResult1.y * sin(angle);
        temp1.z = l.z * cos(angle) + crossResult1.z * sin(angle);
        setL(temp1);
        crossResult2 = crossProduct(r, u);
        temp1.x = u.x * cos(angle) + crossResult2.x * sin(angle);
        temp1.y = u.y * cos(angle) + crossResult2.y * sin(angle);
        temp1.z = u.z * cos(angle) + crossResult2.z * sin(angle);
        setU(temp1);
    }
    else if(key == '5' || key == '6'){
        crossResult1 = crossProduct(l, r);
        temp1.x = r.x * cos(angle) + crossResult1.x * sin(angle);
        temp1.y = r.y * cos(angle) + crossResult1.y * sin(angle);
        temp1.z = r.z * cos(angle) + crossResult1.z * sin(angle);
        setR(temp1);
        crossResult2 = crossProduct(l, u);
        temp1.x = u.x * cos(angle) + crossResult2.x * sin(angle);
        temp1.y = u.y * cos(angle) + crossResult2.y * sin(angle);
        temp1.z = u.z * cos(angle) + crossResult2.z * sin(angle);
        setU(temp1);
    }
}


void capture(){
	bitmap_image image(imageHeight,	imageWidth);
	planeDistance = (windowHeight / 2)/tan(fovY / 2.0 * (pi / 180.0));
	Vector3D topleft = pos + (l * planeDistance) - (r * (windowWidth/2)) + (u * (windowHeight/2));
	double du = (windowWidth * 1.0) / (imageWidth * 1.0);
	double dv = (windowHeight * 1.0) / (imageHeight * 1.0);

	// Choose middle of the grid cell
	topleft = topleft + r * (0.5*du) - u * (0.5*dv);
	topleft.printVector();

	int nearest;
	for(int i = 0; i < imageHeight; i++){
		for(int j = 0; j < imageWidth; j++){
		/* 	Vector3D second = ;
			Vector3D first = topleft + r * ((i + 0.5)*du); */
			Vector3D curPixel = topleft + r * ((i + 0.5)*du) - u * ((j + 0.5)*dv);
			Vector3D eyeTocurPixel = pos.p2pVector(curPixel);
			eyeTocurPixel.normalize();
			Color *color = new Color(0, 0, 0);
			Ray *ray = new Ray(pos, eyeTocurPixel);
			Object *keepTrack = nullptr;
			int nearest = -1;
			double t_min = 1e10;
			for(Object *obj : objects){
				double tVal = obj->intersect(ray, color,  0);
				if(tVal > 0){
					if(tVal < t_min){
						t_min = tVal;
						keepTrack = obj;
						nearest = 1;
					}
				}
			}
			if(nearest == -1){
				continue;
			}
			keepTrack->intersect(ray, color, 1);
			color->clip();
			image.set_pixel(i, j, color->r * 255, color->g * 255, color->b * 255);
			delete ray;
			delete color;
		}
	}
	for(int i = 0; i < objects.size(); i++){
		delete objects[i];
	}
	lights.clear();
	image.save_image("G://RayTracing//test.bmp");
	cout<<"done"<<endl;
/* 	double t, tMin;
	for i=1:imageWidth
		for j=1:imageHeight
			calculate curPixel using topleft,r,u,i,j,du,dv
			cast ray from eye to (curPixel-eye) direction
			double *color = new double[3]
			for each object, o in objects
			t = o.intersect(ray, dummyColor, 0)
			update t so that it stores min +ve value
			save the nearest object, on
			tmin = on->intersect(ray, color, 1)
			update image pixel (i,j)
	save image */
}

void keyboardListener(unsigned char key, int x,int y){
//    if(key == '1')  RotateAfterKeyListening(key, pi / 180);
//    else if(key == '2') RotateAfterKeyListening(key, -pi / 180);
    if(key == '1'){
        //Rotate(&l, &r, &u, pi/180);
        RotateAfterKeyListening(key, pi / 180);
    }
    else if(key == '2'){
        RotateAfterKeyListening(key, -pi / 180);
    }
    else if(key == '3'){
        RotateAfterKeyListening(key, pi / 180);
    }
    else if(key == '4'){
        RotateAfterKeyListening(key, -pi / 180);
    }
    else if(key == '5'){
        RotateAfterKeyListening(key, -pi / 180);
    }
    else if(key == '6'){
        RotateAfterKeyListening(key, pi / 180);
    }
		else if(key == '0'){
			cout<<"yo"<<endl;
			capture();
		}

    switch(key){

		case 'c':
			capture();
			break;

		case 'e':
			if(angle <= 70){
                angle += 5;
			}
			break;
        case 'r':
            if(angle >= -70){
                angle -= 5;
            }
            break;
        case 'q':
            if(angleForQW <= 60){
                angleForQW += 5;
            }
            break;
        case 'w':
            if(angleForQW >= -60){
                angleForQW -= 5;
            }
            break;
        case 'a':
            if(angleForAS <= 70){
                angleForAS += 5;
            }
            break;
        case 's':
            if(angleForAS >= -70){
                angleForAS -= 5;
            }
            break;
        case 'd':
            angleForDF += 5;
            break;
        case 'f':
            angleForDF -= 5;
            break;
		default:
			break;
	}

}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			forwardAlongLaxis(-5);
			break;
		case GLUT_KEY_UP:		// up arrow key
			forwardAlongLaxis(5);
			break;

		case GLUT_KEY_RIGHT:
			forwardAlongRaxis(5);
			break;
		case GLUT_KEY_LEFT:
			forwardAlongRaxis(-5);
			break;

		case GLUT_KEY_PAGE_UP:
		    forwardAlongUaxis(5);
			break;
		case GLUT_KEY_PAGE_DOWN:
		    forwardAlongUaxis(-5);
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}

void multiplyMatrix(int x){
    if(x == 1){
        for(int i = 0;i < 16; i++){
            int div = i / 4;
            int rem = i % 4;
            twod[div][rem] = p1[i];
        }
//        for(int i = 0; i < 4; i++){
//            for(int j = 0; j < 4; j++){
//                cout<<twod[i][j]<<" ";
//            }
//            cout<<endl;
//        }
        double oned[4] = {45, 0, 0, 1};

        for(int i = 0; i < 4; i++){
            double d = 0;
            for(int j = 0; j < 4; j++){
                d += twod[i][j] * oned[j];
            }
            if(i == 0) p1prime.x = d;
            else if(i == 1) p1prime.y = d;
            else if(i == 2) p1prime.z = d;
        }
    }
    else{
        for(int i = 0;i < 16; i++){
            int div = i / 4;
            int rem = i % 4;
            twod[div][rem] = p2[i];
        }

        double oned[4] = {60, 0, 0, 1};

        for(int i = 0; i < 4; i++){
            double d = 0;
            for(int j = 0; j < 4; j++){
                d += twod[i][j] * oned[j];
            }
            if(i == 0) p2prime.x = d;
            else if(i == 1) p2prime.y = d;
            else if(i == 2) p2prime.z = d;
        }
    }
}

void printpoint(point p){
    cout<<p.x <<" "<<p.y <<" "<<p.z<<endl;
}

void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_RIGHT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_LEFT_BUTTON:
			if(state == GLUT_UP) {

                double r = 30 * cos(pi * angle / 180);
                p1prime.y = 30 * sin(pi * angle / 180);
                p1prime.x = r * cos(pi * -angleForQW / 180);
                p1prime.z = r * sin(pi * -angleForQW / 180);

                Vector3D vect;
                r = 30 * cos(pi * (angle + angleForAS) / 180);
                vect.y = 30 * sin(pi * (angle + angleForAS) / 180);
                vect.x = r * cos(pi * -angleForQW / 180);
                vect.z = r * sin(pi * -angleForQW / 180);

                point result;
                printpoint(p1prime);
                double t = (249 - p1prime.x) / vect.x;
                result.x = 249;
                result.y = p1prime.y + t * vect.y;
                result.z = p1prime.z + t * vect.z;

                printpoint(result);
                if(abs(result.y) < 75 and abs(result.z) < 75) pointlist.push_back(result);

			}
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

  drawSS();

    //drawCircle(30,24);

    //drawCone(20,50,24);

	//drawSphere(30,24,20);




	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	//angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void loadTestData(){


	Object *temp;
	ifstream MyReadFile1("G://RayTracing//scene.txt");
	MyReadFile1>>recursion>>pixel;

	int numOfObjects;
	MyReadFile1>>numOfObjects;

	for(int i = 0; i < numOfObjects; i++){
		string nameOfObjects;
		MyReadFile1>>nameOfObjects;
		if(nameOfObjects == "sphere"){
			double centerX, centerY, centerZ, radius, red, green, blue, ambient, diffuse, specular, recursive;
			int shininess;
			MyReadFile1>>centerX>>centerY>>centerZ>>radius>>red>>green>>blue>>ambient>>diffuse>>specular>>recursive>>shininess;
			temp=new Sphere(Vector3D(centerX, centerY, centerZ), radius);
			temp->setColor(red, green, blue);
			temp->setCoEfficients(ambient, diffuse, specular, recursive);
			temp->setShine(shininess);
			objects.push_back(temp);
		}
		else if(nameOfObjects == "triangle"){
			double x1, y1, z1, x2, y2, z2, x3, y3, z3, red, green, blue, ambient, diffuse, specular, recursive;
			int shininess;
			MyReadFile1>>x1>>y1>>z1>>x2>>y2>>z2>>x3>>y3>>z3>>red>>green>>blue>>ambient>>diffuse>>specular>>recursive>>shininess;
			temp = new Triangle(Vector3D(x1, y1, z1), Vector3D(x2, y2, z2), Vector3D(x3, y3, z3));
			temp->setColor(red, green, blue);
			temp->setCoEfficients(ambient, diffuse, specular, recursive);
			temp->setShine(shininess);
			objects.push_back(temp);
		}
		else if(nameOfObjects == "general"){
			double a, b, c, d, e, f, g, h, i, j, rx, ry, rz, length, width, height, red, green, blue, ambient, diffuse, specular, recursive;
			int shininess;
			MyReadFile1>>a>>b>>c>>d>>e>>f>>g>>h>>i>>j>>rx>>ry>>rz>>length>>width>>height>>red>>green>>blue>>ambient>>diffuse>>specular>>recursive>>shininess;
			temp = new Quadric(Vector3D(rx, ry, rz), length, height, width);
			temp->setColor(red, green, blue);
			temp->setCoEfficients(ambient, diffuse, specular, recursive);
			temp->setShine(shininess);
			temp->setConstants(a, b, c, d, e, f, g, h, i, j);
			objects.push_back(temp);
		}

	}

	int numOfLights;
	MyReadFile1>>numOfLights;

	for(int i = 0; i < numOfLights; i++){
		double centerX, centerY, centerZ, red, green, blue;
		MyReadFile1>>centerX>>centerY>>centerZ>>red>>green>>blue;
		lights.push_back(Light(Vector3D(centerX, centerY, centerZ), red, green, blue));
	}

	temp = new Floor(1000, 20);
	temp->setCoEfficients(.4, .3, .2, .1);
	temp->setShine(3);
	objects.push_back(temp);
	cout<<"object length "<<objects.size()<<endl;

}

void init(){

	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;
  angleForQW = 0;
	angleForAS = 0;
	angleForDF = 0;
	fovY = 100, aspectRatio = 1;
	windowHeight = 500, windowWidth = 500;
	imageHeight = 1024, imageWidth = 1024;
	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(fovY,	aspectRatio, 1, 	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance

	pos.x = 100, pos.y = 0, pos.z = 100;
	u.x = 0, u.y = 0, u.z = 1;
	l.x = -1/sqrt(2), l.y = -1/sqrt(2), l.z = 0;
	r.x = -1/sqrt(2), r.y = 1/sqrt(2), r.z = 0;

	loadTestData();
}


int main(int argc, char **argv){

	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
