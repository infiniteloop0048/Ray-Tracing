#include<iostream>
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
#include<bits/stdc++.h>
#include <windows.h>
#include <GL/glut.h>
#include "bitmap_image.hpp"

using namespace std;
#define pi (2*acos(0.0))

int imageHeight, imageWidth, windowHeight, windowWidth;
double fovY, fovX, aspectRatio, planeDistance;
int recursion, pixel;

class Vector3D {
  public:
    double x, y, z;
    Vector3D(){

    }
    Vector3D(double a, double b, double c){
      x = a;
      y = b;
      z = c;
    }
    void normalize(){
      double d = sqrt((x * x) + (y * y) + (z * z));
      x = x / d;
      y = y / d;
      z = z / d;
    }
    Vector3D p2pVector(Vector3D point){
      Vector3D temp;
      temp.x = point.x - x;
      temp.y = point.y - y;
      temp.z = point.z - z;
      return temp;
    }
    Vector3D operator*(double val) {
        Vector3D temp;
        temp.x = x * val;
        temp.y = y * val;
        temp.z = z * val;
        return temp;
    }
    Vector3D operator+(Vector3D vect) {
        Vector3D temp;
        temp.x = vect.x + x;
        temp.y = vect.y + y;
        temp.z = vect.z + z;
        return temp;
    }
    Vector3D operator-(Vector3D vect) {
        Vector3D temp;
        temp.x = vect.x - x;
        temp.y = vect.y - y;
        temp.z = vect.z - z;
        return temp * -1;
    }

    double dist(Vector3D v){
      return sqrt(pow(x - v.x, 2) + pow(y - v.y, 2) + pow(z - v.z, 2));
    }

    double dot(Vector3D vect){
      return (x * vect.x) + (y * vect.y) + (z * vect.z);
    }

    Vector3D crossProduct(Vector3D v2){
      Vector3D cross_p;
      cross_p.x = (y*v2.z) - (z*v2.y);
      cross_p.y = -((x*v2.z) - (z*v2.x));
      cross_p.z = (x*v2.y) - (y*v2.x);
      return cross_p;
    }
    void printVector(){
      cout<<x<<" "<<y<<" "<<z<<endl;
    }

};

Vector3D l, r, u, pos;

class Color {
public:
    double r, g, b;

    Color() {
        r = 0;
        g = 0;
        b = 0;
    }

    Color(double R, double G, double B) {
        r = R;
        g = G;
        b = B;
    }

    void setRGB(double R, double G, double B) {
        r = r;
        g = g;
        b = b;
    }

    void print() {
        cout << "r " << (int) r << " g " << (int) g << " b " << (int) b << endl;
    }

    void clip(){
      r = min(r, 1.0);
      g = min(g, 1.0);
      b = min(b, 1.0);
    }
};


class Ray {
public:
    Vector3D starPoint;
    Vector3D dir;
    double tVal;
    Color color;

    Ray() {
        tVal = -1;
    }

    Ray(Vector3D start,  Vector3D direction) {
        starPoint = start;
        dir = direction;
        tVal = -1;
    }
};


struct point
{
	double x,y,z;
};

class Light{
  public:
    Vector3D center;
    double color[3];
    void setColor(double r, double g, double b){
      color[0] = r;
      color[1] = g;
      color[2] = b;
    }
    void setCenter(Vector3D c){
      center = c;
    }
    Light(Vector3D c, double r, double g, double b){
      center = c;
      color[0] = r, color[1] = g, color[2] = b;
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
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
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
    void draw(){
      glColor3f(color[0], color[1], color[2]);
      glPushMatrix();
			glTranslated(center.x, center.y, center.z);
			drawSphere(5, 30, 20);
			glPopMatrix();
    }
};

class Object {
  public:
    Vector3D reference_point;
    double height, width, length;
    int Shine;
    double color[3];
    double co_efficients[4];

    Object(){

    }

    virtual void draw(){

    }
    virtual double getParamVal(Ray *ray){

    }
    void setColor(double c1, double c2, double c3){
      color[0] = c1, color[1] = c2, color[2] = c3;
    }
    void setShine(int shine){
      Shine = shine;
    }
    void setCoEfficients(double c1, double c2, double c3, double c4){
      co_efficients[0] = c1;
      co_efficients[1] = c2;
      co_efficients[2] = c3;
      co_efficients[3] = c4;
    }
    virtual double intersect(Ray *r, Color *color, int level){
      return -1.0;
    }
    virtual void setConstants(double a, double b, double c,
                      double d, double e, double f,
                      double g, double h, double i, double j){

    }
};

vector <Object *> objects;
vector <Light> lights;

class Sphere : public Object{
  public:
    Sphere(Vector3D point3D, double radius){
      reference_point = point3D;
      length = radius;
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
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
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

    void draw(){
      glColor3f(color[0], color[1], color[2]);
      glPushMatrix();
			glTranslated(reference_point.x, reference_point.y, reference_point.z);
			drawSphere(length, 30, 20);
			glPopMatrix();
    }
    double getParamVal(Ray *ray){
      double a = ray->dir.dot(ray->dir);
      double b = 2 * ray->dir.dot(ray->starPoint - reference_point);
      double c = (ray->starPoint - reference_point).dot(ray->starPoint - reference_point) - (length * length);
      double discriminant = b * b - (4 * a * c);
      if(discriminant < 0) return -1;
      else if(discriminant == 0){
        return -b/(2 * a);
      }
      else{
        double t1 = (-b + sqrt(discriminant)) / (2 * a);
        double t2 = (-b - sqrt(discriminant)) / (2 * a);
        if(t1 >= 0 && t2 >= 0) return min(t1, t2);
        if(t1 < 0 && t2 < 0) return -1;
        if(t1 < 0) return t2;
        if(t2 < 0) return t1;      }
    }
    double intersect(Ray *ray, Color *col, int level){
      double t = getParamVal(ray);
      if(level == 0) return t;

      col->r = color[0] * co_efficients[0], col->g = color[1] * co_efficients[0], col->b = color[2] * co_efficients[0];

      Vector3D ip = ray->starPoint + (ray->dir * t);
      Vector3D N = reference_point.p2pVector(ip);
      N.normalize();
      Vector3D V = ray->dir * -1;
      //light effect
      for(Light light: lights){

        Vector3D L = ip.p2pVector(light.center);
        double L_length = sqrt(L.dot(L));
        L.normalize();
        Vector3D point = ip + L;
        Ray *newRay = new Ray(point, L);

        bool obstacle = false;

        for(Object *obj: objects){
          double val = obj->getParamVal(newRay);
          if(val < 0 || val > L_length){
            continue;
          }
          obstacle = true;
          break;
        }

        if(!obstacle){
          Vector3D R = N * (2 * (N.dot(L))) - L;
          col->r += light.color[0] * ((co_efficients[1] * color[0] * max(L.dot(N), 0.0)) + (co_efficients[2] * max(pow(R.dot(V), Shine), 0.0)));
          col->g += light.color[1] * ((co_efficients[1] * color[1] * max(L.dot(N), 0.0)) + (co_efficients[2] * max(pow(R.dot(V), Shine), 0.0)));
          col->b += light.color[2] * ((co_efficients[1] * color[2] * max(L.dot(N), 0.0)) + (co_efficients[2] * max(pow(R.dot(V), Shine), 0.0)));
        }
        delete newRay;
      
      }

      if(level < recursion){
        Vector3D reflected_ray = N * (2 * N.dot(V)) - V;
        reflected_ray.normalize();
        Vector3D point = ip + reflected_ray;
        Ray *refRay = new Ray(point, reflected_ray);

        double minimum = INT_MAX;
        Color *reflectedColor = new Color();

        Object *track = nullptr;
        int nearer = -1;

        for(Object *obj: objects){
          double val = obj->intersect(refRay, reflectedColor, 0);
          if(val > 0)
          if(val < minimum) {
            minimum = val;
            track = obj;
            nearer = 1;
          }
        }

        if(nearer != -1){
          track->intersect(refRay, reflectedColor, level + 1);
          col->r += reflectedColor->r * co_efficients[3];
          col->g += reflectedColor->g * co_efficients[3];
          col->b += reflectedColor->b * co_efficients[3];
        }
        delete refRay;
        delete reflectedColor;
      }

      /* col->r = (col->r)/max(col->r, max(col->g, col->b));
      col->g = (col->g)/max(col->r, max(col->g, col->b));
      col->b = (col->b)/max(col->r, max(col->g, col->b)); */
    }
};

class Triangle: public Object{
  public:
    Vector3D point[3];
    Triangle(Vector3D p1, Vector3D p2, Vector3D p3){
      point[0] = p1, point[1] = p2, point[2] = p3;
    }
    void drawTriangle(){
      glColor3f(color[0], color[1], color[2]);
      glBegin(GL_TRIANGLES);{
        glVertex3f( point[0].x, point[0].y, point[0].z);
        glVertex3f( point[1].x, point[1].y, point[1].z);
        glVertex3f( point[2].x, point[2].y, point[2].z);
      }glEnd();
    }
    void draw(){
      glColor3f(color[0], color[1], color[2]);
      glPushMatrix();
      drawTriangle();
      glPopMatrix();
    }


    double getParamVal(Ray *ray){
      const float EPSILON = 0.0000001;
      Vector3D vertex0 = point[0];
      Vector3D vertex1 = point[1];
      Vector3D vertex2 = point[2];
      Vector3D edge1, edge2, h, s, q;
      float a,f,u,v;
      edge1 = vertex1 - vertex0;
      edge2 = vertex2 - vertex0;
      h = ray->dir.crossProduct(edge2);
      a = (float)edge1.dot(h);
      if (a > -EPSILON && a < EPSILON)
          return -1;    // This ray is parallel to this triangle.
      f = 1.0/a;
      s = ray->starPoint - vertex0;
      u = f * s.dot(h);
      if (u < 0.0 || u > 1.0)
          return -1;
      q = s.crossProduct(edge1);
      v = (float)f * ray->dir.dot(q);
      if (v < 0.0 || u + v > 1.0)
          return -1;
      // At this stage we can compute t to find out where the intersection point is on the line.
      float t = f * edge2.dot(q);
      if (t > EPSILON) // ray intersection
      {
          return (double)t;
      }
      else // This means that there is a line intersection but not a ray intersection.
          return -1;
    }

    double intersect(Ray *ray, Color *col, int level){
      double tVal = getParamVal(ray);
      if(level == 0) return tVal;
      col->r = color[0] * co_efficients[0], col->g = color[1] * co_efficients[0], col->b = color[2] * co_efficients[0];

      Vector3D ip = ray->starPoint + (ray->dir * tVal);
      Vector3D N = (point[1] - point[0]).crossProduct(point[2] - point[0]);
      N.normalize();
      Vector3D V = ray->dir * -1;
      //light effect
      for(Light light: lights){

        Vector3D L = ip.p2pVector(light.center);
        double L_length = sqrt(L.dot(L));
        L.normalize();
        Vector3D point = ip + L;
        Ray *newRay = new Ray(point, L);

        bool obstacle = false;

        for(Object *obj: objects){
          double val = obj->getParamVal(newRay);
          if(val < 0 || val > L_length){
            continue;
          }
          obstacle = true;
          break;
        }
        if(!obstacle){
          Vector3D R = N * (2 * (N.dot(L))) - L;
          col->r += light.color[0] * ((co_efficients[1] * color[0] * max(L.dot(N), 0.0)) + (co_efficients[2] * max(pow(R.dot(V), Shine), 0.0)));
          col->g += light.color[1] * ((co_efficients[1] * color[1] * max(L.dot(N), 0.0)) + (co_efficients[2] * max(pow(R.dot(V), Shine), 0.0)));
          col->b += light.color[2] * ((co_efficients[1] * color[2] * max(L.dot(N), 0.0)) + (co_efficients[2] * max(pow(R.dot(V), Shine), 0.0)));
        }
        delete newRay;
      }
      if(level < recursion){
        Vector3D reflected_ray = N * (2 * N.dot(V)) - V;
        reflected_ray.normalize();
        Vector3D point = ip + reflected_ray;
        Ray *refRay = new Ray(point, reflected_ray);

        double minimum = INT_MAX;
        Color *reflectedColor = new Color();

        Object *track = nullptr;
        int nearer = -1;

        for(Object *obj: objects){
          double val = obj->intersect(refRay, reflectedColor, 0);
          if(val > 0)
          if(val < minimum) {
            minimum = val;
            track = obj;
            nearer = 1;
          }
        }

        if(nearer != -1){
          track->intersect(refRay, reflectedColor, level + 1);
          col->r += reflectedColor->r * co_efficients[3];
          col->g += reflectedColor->g * co_efficients[3];
          col->b += reflectedColor->b * co_efficients[3];
        }
        delete refRay;
        delete reflectedColor;
      }
    }

};


class Quadric: public Object{
  public:
    double A, B, C, D, E, F, G, H, I, J;
    Quadric(Vector3D point, double l, double h, double w){
      reference_point = point;
      height = h;
      width = w;
      length = l;
    }
    void setConstants(double a, double b, double c,
                      double d, double e, double f,
                      double g, double h, double i, double j){
      A = a, B = b, C = c, D = d, E = e,
      F = f, G = g, H = h, I = i, J = j;
    }
    void draw(){
      return ;
    }
    double getParamVal(Ray *ray){
      double x0 = ray->starPoint.x;
      double xd = ray->dir.x;
      double y0 = ray->starPoint.y;
      double yd = ray->dir.y;
      double z0 = ray->starPoint.z;
      double zd = ray->dir.z;

      double Aq = (A * xd * xd) + (B * yd * yd) + (C * zd * zd) +
      (D * xd * yd) + (E * xd * zd) + (F * yd * zd);

      double Bq = (2 * A * x0 * xd) + (2 * B * y0 * yd) +
          (2 * C * z0 * zd) + D * ((x0 * yd) + (y0 * xd)) +
          E * ((x0 * zd) + (z0 * xd)) +
          F * ((y0 * zd) + (yd * z0)) +
          (G * xd) + (H * yd) + (I * zd);
      double Cq = (A * x0 * x0) + (B * y0 * y0) + (C * z0 * z0) + (D * x0 * y0) + (E * x0 * z0) + (F * y0 * z0) + (G * x0) + (H * y0) + (I * z0) + J;
      double t0, t1;
      if(Aq == 0){
        return -Cq / Bq;
      }else{
        double det = Bq * Bq - 4 * Aq * Cq;
        if(det < 0){
          return -1;
        }
        else{

          bool flag1 = true, flag2 = true;
          t0 =( - Bq - pow(det, .5)) / (2 * Aq);
          t1 =( - Bq + pow(det, .5)) / (2 * Aq);
          if(t0 > 0) {
            Vector3D ip = ray->starPoint + (ray->dir * t0);

            if(length > 0){
              if(ip.x > (reference_point.x + length) || ip.x < reference_point.x) flag1 = false;
            }
            if(width > 0){
              if(ip.y > (reference_point.y + width) || ip.y < reference_point.y) flag1 = false;
            }
            if(height > 0){
              if(ip.z > (reference_point.z + height) || ip.z < reference_point.z) flag1 = false;
            }
          }
          Vector3D ip = ray->starPoint + (ray->dir * t1);
          if(length > 0){
            if(ip.x > (reference_point.x + length) || ip.x < reference_point.x) flag2 = false;
          }
          if(width > 0){
            if(ip.y > (reference_point.y + width) || ip.y < reference_point.y) flag2 = false;
          }
          if(height > 0){
            if(ip.z > (reference_point.z + height) || ip.z < reference_point.z) flag2 = false;
          }
          if(flag1 && flag2) return min(t0, t1);
          else if(flag1) return t0;
          else if(flag2) return t1;
          else return -1;
          /* if(t0 < 0 && t1 < 0) return -1;
          if(t0 >= 0 && t1 >= 0) return min(t0, t1);
          if(t0 >= 0) return t0;
          if(t1 >= 0) return t1; */
        }
      }
    }
    double intersect(Ray *ray, Color *col, int level){
      double tVal = getParamVal(ray);
      if(level == 0) return tVal;
      col->r = color[0] * co_efficients[0], col->g = color[1] * co_efficients[0], col->b = color[2] * co_efficients[0];

      Vector3D ip = ray->starPoint + (ray->dir * tVal);

      //normal code
      Vector3D N;
      N.x = 2*A*ip.x + D*ip.y + E*ip.z + G;
      N.y = 2*B*ip.y + D*ip.x + F*ip.z + H;
      N.z = 2*C*ip.z + E*ip.x + F*ip.y + I;
      N.normalize();

      Vector3D V = ray->dir * -1;
      //light effect
      for(Light light: lights){

        Vector3D L = ip.p2pVector(light.center);
        double L_length = sqrt(L.dot(L));
        L.normalize();
        Vector3D point = ip + L;
        Ray *newRay = new Ray(point, L);

        bool obstacle = false;

        for(Object *obj: objects){
          double val = obj->getParamVal(newRay);
          if(val < 0 || val > L_length){
            continue;
          }
          obstacle = true;
          break;
        }
        if(!obstacle){
          Vector3D R = N * (2 * (N.dot(L))) - L;
          col->r += light.color[0] * ((co_efficients[1] * color[0] * max(L.dot(N), 0.0)) + (co_efficients[2] * max(pow(R.dot(V), Shine), 0.0)));
          col->g += light.color[1] * ((co_efficients[1] * color[1] * max(L.dot(N), 0.0)) + (co_efficients[2] * max(pow(R.dot(V), Shine), 0.0)));
          col->b += light.color[2] * ((co_efficients[1] * color[2] * max(L.dot(N), 0.0)) + (co_efficients[2] * max(pow(R.dot(V), Shine), 0.0)));
        }
        delete newRay;
      }
      if(level < recursion){
        Vector3D reflected_ray = N * (2 * N.dot(V)) - V;
        reflected_ray.normalize();
        Vector3D point = ip + reflected_ray;
        Ray *refRay = new Ray(point, reflected_ray);

        double minimum = INT_MAX;
        int nearer = -1;
        Color *reflectedColor = new Color();

        Object *track = nullptr;
        for(Object *obj: objects){
          double val = obj->intersect(refRay, reflectedColor, 0);
          if(val > 0)
          if(val < minimum) {
            minimum = val;
            track = obj;
            nearer = 1;
          }
        }

        if(nearer != -1){
          track->intersect(refRay, reflectedColor, level + 1);
          col->r += reflectedColor->r * co_efficients[3];
          col->g += reflectedColor->g * co_efficients[3];
          col->b += reflectedColor->b * co_efficients[3];
        }
        delete refRay;
        delete reflectedColor;
      }
    }

};

class Floor: public Object{
  public:
    int iterations, floorWidth;
    Floor(double Floorwidth, double Tilewidth){
      floorWidth = Floorwidth;
      reference_point = Vector3D(-Floorwidth/2, -Floorwidth/2, 0);
      length = Tilewidth;
      iterations = Floorwidth / Tilewidth;
    }

    double getParamVal(Ray *ray){
      double a = ray->starPoint.z;
      double b = ray->dir.z;
      if(b == 0) return -1;
      double tVal = -a/b;
      return tVal;
    }

    double intersect(Ray *ray, Color *col, int level){

      double tVal = getParamVal(ray);

      if(level == 0)
        return tVal;

      double x = (ray->starPoint.x) + (tVal * ray->dir.x);
      double y = (ray->starPoint.y) + (tVal * ray->dir.y);
      double rHere, gHere, bHere;
      if(x >= -floorWidth/2 && x <= floorWidth/2 && y >= -floorWidth/2 && y <= floorWidth/2){
        int floorX = round((x + floorWidth/2)/length);
        int floorY = round((y + floorWidth/2)/length);
        if(floorX % 2 == 0){
          if(floorY % 2 == 1){
            col->r = 0, col->g = 0, col->b = 0;
            color[0] = 0, color[1] = 0, color[2] = 0;
          }
          else {
            col->r = co_efficients[0], col->g = co_efficients[0], col->b = co_efficients[0];
            color[0] = 1, color[1] = 1, color[2] = 1;
          }
        }
        else{
          if(floorY % 2 == 1){
            col->r = co_efficients[0], col->g = co_efficients[0], col->b = co_efficients[0];
            color[0] = 1, color[1] = 1, color[2] = 1;
          }
          else {
            col->r = 0, col->g = 0, col->b = 0;
            color[0] = 0, color[1] = 0, color[2] = 0;
          }
        }
      }
      //col->r = color[0] * co_efficients[0], col->g = color[1] * co_efficients[0], col->b = color[2] * co_efficients[0];

      Vector3D ip = ray->starPoint + (ray->dir * tVal);
      Vector3D N(0, 0, 1);
      if(pos.z < 0) N =  N * -1;

      //N.normalize();
      Vector3D V = ray->dir * -1;
      //light effect

      for(Light light: lights){

        Vector3D L = ip.p2pVector(light.center);
        //L.normalize();
        double L_length = sqrt(L.dot(L));
        L.normalize();
        Vector3D point = ip + L;
        Ray *newRay = new Ray(point, L);

        bool obstacle = false;

        for(Object *obj: objects){
          double val = obj->getParamVal(newRay);
          if(val < 0 || val > L_length){
            continue;
          }
          obstacle = true;
          break;
        }

        if(!obstacle){
          Vector3D R = N * (2 * (N.dot(L))) - L;
          col->r += light.color[0] * ((co_efficients[1] * color[0] * max(L.dot(N), 0.0)) + (co_efficients[2] * max(pow(R.dot(V), Shine), 0.0)));
          col->g += light.color[1] * ((co_efficients[1] * color[1] * max(L.dot(N), 0.0)) + (co_efficients[2] * max(pow(R.dot(V), Shine), 0.0)));
          col->b += light.color[2] * ((co_efficients[1] * color[2] * max(L.dot(N), 0.0)) + (co_efficients[2] * max(pow(R.dot(V), Shine), 0.0)));
        }
        delete newRay;

      }
      
      if(level < recursion){
        Vector3D reflected_ray = N * (2 * N.dot(V)) - V;
        reflected_ray.normalize();
        Vector3D point = ip + reflected_ray;
        Ray *refRay = new Ray(point, reflected_ray);

        double minimum = INT_MAX;
        Color *reflectedColor = new Color();

        Object *track = nullptr;
        int nearer = -1;

        for(Object *obj: objects){
          double val = obj->intersect(refRay, reflectedColor, 0);
          if(val > 0)
          if(val < minimum) {
            minimum = val;
            track = obj;
            nearer = 1;
          }
        }

        if(nearer != -1){
          track->intersect(refRay, reflectedColor, level + 1);
          col->r += reflectedColor->r * co_efficients[3];
          col->g += reflectedColor->g * co_efficients[3];
          col->b += reflectedColor->b * co_efficients[3];
        }
        delete refRay;
        delete reflectedColor;
      }
    }

    void draw(){
      double startX = reference_point.x;
      double startY = reference_point.y;
      glColor3f(1, 1, 1);
      bool first = true, second = true;
      for(int i = 0; i < iterations; i++){
        first = !first;
        startY = reference_point.y;
        for(int j = 0; j < iterations; j++){
          if(first){
            if(second){
              glColor3f(1, 1, 1);
            }
            else glColor3f(0, 0, 0);
          }
          else{
            if(second){
              glColor3f(0, 0, 0);
            }
            else glColor3f(1, 1, 1);
          }
          glBegin(GL_QUADS);{

            glVertex3f(startX, startY, 0);
            glVertex3f(startX + length, startY, 0);
            glVertex3f(startX + length, startY + length, 0);
            glVertex3f(startX, startY + length, 0);


          }glEnd();
          startY += length;
          second = !second;
        }
        startX += length;
      }
    }
};




