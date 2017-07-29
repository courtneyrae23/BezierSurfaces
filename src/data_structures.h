#ifndef __DATA_H_INCLUDED__   // if x.h hasn't been included yet...
#define __DATA_H_INCLUDED__   //   #define this so the compiler knows it has been included


#include <math.h>
#include <vector>


class Point {
  public:
    float x, y, z;
    Point(float a = 0.0, float b = 0.0, float c = 0.0) {
        x = a;
        y = b;
        z = c;
    };
    Point operator*(const float);
    Point operator+(const Point&);
    Point operator-(const Point&);
    void print();
};

class UV {
  public:
    float u, v;
    UV(float a = 0.0, float b = 0.0) {
        u = a;
        v = b;
    };
    void print();
};

class Curve {
  public:
    Point* p1;
    Point* p2;
    Point* p3;
    Point* p4;
    Curve(std::vector<Point*> points) {
      p1 = points[0];
      p2 = points[1];
      p3 = points[2];
      p4 = points[3];
    };
    Curve(std::vector<Point> points) {
      p1 = &points[0];
      p2 = &points[1];
      p3 = &points[2];
      p4 = &points[3];
    };
    Curve(Point* a, Point* b, Point* c, Point* d) {
      p1 = a;
      p2 = b;
      p3 = c;
      p4 = d;
    };
};

class Patch {
  public:
    Curve* c1;
    Curve* c2;
    Curve* c3;
    Curve* c4;
    Patch(std::vector<Curve*> curves) {
      c1 = curves[0];
      c2 = curves[1];
      c3 = curves[2];
      c4 = curves[3];
    };
};

class Normal {
  public:
    float x, y, z, len;
    Normal(float a = 0.0, float b = 0.0, float c = 0.0) {
        float length = sqrt(a*a + b*b + c*c);
        this->len = length;
        if (length == 0.0) {
            x = 0.0;
            y = 0.0;
            z = 0.0;  
        } else {
            x = a / length;
            y = b / length;
            z = c / length;
        }
    };
    void print();
};

class Triangle {
  public:
    Point *v1, *v2, *v3;
    Normal *n1, *n2, *n3;
    UV *uv1, *uv2, *uv3;
    Triangle(Point* a, Point* b, Point* c, Normal* a1, Normal* b1, Normal* c1, UV* a2, UV* b2, UV* c2) {
        v1 = a;
        v2 = b;
        v3 = c;
        n1 = a1;
        n2 = b1;
        n3 = c1;
        uv1 = a2;
        uv2 = b2;
        uv3 = c2;
    };
    Point* vert_midpoint1();
    Point* vert_midpoint2();
    Point* vert_midpoint3();
    UV* uv_midpoint1();
    UV* uv_midpoint2();
    UV* uv_midpoint3();
    void print();
};

#endif