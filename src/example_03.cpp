#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string.h>

//include header file for glfw library so that we can use OpenGL
#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "data_structures.h"

#ifdef _WIN32
static DWORD lastTime;
#else
static struct timeval lastTime;
#endif

#define PI 3.14159265 // Should be used from mathlib

using namespace std;

/*
For UC Berkeley's CS184 Fall 2016 course, assignment 3 (Bezier surfaces)
*/

//****************************************************
// Global Variables
//****************************************************
GLfloat translation[3] = {0.0f, 0.0f, 0.0f};
GLfloat scale[3] = {1.0f, 1.0f, 1.0f};
vector<vector<float> > rotations = *(new vector<vector<float> >());\
bool wiref = false;
bool doub = false;
bool hidd = false;
bool doubh = false;
bool auto_strech = false;
int Width_global = 400;
int Height_global = 400;
int Z_buffer_bit_depth = 128;
float numPatches;
float EPSILON = .001;
float numVertices = 0;
bool adaptive = false;
bool output = false;
string to_output = "";
bool smooth = true;
bool doubs = false;
std::vector<Patch*> patches;
std::vector<std::vector<Point*>> patch_v;
std::vector<std::vector<Normal*>> patch_n;
std::vector<Triangle*> triangles;
std::vector<std::vector<Point>> tris;
std::vector<std::vector<Normal>> tris_n;
bool obj = false;
bool normals_exist = false;
int numTriangles = 0;

inline float sqr(float x) { return x*x; }


//****************************************************
// Simple init function
//****************************************************
void initializeRendering()
{
    glfwInit();
}

//****************************************************
// A routine to set a pixel by drawing a GL point.  This is not a
// general purpose routine as it assumes a lot of stuff specific to
// this example.
//****************************************************
void setPixel(float x, float y, GLfloat r, GLfloat g, GLfloat b) {
    glColor3f(r, g, b);
    glVertex2f(x+0.5, y+0.5);  // The 0.5 is to target pixel centers
    // Note: Need to check for gap bug on inst machines.
}

//****************************************************
// Keyboard inputs. Add things to match the spec! 
//****************************************************
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    switch (key) {
        

        case GLFW_KEY_EQUAL:
          if (action && mods == GLFW_MOD_SHIFT){
            scale[0] += 0.01f;
            scale[1] += 0.01f;
            scale[2] += 0.01f;
          } break;
        case GLFW_KEY_MINUS:
          scale[0] -= 0.01f;
          scale[1] -= 0.01f;
          scale[2] -= 0.01f;
          break;
        case GLFW_KEY_ESCAPE: glfwSetWindowShouldClose(window, GLFW_TRUE); break;
        case GLFW_KEY_Q: glfwSetWindowShouldClose(window, GLFW_TRUE); break;
        case GLFW_KEY_LEFT :
          if (action && mods == GLFW_MOD_SHIFT){
            translation[0] -= 0.001f * Width_global;
            break;
          } else{
            vector<float> to_push = vector<float> {10, 0, 1, 0};
            rotations.push_back(to_push);
            break;
          }
        case GLFW_KEY_RIGHT:
          if (action && mods == GLFW_MOD_SHIFT){ 
            translation[0] += 0.001f * Width_global;
            break;
          }else{
            vector<float> to_push = vector<float> {-5,0,1,0};
            rotations.push_back(to_push);
            //vector<float> to_push2 = vector<float> {-5, 1, 0, 0};
            //rotations.push_back(to_push2);
            break;
          }
        case GLFW_KEY_UP:
          if (action && mods == GLFW_MOD_SHIFT){ 
            translation[1] += 0.001f * Height_global;
            break;
          } else{
            vector<float> to_push = vector<float> {10, 1, 0, 0};
            rotations.push_back(to_push);
            break;
          }
        case GLFW_KEY_DOWN :
          if (action && mods == GLFW_MOD_SHIFT){ 
            translation[1] -= 0.001f * Height_global;
            break;
          }
          else{
            vector<float> to_push = vector<float> {-5,1,0,0};
            rotations.push_back(to_push);
            break;
          }
        case GLFW_KEY_F:
          if (action && mods == GLFW_MOD_SHIFT) auto_strech = !auto_strech; break;
        case GLFW_KEY_S:
          if (smooth and doubs) {
            smooth = false;
            doubs = false;
            break;
          } else if (doubs) {
            smooth = true;
            doubs = false;
            break;
          } else {
            doubs = true;
            break;
          }
        case GLFW_KEY_W:
          hidd = false;
          if(not wiref and doub){
            wiref = true;
            doub = false;
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            break;
          } else if(doub) {
            wiref = false;
            doub = false;
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            break;
          } else {
            doub = true;
            break;
          }
        case GLFW_KEY_H:
          if(not hidd and doubh){
            hidd = true;
            doubh = false;
            break;
          } else if (doubh){
            hidd = false;
            doubh = false;
            break;
          } else {
            doubh = true;
            break;
          }
        case GLFW_KEY_SPACE: break;
            
        default: break;
    }
    
}

void drawSurface() {
  glBegin(GL_QUADS);

  float div = sqrt(numVertices/numPatches);
  int offset;

  for (int k = 0; k < numPatches; k++) {
    for (int i = 0; i < div-1; i++) {
      for (int j = 0; j < (div-1)*(div-1); j+=div) {

        glNormal3f(patch_n[k][j+i]->x, patch_n[k][j+i]->y, patch_n[k][j+i]->z);
        glVertex3f(patch_v[k][j+i]->x, patch_v[k][j+i]->y, patch_v[k][j+i]->z);

        glNormal3f(patch_n[k][j+i+1]->x, patch_n[k][j+i+1]->y, patch_n[k][j+i+1]->z);
        glVertex3f(patch_v[k][j+i+1]->x, patch_v[k][j+i+1]->y, patch_v[k][j+i+1]->z);

        glNormal3f(patch_n[k][j+i+div+1]->x, patch_n[k][j+i+div+1]->y, patch_n[k][j+i+div+1]->z);
        glVertex3f(patch_v[k][j+i+div+1]->x, patch_v[k][j+i+div+1]->y, patch_v[k][j+i+div+1]->z);

        glNormal3f(patch_n[k][j+i+div]->x, patch_n[k][j+i+div]->y, patch_n[k][j+i+div]->z);
        glVertex3f(patch_v[k][j+i+div]->x, patch_v[k][j+i+div]->y, patch_v[k][j+i+div]->z);

      }
    }
  }

  glEnd();
}

void drawSurfaceAdaptive() {
  glBegin(GL_TRIANGLES);

  for (int i = 0; i < numTriangles; i++) {
    glNormal3f(triangles[i]->n1->x, triangles[i]->n1->y, triangles[i]->n1->z);
    glVertex3f(triangles[i]->v1->x, triangles[i]->v1->y, triangles[i]->v1->z);

    glNormal3f(triangles[i]->n2->x, triangles[i]->n2->y, triangles[i]->n2->z);
    glVertex3f(triangles[i]->v2->x, triangles[i]->v2->y, triangles[i]->v2->z);

    glNormal3f(triangles[i]->n3->x, triangles[i]->n3->y, triangles[i]->n3->z);
    glVertex3f(triangles[i]->v3->x, triangles[i]->v3->y, triangles[i]->v3->z);
  }

  glEnd();
}

void outputFile() {
  ofstream obj_file;
  obj_file.open(to_output);
  if(adaptive){
    for(int i = 0; i < numTriangles; i++){
      obj_file << "v " << triangles[i]->v1->x << " " << triangles[i]->v1->y << " " << triangles[i]->v1->z << endl;
      obj_file << "v " << triangles[i]->v2->x << " " << triangles[i]->v2->y << " " << triangles[i]->v2->z << endl;
      obj_file << "v " << triangles[i]->v3->x << " " << triangles[i]->v3->y << " " << triangles[i]->v3->z << endl;
      obj_file << "vn " << triangles[i]->n1->x << " " << triangles[i]->n1->y << " " << triangles[i]->n1->z << endl;
      obj_file << "vn " << triangles[i]->n2->x << " " << triangles[i]->n2->y << " " << triangles[i]->n2->z << endl;
      obj_file << "vn " << triangles[i]->n3->x << " " << triangles[i]->n3->y << " " << triangles[i]->n3->z << endl;
      obj_file << "f " << i*3+1 << "//" << i*3+1 << " " << i*3+2 << "//" << i*3+2 << " " << i*3+3 << "//" << i*3+3 << endl;
    }
  } else {
    float div = sqrt(numVertices/numPatches);
    int offset;
    for (int k = 0; k < numPatches; k++) {
      for (int i = 0; i < div-1; i++) {
        for (int j = 0; j < (div-1)*(div-1); j+=div) {

          obj_file << "v " << patch_v[k][j+i]->x << " " << patch_v[k][j+i]->y << " " << patch_v[k][j+i]->z << endl;
          obj_file << "v " << patch_v[k][j+i+1]->x << " " << patch_v[k][j+i+1]->y << " " << patch_v[k][j+i+1]->z << endl;
          obj_file << "v " << patch_v[k][j+i+div+1]->x << " " << patch_v[k][j+i+div+1]->y << " " << patch_v[k][j+i+div+1]->z << endl;
          obj_file << "v " << patch_v[k][j+i+div]->x << " " << patch_v[k][j+i+div]->y << " " << patch_v[k][j+i+div]->z << endl;
          obj_file << "vn " << patch_n[k][j+i]->x << " " << patch_n[k][j+i]->y << " " << patch_n[k][j+i]->z << endl;
          obj_file << "vn " << patch_n[k][j+i+1]->x << " " << patch_n[k][j+i+1]->y << " " << patch_n[k][j+i+1]->z << endl;
          obj_file << "vn " << patch_n[k][j+i+div+1]->x << " " << patch_n[k][j+i+div+1]->y << " " << patch_n[k][j+i+div+1]->z << endl;
          obj_file << "vn " << patch_n[k][j+i+div]->x << " " << patch_n[k][j+i+div]->y << " " << patch_n[k][j+i+div]->z << endl;
          obj_file << "f" << k*4+1 << "//" << k*4+1 << " " << k*4+2 << "//" << k*4+2 << " " << k*4+3 << "//" << k*4+3 << " " << k*4+4 << "//" << k*4+4 << endl;

        }
      }
    }
  }
}


void drawSurfaceObj() {
  glBegin(GL_TRIANGLES);

  for (int i = 0; i < numTriangles; i++) {
    if (normals_exist) {
      glNormal3f(tris_n[i][0].x, tris_n[i][0].y, tris_n[i][0].z);
      glVertex3f(tris[i][0].x, tris[i][0].y, tris[i][0].z);

      glNormal3f(tris_n[i][1].x, tris_n[i][1].y, tris_n[i][1].z);
      glVertex3f(tris[i][1].x, tris[i][1].y, tris[i][1].z);

      glNormal3f(tris_n[i][2].x, tris_n[i][2].y, tris_n[i][2].z);
      glVertex3f(tris[i][2].x, tris[i][2].y, tris[i][2].z);
    } else {
      Point v2_v1 = *(new Point(tris[i][1].x - tris[i][0].x, tris[i][1].y - tris[i][0].y, tris[i][1].z - tris[i][0].z));
      Point v3_v1 = *(new Point(tris[i][2].x - tris[i][0].x, tris[i][2].y - tris[i][0].y, tris[i][2].z - tris[i][0].z));
      Normal n = *(new Normal (v2_v1.y*v3_v1.z-v2_v1.z*v3_v1.y, v2_v1.z*v3_v1.x-v2_v1.x*v3_v1.z, v2_v1.x*v3_v1.y-v2_v1.y*v3_v1.x));
      glNormal3f(n.x, n.y, n.z);
      glVertex3f(tris[i][0].x, tris[i][0].y, tris[i][0].z);
      glNormal3f(n.x, n.y, n.z);
      glVertex3f(tris[i][1].x, tris[i][1].y, tris[i][1].z);
      glNormal3f(n.x, n.y, n.z);
      glVertex3f(tris[i][2].x, tris[i][2].y, tris[i][2].z);
    }
  }

  glEnd();
}


//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void display( GLFWwindow* window )
{
    
    glClearColor( 0.0f, 0.0f, 0.0f, 0.0f ); //clear background screen to black
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                // clear the color buffer (sets everything to black)
    glMatrixMode(GL_MODELVIEW);                  // indicate we are specifying camera transformations
    glLoadIdentity();                            // make sure transformation is "zero'd"

    GLfloat mat_specular[] = { 0.0, 0.5, 0.8, 1.0 };
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat light_position[] = { 0.0, 0.0, -1.0, 0.0 };
    
    if (!smooth) {
     glShadeModel(GL_FLAT);
    } else {
     glShadeModel(GL_SMOOTH);
    }

    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    
    //----------------------- code to draw objects --------------------------
    glPushMatrix();
    glTranslatef (translation[0], translation[1], translation[2]);
    glScalef( scale[0], scale[1], scale[2] );
    for(int i = rotations.size()-1; i > -1; i--){
      glRotated(rotations[i][0], rotations[i][1], rotations[i][2], rotations[i][3]);
    }

    if(hidd){
      /*glCullFace(GL_FRONT);
      glPolygonMode( GL_BACK, GL_LINE);

      glEnable( GL_CULL_FACE);
      glEnable(GL_DEPTH_TEST);
      glPolygonMode( GL_FRONT, GL_FILL);
     */
     // Second pass states
      glCullFace(GL_FRONT);
      glDisable(GL_DEPTH_TEST);
      glEnable( GL_CULL_FACE);
      glPolygonMode( GL_FRONT_AND_BACK, GL_LINE);
      glDisable( GL_CULL_FACE);
    } else {
      glDisable(GL_CULL_FACE);
      glEnable(GL_DEPTH_TEST);
      if(wiref){
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      } else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      }
    }

    if (adaptive) {
      drawSurfaceAdaptive();
    } else if (obj) {
      drawSurfaceObj();
    } else {
      drawSurface();
    }

    glPopMatrix();
    
    glfwSwapBuffers(window);

    // note: check out glPolygonMode and glShadeModel 
    // for wireframe and shading commands

    
}

//****************************************************
// function that is called when window is resized
//***************************************************
void size_callback(GLFWwindow* window, int width, int height)
{
    // Get the pixel coordinate of the window
    // it returns the size, in pixels, of the framebuffer of the specified window
    glfwGetFramebufferSize(window, &Width_global, &Height_global);
    
    glViewport(0, 0, Width_global, Height_global);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-3.5, 3.5, -3.5, 3.5, 5, -5);
    
    display(window);
}

//****************************************************
// function to parse the .bez file
//***************************************************
void parse(string filename) {

  ifstream bezfile(filename);
  string line;

  if(bezfile.is_open()){

    getline(bezfile, line);
    istringstream iss(line);

    iss >> numPatches;

    for (int i = 0; i < numPatches; i++) {

      std::vector<Curve*> patch;
      for (int j = 0; j < 4; j++) {
        
        if(!getline(bezfile, line)) {
          cerr << "Missing a line!" << endl;
          glfwTerminate();
        };

        istringstream iss(line);
        std::vector<Point*> points;

        for (int k = 0; k < 4; k++) {
          
          float x, y, z;
          iss >> x;
          iss >> y;
          iss >> z;

          points.push_back(new Point(x, y, z));

        }
        patch.push_back(new Curve(points));
      }

      getline(bezfile, line);

      patches.push_back(new Patch(patch));
    }
  }
}


void parseObj(string filename) {

  std::vector<Point> vertex_list;
  std::vector<Normal> normal_list;

  ifstream objfile(filename);
  string line;

  if(objfile.is_open()){

    while(getline(objfile, line)){
      istringstream iss(line);

      string arg;
      iss >> arg;
      
      if(arg == "v"){
        float x, y, z;
        iss >> x;
        iss >> y;
        iss >> z;
        vertex_list.push_back(*(new Point(x, y, z)));
      
      } else if (arg == "vn") {
        normals_exist = true;
        float x, y, z;
        iss >> x;
        iss >> y;
        iss >> z;
        normal_list.push_back(*(new Normal(x, y, z)));
      
      } else if(arg == "f"){
        std::vector<Point*> vertices;
        std::vector<Normal*> normals;
        string val;
        float val1, val2;
        while(iss >> val){
          if(normals_exist){
            int ind = val.find("/");
            val1 = atof(val.substr(0, ind).c_str());
            string int_val = val.substr(ind+1);
            int ind1 = int_val.find("/");
            val2 = atof(int_val.substr(ind1+1).c_str());
            vertices.push_back(&(vertex_list[val1-1]));
            normals.push_back(&(normal_list[val2-1]));
          } else {
            vertices.push_back(&(vertex_list[atof(val.c_str())-1]));
          }
        }

        for(int i = 1; i < vertices.size()-1; i++){
          std::vector<Point> tri = {*(vertices[0]), *(vertices[i]), *(vertices[i+1])};
          tris.push_back(tri);
          numTriangles += 1;
          if (normals_exist) {
            std::vector<Normal> tri_n = {*(normals[0]), *(normals[i]), *(normals[i+1])};
            tris_n.push_back(tri_n);
          }
        }
      }
    }
  }
}


//****************************************************
// Given the control points of a Bezier curve and a parametric
// value, return the curve point and derivative
//***************************************************
std::vector<Point> bezCurveInterp(Curve* curve, float p) {
  // first, split each of the three segments to form two new ones AB and BC
  Point A = *(curve->p1)*(1.0f-p) + *(curve->p2)*p;
  Point B = *(curve->p2)*(1.0f-p) + *(curve->p3)*p;
  Point C = *(curve->p3)*(1.0f-p) + *(curve->p4)*p;

  // now split AB and BC to form a new segment DE
  Point D = A*(1.0f-p) + B*p;
  Point E = B*(1.0f-p) + C*p;

  // finally, pick the right point on DE, this is the point on the curve
  Point P = D*(1.0f-p) + E*p;

  //compute derivative also
  Point dPdp = (E - D) * 3.0f;

  //cout << P[0] << " " << P[1] << " " << P[2] << endl;
  return std::vector<Point> {P, dPdp};

}


void bezPatchInterp(Patch* patch, float u, float v, Point& p, Normal& n) {
  // build control points for the Bezier curve in v
  std::vector<Point> vcurve;
  
  vcurve.push_back(bezCurveInterp(patch->c1, u)[0]);
  vcurve.push_back(bezCurveInterp(patch->c2, u)[0]);
  vcurve.push_back(bezCurveInterp(patch->c3, u)[0]);
  vcurve.push_back(bezCurveInterp(patch->c4, u)[0]);

  // build control points for the Bezier curve in u
  std::vector<Point> ucurve;
  ucurve.push_back(bezCurveInterp(new Curve(patch->c1->p1, patch->c2->p1, patch->c3->p1, patch->c4->p1), v)[0]);
  ucurve.push_back(bezCurveInterp(new Curve(patch->c1->p2, patch->c2->p2, patch->c3->p2, patch->c4->p2), v)[0]);
  ucurve.push_back(bezCurveInterp(new Curve(patch->c1->p3, patch->c2->p3, patch->c3->p3, patch->c4->p3), v)[0]);
  ucurve.push_back(bezCurveInterp(new Curve(patch->c1->p4, patch->c2->p4, patch->c3->p4, patch->c4->p4), v)[0]);

  //evaluate surface and derivative for u and v
  std::vector<Point> bezCurveV = bezCurveInterp(new Curve(vcurve), v);
  std::vector<Point> bezCurveU = bezCurveInterp(new Curve(ucurve), u);
  
  p = bezCurveV[0];

  Point dPdv = bezCurveV[1];
  Point dPdu = bezCurveU[1];

  //take cross product of partials to find normal
  n = *(new Normal (dPdu.y*dPdv.z-dPdu.z*dPdv.y, dPdu.z*dPdv.x-dPdu.x*dPdv.z, dPdu.x*dPdv.y-dPdu.y*dPdv.x));

  if (n.len == 0) {
    u += 1e-5;
    v += 1e-5;
    bezPatchInterp(patch, u, v, p, n);
  }

}

//****************************************************
// function to subdivide a patch
//***************************************************

void subdividePatch(Patch* patch, float step) {
  // Compute how many subdivisions there are for this step size
  std::vector<Point*> vertices;
  std::vector<Normal*> normals;
  
  float numdiv = ((1 + EPSILON) / step);

  //for each parametric value of u
  for (int iu = 0; iu < numdiv; iu++) {
    float u = iu * step;

    //for each parametric value of v
    for (int iv = 0; iv < numdiv; iv++) {
      float v = iv * step;

      //evaluate surface
      Point* p = new Point();
      Normal* n = new Normal();
      bezPatchInterp(patch, u, v, *p, *n);

      vertices.push_back(p);
      normals.push_back(n);
      numVertices += 1;
    }
  }

  patch_v.push_back(vertices);
  patch_n.push_back(normals);

}



void adaptiveSubdividePatch(Patch* patch, float error) {

    bool subdiv = true;

  // Split patch into two triangles
  
    Point* pt1 = new Point();
    Point* pt2 = new Point();
    Point* pt3 = new Point();
    Point* pt4 = new Point();

    Normal* n1 = new Normal();
    Normal* n2 = new Normal();
    Normal* n3 = new Normal();
    Normal* n4 = new Normal();

    UV* uv1 = new UV(0, 0);
    UV* uv2 = new UV(1, 0);
    UV* uv3 = new UV(1, 1);
    UV* uv4 = new UV(0, 1);

    bezPatchInterp(patch, 0, 0, *pt1, *n1);
    bezPatchInterp(patch, 1, 0, *pt2, *n2);
    bezPatchInterp(patch, 1, 1, *pt3, *n3);
    bezPatchInterp(patch, 0, 1, *pt4, *n4);

    std::vector<Triangle*> initial_triangles;
    int initial_numTriangles = 2;

    initial_triangles.push_back(new Triangle(pt1, pt2, pt3, n1, n2, n3, uv1, uv2, uv3));
    initial_triangles.push_back(new Triangle(pt1, pt3, pt4, n1, n3, n4, uv1, uv3, uv4));

    while(subdiv) {
    
      std::vector<Triangle*> new_triangles;
      int new_numTriangles = 0;
      subdiv = false;

      for (int i = 0; i < initial_numTriangles; i++) {
        bool edge1, edge2, edge3 = false;

        Triangle* triangle = initial_triangles[i];
        
        Point* real_pt1 = triangle->vert_midpoint1();
        UV* uv1 = triangle->uv_midpoint1();
        Point* mdpt1 = new Point();
        Normal* mdn1 = new Normal();
        bezPatchInterp(patch, uv1->u, uv1->v, *mdpt1, *mdn1);
        edge1 = (sqrt(sqr((*mdpt1-*real_pt1).x) + sqr((*mdpt1-*real_pt1).y) + sqr((*mdpt1-*real_pt1).z)) > error);

        Point* real_pt2 = triangle->vert_midpoint2();
        UV* uv2 = triangle->uv_midpoint2();
        Point* mdpt2 = new Point();
        Normal* mdn2 = new Normal();
        bezPatchInterp(patch, uv2->u, uv2->v, *mdpt2, *mdn2);
        edge2 = (sqrt(sqr((*mdpt2-*real_pt2).x) + sqr((*mdpt2-*real_pt2).y) + sqr((*mdpt2-*real_pt2).z)) > error);
  
        Point* real_pt3 = triangle->vert_midpoint3();
        UV* uv3 = triangle->uv_midpoint3();
        Point* mdpt3 = new Point();
        Normal* mdn3 = new Normal();
        bezPatchInterp(patch, uv3->u, uv3->v, *mdpt3, *mdn3);
        edge3 = (sqrt(sqr((*mdpt3-*real_pt3).x) + sqr((*mdpt3-*real_pt3).y) + sqr((*mdpt3-*real_pt3).z)) > error);

        if (edge1 && edge2 && edge3) {
          new_triangles.push_back(new Triangle(triangle->v1, mdpt1, mdpt3, triangle->n1, mdn1, mdn3, triangle->uv1, uv1, uv3));
          new_triangles.push_back(new Triangle(mdpt1, triangle->v2, mdpt2, mdn1, triangle->n2, mdn2, uv1, triangle->uv2, uv2));
          new_triangles.push_back(new Triangle(mdpt3, mdpt2, triangle->v3, mdn3, mdn2, triangle->n3, uv3, uv2, triangle->uv3));
          new_triangles.push_back(new Triangle(mdpt2, mdpt1, mdpt3, mdn2, mdn1, mdn3, uv2, uv1, uv3));
          new_numTriangles += 4;
          subdiv = true;
        } else if (edge1 && edge2) {
          new_triangles.push_back(new Triangle(triangle->v1, triangle->v3, mdpt2, triangle->n1, triangle->n3, mdn2, triangle->uv1, triangle->uv3, uv2));
          new_triangles.push_back(new Triangle(triangle->v1, mdpt2, mdpt1, triangle->n1, mdn2, mdn1, triangle->uv1, uv2, uv1));
          new_triangles.push_back(new Triangle(triangle->v2, mdpt2, mdpt1, triangle->n2, mdn2, mdn1, triangle->uv2, uv2, uv1));
          new_numTriangles += 3;
          subdiv = true;
        } else if (edge2 && edge3) {
          new_triangles.push_back(new Triangle(triangle->v1, mdpt3, triangle->v2, triangle->n1, mdn3, triangle->n2, triangle->uv1, uv3, triangle->uv2));
          new_triangles.push_back(new Triangle(mdpt2, mdpt3, triangle->v2, mdn2, mdn3, triangle->n2, uv2, uv3, triangle->uv2));
          new_triangles.push_back(new Triangle(mdpt2, mdpt3, triangle->v3, mdn2, mdn3, triangle->n3, uv2, uv3, triangle->uv3));
          new_numTriangles += 3;
          subdiv = true;
        } else if (edge1 && edge3) {
          new_triangles.push_back(new Triangle(triangle->v1, mdpt3, mdpt1, triangle->n1, mdn3, mdn1, triangle->uv1, uv3, uv1));
          new_triangles.push_back(new Triangle(triangle->v3, mdpt3, mdpt1, triangle->n3, mdn3, mdn1, triangle->uv3, uv3, uv1));
          new_triangles.push_back(new Triangle(triangle->v2, triangle->v3, mdpt1, triangle->n2, triangle->n3, mdn1, triangle->uv2, triangle->uv3, uv1));
          new_numTriangles += 3;
          subdiv = true;
        } else if (edge1) {
          new_triangles.push_back(new Triangle(triangle->v2, triangle->v3, mdpt1, triangle->n2, triangle->n3, mdn1, triangle->uv2, triangle->uv3, uv1));
          new_triangles.push_back(new Triangle(triangle->v1, triangle->v3, mdpt1, triangle->n1, triangle->n3, mdn1, triangle->uv1, triangle->uv3, uv1));
          new_numTriangles += 2;
          subdiv = true;
        } else if (edge2) {
          new_triangles.push_back(new Triangle(triangle->v1, triangle->v3, mdpt2, triangle->n1, triangle->n3, mdn2, triangle->uv1, triangle->uv3, uv2));
          new_triangles.push_back(new Triangle(triangle->v1, triangle->v2, mdpt2, triangle->n1, triangle->n2, mdn2, triangle->uv1, triangle->uv2, uv2));
          new_numTriangles += 2;
          subdiv = true;
        } else if (edge3) {
          new_triangles.push_back(new Triangle(triangle->v1, mdpt3, triangle->v2, triangle->n1, mdn3, triangle->n2, triangle->uv1, uv3, triangle->uv2));
          new_triangles.push_back(new Triangle(triangle->v3, mdpt3, triangle->v2, triangle->n3, mdn3, triangle->n2, triangle->uv3, uv3, triangle->uv2));
          new_numTriangles += 2;
          subdiv = true;
        } else {
          triangles.push_back(triangle);
          numTriangles += 1;
        }
      }

    initial_triangles.clear();
    initial_triangles = new_triangles;
    initial_numTriangles = new_numTriangles;

    }
}


//****************************************************
// the usual stuff, nothing exciting here
//****************************************************
int main(int argc, char *argv[]) {
    //This initializes glfw
    initializeRendering();

    if (argc == 1) {
        cerr << "Need input file" << endl;
        glfwTerminate();
        return -1;
    } else if (argc == 2) {
        cerr << "Need subdivision parameter" << endl;
        glfwTerminate();
        return -1;
    } else {
      
      string filename = argv[1];
      const char* o = "obj";
      const char* ext = filename.substr(filename.size()-3, 3).c_str();
      if (!strcmp(ext, o)) {
        parseObj(filename);
        obj = true;
      } else {
        parse(filename);
        float step = atof(argv[2]);
        if (argc == 4) {
          string flag = argv[3];
          if (flag == "-a") {
            adaptive = true;
            for (int i = 0; i < numPatches; i++) {
              adaptiveSubdividePatch(patches[i], step);
            }
          } else {
            cerr << "Incorrect Subdivision Flag. Valid options: -a" << endl;
            glfwTerminate();
            return -1;
          }
        } else if (argc == 5){
          string flag1 = argv[3];
          string flag2 = argv[4];
          if(flag1 == "-o"){
            output = true;
            to_output = flag2;
          } else {
            cerr << "Incorrect Output Flag. Valid Options: -o" << endl;
            glfwTerminate();
            return -1;
          }
        } else if (argc == 6){
          string flag1 = argv[3];
          string flag2 = argv[4];
          string flag3 = argv[5];
          if(flag1 == "-a"){
            adaptive = true;
            for(int i = 0; i < numPatches; i++){
              adaptiveSubdividePatch(patches[i], step);
            }
          } else {
            cerr << "Incorrect Subdivision Flag. Valid options: -a" << endl;
            glfwTerminate();
            return -1;
          } if (flag2 == "-o"){
            output = true;
            to_output = flag3;
          } else {
            cerr << "Incorrect Output Flag. Valid Options: -o" << endl;
            glfwTerminate();
            return -1;
          }

        } else {
          for (int i = 0; i < numPatches; i++) {
            subdividePatch(patches[i], step);
          }
        }
      }
    }
    
    GLFWwindow* window = glfwCreateWindow( Width_global, Height_global, "CS184", NULL, NULL );
    if ( !window )
    {
        cerr << "Error on window creating" << endl;
        glfwTerminate();
        return -1;
    }

    const GLFWvidmode * mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    if ( !mode )
    {
        cerr << "Error on getting monitor" << endl;
        glfwTerminate();
        return -1;
    }
    
    glfwMakeContextCurrent( window );
    
    // Get the pixel coordinate of the window
    // it returns the size, in pixels, of the framebuffer of the specified window
    glfwGetFramebufferSize(window, &Width_global, &Height_global);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
  	glOrtho(-3.5, 3.5, -3.5, 3.5, 5, -5);
  	glMatrixMode(GL_MODELVIEW);
  	glLoadIdentity();

    glEnable(GL_DEPTH_TEST);	// enable z-buffering
    glDepthFunc(GL_LESS);

    glfwSetWindowTitle(window, "CS184");
    glfwSetWindowSizeCallback(window, size_callback);
    glfwSetKeyCallback(window, key_callback);

    if(output){
      outputFile();
    }

    while( !glfwWindowShouldClose( window ) ) // infinite loop to draw object again and again
    {   // because once object is draw then window is terminated
        display( window );
        
        if (auto_strech){
            glfwSetWindowSize(window, mode->width, mode->height);
            glfwSetWindowPos(window, 0, 0);
        }
        
        glfwPollEvents();
        
    }

    return 0;
}