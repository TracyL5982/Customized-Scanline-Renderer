#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include "SETTINGS.h"

using namespace std;

int xScreenRes = 800;
int yScreenRes = 600;

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
VEC3 truncate(const VEC4& v)
{
  return VEC3(v[0], v[1], v[2]);
}
VEC4 extend(const VEC3& v)
{
  return VEC4(v[0], v[1], v[2], 1.0);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void readPPM(const string& filename, int& xRes, int& yRes, float*& values)
{
  // try to open the file
  FILE *fp;
  fp = fopen(filename.c_str(), "rb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for reading." << endl;
    cout << " Make sure you're not trying to read from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  // get the dimensions
  unsigned char newline;
  fscanf(fp, "P6\n%d %d\n255%c", &xRes, &yRes, &newline);
  if (newline != '\n') {
    cout << " The header of " << filename.c_str() << " may be improperly formatted." << endl;
    cout << " The program will continue, but you may want to check your input. " << endl;
  }
  int totalCells = xRes * yRes;

  // grab the pixel values
  unsigned char* pixels = new unsigned char[3 * totalCells];
  fread(pixels, 1, totalCells * 3, fp);

  // copy to a nicer data type
  values = new float[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = pixels[i];

  // clean up
  delete[] pixels;
  fclose(fp);
  cout << " Read in file " << filename.c_str() << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void writePPM(const string& filename, int& xRes, int& yRes, float*& values)
{
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    pixels[i] = values[i];

  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for writing." << endl;
    cout << " Make sure you're not trying to write from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
  fwrite(pixels, 1, totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}

//////////////////////////////////////////////////////////////////////////////////
// build out a single square
//////////////////////////////////////////////////////////////////////////////////
void buildBigSquare(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
{
  vertices.push_back(VEC3(1.0, 1.0, 1.5));
  vertices.push_back(VEC3( 11.0, 1.0, 1.5));
  vertices.push_back(VEC3(1.0,  11.0, 1.5));
  vertices.push_back(VEC3( 11.0,  11.0, 1.5));

  // front face
  indices.push_back(VEC3I(0, 1, 2));
  indices.push_back(VEC3I(2, 1, 3));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(0.0, 0.0, 1.0));
}

//////////////////////////////////////////////////////////////////////////////////
// build out a single square
//////////////////////////////////////////////////////////////////////////////////
void buildSquare(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
{
  vertices.push_back(VEC3(-0.5, -0.5,  0.5));
  vertices.push_back(VEC3( 0.5, -0.5,  0.5));
  vertices.push_back(VEC3(-0.5,  0.5,  0.5));
  vertices.push_back(VEC3( 0.5,  0.5,  0.5));

  // front face
  indices.push_back(VEC3I(0, 1, 2));
  indices.push_back(VEC3I(2, 1, 3));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(0.0, 1.0, 0.0));
}

//////////////////////////////////////////////////////////////////////////////////
// build out a cube
//////////////////////////////////////////////////////////////////////////////////
void buildCube(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
{
  vertices.push_back(VEC3(-0.5, -0.5,  0.5));
  vertices.push_back(VEC3( 0.5, -0.5,  0.5));
  vertices.push_back(VEC3(-0.5,  0.5,  0.5));
  vertices.push_back(VEC3( 0.5,  0.5,  0.5));
  vertices.push_back(VEC3(-0.5, -0.5, -0.5));
  vertices.push_back(VEC3( 0.5, -0.5, -0.5));
  vertices.push_back(VEC3(-0.5,  0.5, -0.5));
  vertices.push_back(VEC3( 0.5,  0.5, -0.5));

  // front face
  indices.push_back(VEC3I(0, 1, 2));
  indices.push_back(VEC3I(2, 1, 3));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(1.0, 0.0, 0.0));

  // back face
  indices.push_back(VEC3I(5, 4, 7));
  indices.push_back(VEC3I(7, 4, 6));
  colors.push_back(VEC3(0.0, 1.0, 0.0));
  colors.push_back(VEC3(0.0, 1.0, 0.0));

  // left face
  indices.push_back(VEC3I(4, 0, 6));
  indices.push_back(VEC3I(6, 0, 2));
  colors.push_back(VEC3(0.0, 0.0, 1.0));
  colors.push_back(VEC3(0.0, 0.0, 1.0));

  // right face
  indices.push_back(VEC3I(1, 5, 3));
  indices.push_back(VEC3I(3, 5, 7));
  colors.push_back(VEC3(0.0, 1.0, 1.0));
  colors.push_back(VEC3(0.0, 1.0, 1.0));

  // top face
  indices.push_back(VEC3I(2, 3, 6));
  indices.push_back(VEC3I(6, 3, 7));
  colors.push_back(VEC3(1.0, 1.0, 0.0));
  colors.push_back(VEC3(1.0, 1.0, 0.0));

  // bottom face
  indices.push_back(VEC3I(4, 5, 0));
  indices.push_back(VEC3I(0, 5, 1));
  colors.push_back(VEC3(1.0, 0.0, 1.0));
  colors.push_back(VEC3(1.0, 0.0, 1.0));
}

//////////////////////////////////////////////////////////////////////////////////
// build out a cube
//////////////////////////////////////////////////////////////////////////////////
void buildCubePerVertexColors(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
{
  vertices.push_back(VEC3(-0.5, -0.5,  0.5));
  vertices.push_back(VEC3( 0.5, -0.5,  0.5));
  vertices.push_back(VEC3(-0.5,  0.5,  0.5));
  vertices.push_back(VEC3( 0.5,  0.5,  0.5));
  vertices.push_back(VEC3(-0.5, -0.5, -0.5));
  vertices.push_back(VEC3( 0.5, -0.5, -0.5));
  vertices.push_back(VEC3(-0.5,  0.5, -0.5));
  vertices.push_back(VEC3( 0.5,  0.5, -0.5));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(0.0, 1.0, 0.0));
  colors.push_back(VEC3(0.0, 1.0, 0.0));
  colors.push_back(VEC3(0.0, 0.0, 1.0));
  colors.push_back(VEC3(0.0, 0.0, 1.0));
  colors.push_back(VEC3(1.0, 1.0, 0.0));
  colors.push_back(VEC3(1.0, 1.0, 0.0));

  // front face
  indices.push_back(VEC3I(0, 1, 2));
  indices.push_back(VEC3I(2, 1, 3));

  // back face
  indices.push_back(VEC3I(5, 4, 7));
  indices.push_back(VEC3I(7, 4, 6));

  // left face
  indices.push_back(VEC3I(4, 0, 6));
  indices.push_back(VEC3I(6, 0, 2));

  // right face
  indices.push_back(VEC3I(1, 5, 3));
  indices.push_back(VEC3I(3, 5, 7));

  // top face
  indices.push_back(VEC3I(2, 3, 6));
  indices.push_back(VEC3I(6, 3, 7));

  // bottom face
  indices.push_back(VEC3I(4, 5, 0));
  indices.push_back(VEC3I(0, 5, 1));
}

//STEPS:
//////////////////////////////////////////////////////////////////////////////////
// 1. Computing Mvp
//////////////////////////////////////////////////////////////////////////////////
MATRIX4 applyViewportMatrix(){
  MATRIX4 mvp;
  mvp.setZero();
  mvp(0,0) = xScreenRes/2;
  mvp(1,1) = yScreenRes/2;
  mvp(2,2) = 1;
  mvp(3,3) = 1;
  mvp(0,3) = (xScreenRes-1)/2;
  mvp(1,3) = (yScreenRes-1)/2;
  return mvp;
}

//////////////////////////////////////////////////////////////////////////////////
// 2. Basic Triangle Rasterization
//////////////////////////////////////////////////////////////////////////////////
double barycentric_calculation (double x0, double y0, double x1, double y1, double x, double y){
  return (y0 - y1) * x + (x1 - x0) * y + x0 * y1 - x1 * y0;
}

void rasterization (float* &values, vector<VEC4>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors, int surface_index){
  VEC3 curr_color = colors[surface_index];
  double x0 = vertices[indices[surface_index][0]][0];
  double y0 = vertices[indices[surface_index][0]][1];
  double x1 = vertices[indices[surface_index][1]][0];
  double y1 = vertices[indices[surface_index][1]][1];
  double x2 = vertices[indices[surface_index][2]][0];
  double y2 = vertices[indices[surface_index][2]][1];

  int coloredPixels = 0;
  double f_alpha = barycentric_calculation(x1,y1,x2,y2,x0,y0);
  double f_beta = barycentric_calculation(x2,y2,x0,y0,x1,y1);
  double f_gamma = barycentric_calculation(x0,y0,x1,y1,x2,y2);

  for (int y = 0; y < yScreenRes; y++){
    for (int x = 0; x < xScreenRes; x++){
      double alpha = barycentric_calculation(x1,y1,x2,y2,x,y)/f_alpha;
      double beta = barycentric_calculation(x2,y2,x0,y0,x,y)/f_beta;
      double gamma = barycentric_calculation(x0,y0,x1,y1,x,y)/f_gamma;
      if((alpha >= 0) && (beta >= 0) && (gamma >= 0)){
        int index = (yScreenRes -1 - y) * xScreenRes + x;
        values[index*3] = curr_color[0] * 255.0f;
        values[index*3+1] = curr_color[1] * 255.0f;
        values[index*3+2] = curr_color[2] * 255.0f;
        coloredPixels++;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
// 3. Computing M_orth
//////////////////////////////////////////////////////////////////////////////////
MATRIX4 applyOrthMatrix(double l, double r, double t, double b, double n, double f){
  MATRIX4 m_orth;
  m_orth.setZero();
  m_orth(0,0) = 2/(r-l);
  m_orth(1,1) = 2/(t-b);
  m_orth(2,2) = 2/(n-f);
  m_orth(3,3) = 1;
  m_orth(0,3) = 0 - (r+l)/(r-l);
  m_orth(1,3) = 0 - (t+b)/(t-b);
  m_orth(2,3) = 0 - (n+f)/(n-f);
  return m_orth;
}

//////////////////////////////////////////////////////////////////////////////////
// 4. Computing M_cam
//////////////////////////////////////////////////////////////////////////////////
MATRIX4 applyCamMatrix(VEC3 e, VEC3 g, VEC3 t){
  VEC3 w = -g.normalized();
  VEC3 u = t.cross(w);
  u.normalize ();
  VEC3 v = w.cross(u);
  MATRIX4 cam_1, cam_2;
  cam_1.setZero();
  cam_2.setZero();

  cam_1(0,0) = u[0];
  cam_1(0,1) = u[1];
  cam_1(0,2) = u[2];
  cam_1(1,0) = v[0];
  cam_1(1,1) = v[1];
  cam_1(1,2) = v[2];
  cam_1(2,0) = w[0];
  cam_1(2,1) = w[1];
  cam_1(2,2) = w[2];
  cam_1(3,3) = 1;

  cam_2(0,0) = 1;
  cam_2(1,1) = 1;
  cam_2(2,2) = 1;
  cam_2(3,3) = 1;
  cam_2(0,3) = -e[0];
  cam_2(1,3) = -e[1];
  cam_2(2,3) = -e[2];
  return cam_1 * cam_2;
}

//////////////////////////////////////////////////////////////////////////////////
// 5. Computing M_pers
//////////////////////////////////////////////////////////////////////////////////
MATRIX4 applyPersMatrix(float fovy, float aspect, float n, float f){
  float scale = tan(fovy * 0.5 * M_PI / 180);
  float t = scale * n;
  float r = scale * aspect * n;
  float l = -r;
  float b = -t;
  MATRIX4 m_pers;
  m_pers.setZero ();
  m_pers(0,0) = 2*n/(r-l);
  m_pers(0,2) = (l+r)/(r-l);
  m_pers(1,1) = 2*n/(t-b);
  m_pers(1,2) = (b+t)/(t-b);
  m_pers(2,2) = (f+n)/(n-f);
  m_pers(2,3) = -(2*f*n)/(f-n);
  m_pers(3,2) = -1;
  return m_pers;
}

//////////////////////////////////////////////////////////////////////////////////
// 7. Rasterization with z buffering
//////////////////////////////////////////////////////////////////////////////////
void z_buffer (float* &values, vector<VEC4>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors, int surface_index, vector<double>& depth_storage){
  VEC3 curr_color = colors[surface_index];
  double x0 = vertices[indices[surface_index][0]][0];
  double y0 = vertices[indices[surface_index][0]][1];
  double x1 = vertices[indices[surface_index][1]][0];
  double y1 = vertices[indices[surface_index][1]][1];
  double x2 = vertices[indices[surface_index][2]][0];
  double y2 = vertices[indices[surface_index][2]][1];
  double z0 = vertices[indices[surface_index][0]][2];
  double z1 = vertices[indices[surface_index][1]][2];
  double z2 = vertices[indices[surface_index][2]][2];

  //int coloredPixels = 0;
  double f_alpha = barycentric_calculation(x1,y1,x2,y2,x0,y0);
  double f_beta = barycentric_calculation(x2,y2,x0,y0,x1,y1);
  double f_gamma = barycentric_calculation(x0,y0,x1,y1,x2,y2);

  for (int y = 0; y < yScreenRes; y++){
    for (int x = 0; x < xScreenRes; x++){
      double alpha = barycentric_calculation(x1,y1,x2,y2,x,y)/f_alpha;
      double beta = barycentric_calculation(x2,y2,x0,y0,x,y)/f_beta;
      double gamma = barycentric_calculation(x0,y0,x1,y1,x,y)/f_gamma;
      if((alpha >= 0) && (beta >= 0) && (gamma >= 0)){
        double depth = alpha * z0 + beta * z1 + gamma * z2;
        int index = (yScreenRes -1 - y) * xScreenRes + x;
        if(depth < depth_storage[index]){
          depth_storage[index] = depth;
          values[index*3] = curr_color[0] * 255.0f;
          values[index*3+1] = curr_color[1] * 255.0f;
          values[index*3+2] = curr_color[2] * 255.0f;
          //coloredPixels++;
        }
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////
// 8. Rasterization with color interpolation
//////////////////////////////////////////////////////////////////////////////////
void color_interp (float* &values, vector<VEC4>& vertices, vector<VEC3I>& indices, vector<VEC3>& vertexColors, int surface_index, vector<double>& depth_storage){
  double x0 = vertices[indices[surface_index][0]][0];
  double y0 = vertices[indices[surface_index][0]][1];
  double x1 = vertices[indices[surface_index][1]][0];
  double y1 = vertices[indices[surface_index][1]][1];
  double x2 = vertices[indices[surface_index][2]][0];
  double y2 = vertices[indices[surface_index][2]][1];
  double z0 = vertices[indices[surface_index][0]][2];
  double z1 = vertices[indices[surface_index][1]][2];
  double z2 = vertices[indices[surface_index][2]][2];
  VEC3 color0 = vertexColors[indices[surface_index][0]];
  VEC3 color1 = vertexColors[indices[surface_index][1]];
  VEC3 color2 = vertexColors[indices[surface_index][2]];

  //int coloredPixels = 0;
  double f_alpha = barycentric_calculation(x1,y1,x2,y2,x0,y0);
  double f_beta = barycentric_calculation(x2,y2,x0,y0,x1,y1);
  double f_gamma = barycentric_calculation(x0,y0,x1,y1,x2,y2);

  for (int y = 0; y < yScreenRes; y++){
    for (int x = 0; x < xScreenRes; x++){
      double alpha = barycentric_calculation(x1,y1,x2,y2,x,y)/f_alpha;
      double beta = barycentric_calculation(x2,y2,x0,y0,x,y)/f_beta;
      double gamma = barycentric_calculation(x0,y0,x1,y1,x,y)/f_gamma;
      if((alpha >= 0) && (beta >= 0) && (gamma >= 0)){
        double depth = alpha * z0 + beta * z1 + gamma * z2;
        VEC3 interpolatedColor = alpha * color0 + beta * color1 + gamma * color2;
        int index = (yScreenRes -1 - y) * xScreenRes + x;
        if(depth < depth_storage[index]){
          depth_storage[index] = depth;
          values[index*3]   = (alpha * color0[0]+ beta * color1[0] + gamma * color2[0]) * 255.0f;
          values[index*3+1] = (alpha * color0[1]+ beta * color1[1] + gamma * color2[1]) * 255.0f;
          values[index*3+2] = (alpha * color0[2]+ beta * color1[2] + gamma * color2[2]) * 255.0f;
          //coloredPixels++;
        }
      }
    }
  }
}

//MAKE PPM
//////////////////////////////////////////////////////////////////////////////////
// Fig.1. producing 1.ppm
//////////////////////////////////////////////////////////////////////////////////
void make_1_ppm(){
  MATRIX4 mvp = applyViewportMatrix();
  vector<VEC3> vertices;
  vector<VEC3I> indices;
  vector<VEC3> colors;
  vector<VEC4> new_vertices;
  buildSquare(vertices, indices, colors);

  for(auto vertex : vertices){
    VEC4 new_vertex = extend(vertex);
    new_vertex = mvp * new_vertex;
    new_vertices.push_back(new_vertex);
  }

  string filename = "1.ppm";
  int totalCells = xScreenRes * yScreenRes;
  float* values = new float[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  for (int i = 0; i < indices.size(); i++){
    for(int j = 0; j < indices[0].size(); j++){
      VEC4 curr_vertex = new_vertices[indices[i][j]];
      int x = curr_vertex[0];
      int y = yScreenRes - curr_vertex[1];
      int index = y * xScreenRes + x;
      values[index * 3] = colors[i][0] * 255.0f;
      values[index * 3+1] = colors[i][1] * 255.0f;
      values[index * 3+2] = colors[i][2] * 255.0f;
    }
  }
  writePPM(filename, xScreenRes, yScreenRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// Fig.2. producing 2.ppm
//////////////////////////////////////////////////////////////////////////////////
void make_2_ppm(){
  MATRIX4 mvp = applyViewportMatrix();
  vector<VEC3> vertices;
  vector<VEC3I> indices;
  vector<VEC3> colors;
  vector<VEC4> new_vertices;
  buildSquare(vertices, indices, colors);

  for(auto vertex : vertices){
    VEC4 new_vertex = extend(vertex);
    new_vertex = mvp * new_vertex;
    new_vertices.push_back(new_vertex);
  }

  string filename = "2.ppm";
  int totalCells = xScreenRes * yScreenRes;
  float* values = new float[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  rasterization(values, new_vertices, indices, colors, 0);
  rasterization(values, new_vertices, indices, colors, 1);

  writePPM(filename, xScreenRes, yScreenRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// Fig.3. producing 3.ppm
//////////////////////////////////////////////////////////////////////////////////
void make_3_ppm(){
  MATRIX4 mvp = applyViewportMatrix();
  MATRIX4 m_orth = applyOrthMatrix(0.0, 12.0, 12.0, 0.0, 12.0, 0.0);
  vector<VEC3> vertices;
  vector<VEC3I> indices;
  vector<VEC3> colors;
  vector<VEC4> new_vertices;
  buildBigSquare(vertices, indices, colors);

  for(auto vertex : vertices){
    VEC4 new_vertex = extend(vertex);
    new_vertex = mvp * m_orth * new_vertex;
    new_vertices.push_back(new_vertex);
  }

  string filename = "3.ppm";
  int totalCells = xScreenRes * yScreenRes;
  float* values = new float[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  rasterization(values, new_vertices, indices, colors, 0);
  rasterization(values, new_vertices, indices, colors, 1);

  writePPM(filename, xScreenRes, yScreenRes, values);

}

//////////////////////////////////////////////////////////////////////////////////
// Fig.4. producing 4.ppm
//e = (0.2, 0.2, 1), lookat = (0, 0, 0), and t = (0, 1, 0)
//////////////////////////////////////////////////////////////////////////////////
void make_4_ppm(){
  //mvp
  MATRIX4 mvp = applyViewportMatrix();

  //m_orth
  MATRIX4 m_orth = applyOrthMatrix(0.0, 12.0, 12.0, 0.0, 12.0, 0.0);

  //m_cam
  VEC3 e(0.2, 0.2, 1.0);
  VEC3 target(0.0, 0.0, 0.0);
  VEC3 g = target - e;
  VEC3 t(0.0, 1.0, 0.0);
  MATRIX m_cam = applyCamMatrix(e, g, t);

  //build big square
  vector<VEC3> vertices;
  vector<VEC3I> indices;
  vector<VEC3> colors;
  vector<VEC4> new_vertices;
  buildBigSquare(vertices, indices, colors);

  for(auto vertex : vertices){
    VEC4 new_vertex = extend(vertex);
    new_vertex = mvp * m_orth * m_cam * new_vertex;
    new_vertices.push_back(new_vertex);
  }

  string filename = "4.ppm";
  int totalCells = xScreenRes * yScreenRes;
  float* values = new float[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  rasterization(values, new_vertices, indices, colors, 0);
  rasterization(values, new_vertices, indices, colors, 1);

  writePPM(filename, xScreenRes, yScreenRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// Fig.5. producing 5.ppm
// M = Mvp * M_orth * P * M_cam
// all the vertices scaled by 0.5
// fovy = 65.0, aspect = 4.0 / 3.0, near = 1.0, far = 100.0
// eye = (1,1,1), lookAt = (0,0,0), up = (0,1,0).
//////////////////////////////////////////////////////////////////////////////////
void make_5_ppm(){
  //mvp
  MATRIX4 mvp = applyViewportMatrix();

  //m_cam
  VEC3 e(1.0, 1.0, 1.0);
  VEC3 target(0.0, 0.0, 0.0);
  VEC3 g = target - e;
  VEC3 t(0.0, 1.0, 0.0);
  MATRIX4 m_cam = applyCamMatrix(e, g, t);

  //m_orth
  double fovy = 65.0;
  double aspect = 4.0/3.0;
  double near = 1.0;
  double far = 100.0;
  double scale = tan(fovy * 0.5 * M_PI / 180);
  double top = scale * near;
  double right = scale * aspect * near;
  double left = -right;
  double bottom = -top;
  MATRIX4 m_orth = applyOrthMatrix(left, right, top, bottom, near, far);

  //build cube
  vector<VEC3> vertices;
  vector<VEC3I> indices;
  vector<VEC3> colors;
  vector<VEC4> new_vertices;
  buildCube(vertices, indices, colors);

  //scaling
  MATRIX4 m_scale;
  m_scale.setZero();
  m_scale(0,0) = 0.5;
  m_scale(1,1) = 0.5;
  m_scale(2,2) = 0.5;
  m_scale(3,3) = 1.0;

  for(auto vertex : vertices){
    VEC4 new_vertex = extend(vertex);
    new_vertex = m_scale * new_vertex;
    new_vertex = mvp * m_orth * m_cam * new_vertex;
    new_vertices.push_back(new_vertex);
  }

  string filename = "5.ppm";
  int totalCells = xScreenRes * yScreenRes;
  float* values = new float[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  for(int i = 0; i < 12; i++){
    rasterization(values, new_vertices, indices, colors, i);
  }

  writePPM(filename, xScreenRes, yScreenRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// Fig.6. producing 6.ppm.
// perspective divide
//////////////////////////////////////////////////////////////////////////////////
void make_6_ppm(){
  //mvp
  MATRIX4 mvp = applyViewportMatrix();

  //m_cam
  VEC3 e(1.0, 1.0, 1.0);
  VEC3 target(0.0, 0.0, 0.0);
  VEC3 g = target - e;
  VEC3 t(0.0, 1.0, 0.0);
  MATRIX4 m_cam = applyCamMatrix(e, g, t);

  //m_pers
  MATRIX4 m_pers = applyPersMatrix(65.0, 4.0/3.0, 1.0, 100.0);

  vector<VEC3> vertices;
  vector<VEC3I> indices;
  vector<VEC3> colors;
  vector<VEC4> new_vertices;
  buildCube(vertices, indices, colors);

  for(auto vertex : vertices){
    VEC4 new_vertex = extend(vertex);
    new_vertex = mvp * m_pers * m_cam * new_vertex;
    // applying perspective divide
    if (new_vertex[3] != 0.0) {
        new_vertex[0] /= new_vertex[3];
        new_vertex[1] /= new_vertex[3];
        new_vertex[2] /= new_vertex[3];
    }
    new_vertices.push_back(new_vertex);
  }

  string filename = "6.ppm";
  int totalCells = xScreenRes * yScreenRes;
  float* values = new float[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  for(int i = 0; i < 12; i++){
    rasterization(values, new_vertices, indices, colors, i);
  }

  writePPM(filename, xScreenRes, yScreenRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// Fig.7. producing 7.ppm.
// Z-buffering
//////////////////////////////////////////////////////////////////////////////////
void make_7_ppm(){
  //mvp
  MATRIX4 mvp = applyViewportMatrix();

  //m_cam
  VEC3 e(1.0, 1.0, 1.0);
  VEC3 target(0.0, 0.0, 0.0);
  VEC3 g = target - e;
  VEC3 t(0.0, 1.0, 0.0);
  MATRIX4 m_cam = applyCamMatrix(e, g, t);

  //m_pers
  MATRIX4 m_pers = applyPersMatrix(65.0, 4.0/3.0, 1.0, 100.0);

  //build cube
  vector<VEC3> vertices;
  vector<VEC3I> indices;
  vector<VEC3> colors;
  vector<VEC4> new_vertices;
  buildCube(vertices, indices, colors);

  for(auto vertex : vertices){
    VEC4 new_vertex = extend(vertex);
    new_vertex = mvp * m_pers * m_cam * new_vertex;

    if (new_vertex[3] != 0.0) {
      new_vertex[0] /= new_vertex[3];
      new_vertex[1] /= new_vertex[3];
      new_vertex[2] /= new_vertex[3];
    }

    new_vertices.push_back(new_vertex);
  }

  string filename = "7.ppm";
  int totalCells = xScreenRes * yScreenRes;
  float* values = new float[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  vector<double> depth_storage (totalCells);

  for (int i = 0; i < totalCells; i++)
    depth_storage[i] = numeric_limits<double>::max();

  for(int i = 0; i < 12; i++){
    z_buffer(values, new_vertices, indices, colors, i, depth_storage);
  }

  writePPM(filename, xScreenRes, yScreenRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// Fig.8. producing 8.ppm.
// Color interpolation
//////////////////////////////////////////////////////////////////////////////////
void make_8_ppm(){
  //mvp
  MATRIX4 mvp = applyViewportMatrix();

  //m_cam
  VEC3 e(1.0, 1.0, 1.0);
  VEC3 target(0.0, 0.0, 0.0);
  VEC3 g = target - e;
  VEC3 t(0.0, 1.0, 0.0);
  MATRIX4 m_cam = applyCamMatrix(e, g, t);

  //m_pers
  MATRIX4 m_pers = applyPersMatrix(65.0, 4.0/3.0, 1.0, 100.0);

  //build cube
  vector<VEC3> vertices;
  vector<VEC3I> indices;
  vector<VEC3> colors;
  vector<VEC4> new_vertices;
  buildCubePerVertexColors(vertices, indices, colors);

  for(auto vertex : vertices){
    VEC4 new_vertex = extend(vertex);
    new_vertex = mvp * m_pers * m_cam * new_vertex;

    if (new_vertex[3] != 0.0) {
      new_vertex[0] /= new_vertex[3];
      new_vertex[1] /= new_vertex[3];
      new_vertex[2] /= new_vertex[3];
    }

    new_vertices.push_back(new_vertex);
  }

  string filename = "8.ppm";
  int totalCells = xScreenRes * yScreenRes;
  float* values = new float[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  vector<double> depth_storage (totalCells);

  for (int i = 0; i < totalCells; i++)
    depth_storage[i] = numeric_limits<double>::max();

  for(int i = 0; i < 12; i++){
    color_interp(values, new_vertices, indices, colors, i, depth_storage);
  }

  writePPM(filename, xScreenRes, yScreenRes, values);
}

//////////////////////////////////////////////////////////////////////////////////
// Fig.custom. producing custom.ppm.
// taking in 9 parameters
//////////////////////////////////////////////////////////////////////////////////
void make_custom_ppm(VEC3 e, VEC3 g, VEC3 t){
  //compute the matrices
  MATRIX4 mvp = applyViewportMatrix();
  MATRIX4 m_cam = applyCamMatrix(e, g, t);
  MATRIX4 m_pers = applyPersMatrix(65.0, 4.0/3.0, 1.0, 100.0);

  //build the cube
  vector<VEC3> vertices;
  vector<VEC3I> indices;
  vector<VEC3> colors;
  vector<VEC4> new_vertices;
  buildCubePerVertexColors(vertices, indices, colors);

  for(auto vertex : vertices){
    VEC4 new_vertex = extend(vertex);
    new_vertex = mvp * m_pers * m_cam * new_vertex;

    if (new_vertex[3] != 0.0) {
      new_vertex[0] /= new_vertex[3];
      new_vertex[1] /= new_vertex[3];
      new_vertex[2] /= new_vertex[3];
    }

    new_vertices.push_back(new_vertex);
  }

  string filename = "custom.ppm";
  int totalCells = xScreenRes * yScreenRes;
  float* values = new float[3 * totalCells];

  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = 0.0f;

  vector<double> depth_storage (totalCells);

  for (int i = 0; i < totalCells; i++)
    depth_storage[i] = numeric_limits<double>::max();

  for(int i = 0; i < 12; i++){
    color_interp(values, new_vertices, indices, colors, i, depth_storage);
  }

  writePPM(filename, xScreenRes, yScreenRes, values);
}

int main(int argc, char** argv)
{
  if(argc == 1){
    make_2_ppm();
    make_3_ppm();
    make_4_ppm();
    make_5_ppm();
    make_6_ppm();
    make_7_ppm();
    make_8_ppm();
  }
  if(argc == 10){
    float eye_x = stof(argv[1]);
    float eye_y = stof(argv[2]);
    float eye_z = stof(argv[3]);
    float lookat_x = stof(argv[4]);
    float lookat_y = stof(argv[5]);
    float lookat_z = stof(argv[6]);
    float up_x = stof(argv[7]);
    float up_y = stof(argv[8]);
    float up_z = stof(argv[9]);
    VEC3 e(eye_x, eye_y, eye_z);
    VEC3 target(lookat_x, lookat_y, lookat_z);
    VEC3 g = target - e;
    VEC3 t(up_x, up_y, up_z);
    make_custom_ppm(e, g, t);
  }
  return 0;
}
