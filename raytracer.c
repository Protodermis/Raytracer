#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "vecmat.h"
#include "msg.h"
#include <string.h>
#include <math.h>

/*
A raytracer program which takes the arguments "reference"
or "custom" to determine which of two predefined scenes to
genreate a .png image of.

Author: Zane Barker

*/


/* Information on the perspective through which the image is drawn 
   width is also height as viewport and final image are square */
typedef struct Perspective {
  float camPos[3];  // location of camera in 3D space
  float screenDist; // The distance from the camera(screen) to the viewport 
  int   worldWidth; // The width of the viewport in world coordinates through which the perspective is drawn.
  int   screenWidth; // The width in pixels of the screen
} Perspective;


typedef struct Ray {
  float startPos[3]; // 3D point of origin
  float vec[3];      // Unit vector for the direction of ray in 3D space
} Ray;

/* Defines surface characteristics of 
   a triangle or sphere */
typedef struct Material {
  int refl;        // A value of one means the material is reflective
  float color[3];  // RGB value. Each value is scaled from 0 to 1
} Material;

typedef struct Sphere {
  float pos[3]; // Center point of sphere in 3D space
  float rad;    // Sphere radius
  Material* mat;
} Sphere;

typedef struct Triangle {
  float vert1[3];  // The points of a triangle's three vertices in 3D space
  float vert2[3];
  float vert3[3];
  Material* mat;
} Triangle;

/* Contains information of the surface a ray intersects with
   as well as the incoming ray. */
typedef struct RayHit {
  float t;         /* The 'time' taken for the ray to reach the hit surface from its
		      orgin point */
  float normal[3]; // Normal vector of the surface hit by ray
  int   isRefl;    // Reflection value of the hit surface
  float loc[3];    // Ray-surface intersection location in 3D space
  Ray* inRay;      // The ray which hit the surface
  float color[3];  // The color value of the hit surface
} RayHit;

Sphere* spheres;
Triangle* triangles;
Material* materials;
float lightPos[3];    // 3D light source location in 3D space
int   numSpheres;
int   numTriangles;
int   numMat;

/* Returns the index in the image array of a pixel's
   RGB value.

   @param row  - the row of the pixel
   @param col  - the column of the pixel
   @param comp - the R, G, or B value of the pixel
                 R = 0
                 G = 1
                 B = 2                                */
int getArrayIndex(int row, int col, int comp) {
  return col*3 + row*1536 + comp;
}


/* 'Builds' the reference scene
   as defined by the instrucotr */
void setSceneRef(void) {
  numMat = 4;
  numSpheres = 3;
  numTriangles = 5;

  // Set light source
  vec3f_set(lightPos, 3,5,-15);

  /* Set up materials */
  materials = malloc(sizeof(*materials) * numMat);
  
  // Reflective material
  materials[0].refl = 1;
  vec3f_set(materials[0].color, 0,0,0);
  
  // Red
  materials[1].refl = 0;
  vec3f_set(materials[1].color, 1,0,0);

  // Blue
  materials[2].refl = 0;
  vec3f_set(materials[2].color, 0,0,1);

  // White
  materials[3].refl = 0;
  vec3f_set(materials[3].color, 1,1,1);


  /* Set up spheres */
  spheres = malloc(sizeof(*spheres) * numSpheres);

  // Sphere 1
  vec3f_set(spheres[0].pos, 0,0,-16);
  spheres[0].rad = 2;
  spheres[0].mat = &materials[0]; // reflective

  // Sphere 2
  vec3f_set(spheres[1].pos, 3,-1,-14);
  spheres[1].rad = 1;
  spheres[1].mat = &materials[0]; // reflective

  // Sphere 3
  vec3f_set(spheres[2].pos, -3,-1,-14);
  spheres[2].rad = 1;
  spheres[2].mat = &materials[1]; // red
  
  /* Set up triangles */
  triangles = malloc(sizeof(*triangles) * numTriangles);
  
  // Back wall triangles
  vec3f_set(triangles[0].vert1, -8,-2,-20);
  vec3f_set(triangles[0].vert2,  8,-2,-20);
  vec3f_set(triangles[0].vert3,  8,10,-20);
  triangles[0].mat = &materials[2];  // blue
  vec3f_set(triangles[1].vert1, -8,-2,-20);
  vec3f_set(triangles[1].vert2,  8,10,-20);
  vec3f_set(triangles[1].vert3, -8,10,-20);
  triangles[1].mat = &materials[2];  // blue

  // Floor triangles
  vec3f_set(triangles[2].vert1, -8,-2,-20);
  vec3f_set(triangles[2].vert2,  8,-2,-10);
  vec3f_set(triangles[2].vert3,  8,-2,-20);
  triangles[2].mat = &materials[3];  // white
  vec3f_set(triangles[3].vert1, -8,-2,-20);
  vec3f_set(triangles[3].vert2, -8,-2,-10);
  vec3f_set(triangles[3].vert3,  8,-2,-10);
  triangles[3].mat = &materials[3];  // white

  // Right triangle
  vec3f_set(triangles[4].vert1,  8,-2,-20);
  vec3f_set(triangles[4].vert2,  8,-2,-10);
  vec3f_set(triangles[4].vert3,  8,10,-20);
  triangles[4].mat = &materials[1];  // red

}

/* Sets up the custom scene */
void setSceneCustom(void) {
  numMat = 6;
  numSpheres = 4;
  numTriangles = 5;

  // Set light source
  vec3f_set(lightPos, 0,8,-6);

  /* Set up materials */
  materials = malloc(sizeof(*materials) * numMat);
  
  // Reflective material
  materials[0].refl = 1;
  vec3f_set(materials[0].color, 0,0,0);
  
  // Green
  materials[1].refl = 0;
  vec3f_set(materials[1].color, 0.13,0.55,0.13);

  // Brown
  materials[2].refl = 0;
  vec3f_set(materials[2].color, 0.66,0.16,0.16);

  // Yellow
  materials[3].refl = 0;
  vec3f_set(materials[3].color, 1,1,0);

  // Red
  materials[4].refl = 0;
  vec3f_set(materials[4].color, 1,0,0);

  // Blue
  materials[5].refl = 0;
  vec3f_set(materials[5].color, 0,0,1);


  /* Set up spheres */
  spheres = malloc(sizeof(*spheres) * numSpheres);

  // Yellow Sphere
  vec3f_set(spheres[0].pos, 0,3.25,-20);
  spheres[0].rad = 0.5;
  spheres[0].mat = &materials[3]; // yellow

  // Other ornaments
  vec3f_set(spheres[1].pos,  1, 0.5,-17);
  spheres[1].rad = 0.5;
  spheres[1].mat = &materials[0]; // reflective
  vec3f_set(spheres[2].pos, -1.5,-1.5,-17);
  spheres[2].rad = 0.5;
  spheres[2].mat = &materials[4]; // red
  vec3f_set(spheres[3].pos, 1.75,-3.5,-17);
  spheres[3].rad = 0.5;
  spheres[3].mat = &materials[5]; // blue
  
  /* Set up triangles */
  triangles = malloc(sizeof(*triangles) * numTriangles);
  
  // Tree triangles
  vec3f_set(triangles[0].vert1,    0,3,-20);
  vec3f_set(triangles[0].vert2, -2.5,0,-17);
  vec3f_set(triangles[0].vert3,  2.5,0,-17);
  triangles[0].mat = &materials[1];  // Green
  vec3f_set(triangles[1].vert1,   0, 1,-20);
  vec3f_set(triangles[1].vert2,-3.5,-2,-17);
  vec3f_set(triangles[1].vert3, 3.5,-2,-17);
  triangles[1].mat = &materials[1];  // Green
  vec3f_set(triangles[2].vert1,    0,-1,-20);
  vec3f_set(triangles[2].vert2, -4.5,-4,-17);
  vec3f_set(triangles[2].vert3,  4.5,-4,-17);
  triangles[2].mat = &materials[1];  // Green
  
  // Trunk triangles
  vec3f_set(triangles[3].vert1, -1,-3,-18);
  vec3f_set(triangles[3].vert2,  1,-3,-18);
  vec3f_set(triangles[3].vert3, -1,-10,-18);
  triangles[3].mat = &materials[2];  // brown
  vec3f_set(triangles[4].vert1,  1,-10,-18);
  vec3f_set(triangles[4].vert2, -1,-10,-18);
  vec3f_set(triangles[4].vert3,  1,-3,-18);
  triangles[4].mat = &materials[2];  // brown

}

/* 
   Checks if a given ray intersects with a given sphere
   in the scene.

   @return A RayHit containing details of the hit location.
           If the Ray does not interesect, the returned RayHit's
	   t value is -1.
*/
RayHit sphereIntersect(Ray * r, Sphere* sph) {
  RayHit hit;
  hit.inRay  = r;
  
  // Set variables from ray and sphere
  float e[3], d[3], c[3];
  vec3f_copy(e, r->startPos);
  vec3f_copy(d, r->vec);

  vec3f_copy(c, sph->pos); // center of sphere
  float rad = sph->rad;

  float eminc[3], ddotd; // (e-c) and d*d
  vec3f_sub_new(eminc, e, c);
  ddotd = vec3f_dot(d,d);
  
  // Find discriminant
  float discrm = pow(vec3f_dot(d,eminc),2) - ddotd*(vec3f_dot(eminc,eminc)-rad*rad);
  float negd[3];
  if(discrm < 0) {
    hit.t=-1;   // Ray does not intersect
  } else {
    vec3f_scalarMult_new(negd, d, -1);
    // Find time
    hit.t = fmin((vec3f_dot(negd, eminc) + sqrt(discrm))/ddotd,
		  (vec3f_dot(negd, eminc) - sqrt(discrm))/ddotd);
    // Add refl value of geometry material
    hit.isRefl = sph->mat->refl;
    
    // Hit location
    float hitVec[3];
    vec3f_scalarMult_new(hitVec, d, hit.t);
    vec3f_add_new(hit.loc,r->startPos, hitVec);

    // Find Normal vec
    vec3f_sub_new(hit.normal, hit.loc, sph->pos);
    vec3f_normalize_new(hit.normal, hit.normal);
    
    // Get color
    vec3f_copy(hit.color, sph->mat->color);
  }

  return hit;
}

/* 
   Checks to see if a ray hits a triangle in the scene 

   @return A RayHit containing details of the hit location.
           If the Ray does not interesect, the returned RayHit's
           t value is negative.
*/
RayHit triangleIntersect(Ray * r, Triangle* tri) {
  RayHit hit;
  hit.inRay = r;
  
  /* Set values for the ray-triangle intersection equations
     provided by the instructor */
  float A,B,C,D,E,F,G,H,I,J,K,L,M,beta,gamma;
  A = tri->vert1[0]-tri->vert2[0]; // xa-xb
  B = tri->vert1[1]-tri->vert2[1]; // ya-yb
  C = tri->vert1[2]-tri->vert2[2]; // za-zb
  D = tri->vert1[0]-tri->vert3[0]; // xa-xc
  E = tri->vert1[1]-tri->vert3[1]; // ya-yc
  F = tri->vert1[2]-tri->vert3[2]; // za-zc
  G = r->vec[0]; // xd
  H = r->vec[1]; // yd
  I = r->vec[2]; // zd
  J = tri->vert1[0] - r->startPos[0]; // xa-xe
  K = tri->vert1[1] - r->startPos[1]; // ya-ye
  L = tri->vert1[2] - r->startPos[2]; // za-ze

  M = A*(E*I-H*F)+B*(G*F-D*I)+C*(D*H-E*G);

  hit.t = -1*(F*(A*K-J*B)+E*(J*C-A*L)+D*(B*L-K*C))/M;
  if(hit.t < 0)
    return hit; // ray does not interesect

  gamma = (I*(A*K-J*B)+H*(J*C-A*L)+G*(B*L-K*C))/M;
  if(gamma < 0 || gamma > 1) {
    hit.t = -1; // ray does not intersect
    return hit;
  }
  
  beta = (J*(E*I-H*F)+K*(G*F-D*I)+L*(D*H-E*G))/M;
  if(beta<0 || beta > 1-gamma) {
    hit.t = -1; // ray does not intersect
    return hit;
  }
  
  // Find if triangle is reflective
  hit.isRefl = tri->mat->refl;
  
  // Get hit location
  float hitVec[3];
  vec3f_scalarMult_new(hitVec, r->vec, hit.t);
  vec3f_add_new(hit.loc,r->startPos, hitVec);

  // Get color
  vec3f_copy(hit.color, tri->mat->color);
  
  // Find Normal
  float vecAB[3], vecAC[3];
  vec3f_sub_new(vecAB, tri->vert2, tri->vert1);
  vec3f_sub_new(vecAC, tri->vert3, tri->vert1);
  vec3f_cross_new(hit.normal, vecAB, vecAC);
  vec3f_normalize_new(hit.normal, hit.normal);
  
  return hit;
}

/* 
   Calculates a diffuse shading value on a hit surface

   @return The minimum value of 0.2 is returned if the point on the surface
   is completely obstructed by another object in the scene. Otherwise a greater
   value is calculated based on the surfaces proximity to the light source.
*/
float diffuseShading(RayHit pixHit) {
  
  // Get shadow vector
  Ray shadRay;
  vec3f_sub_new(shadRay.vec, lightPos, pixHit.loc);
  vec3f_normalize_new(shadRay.vec, shadRay.vec);
  
  // Bump starting location
  float bump[3];
  vec3f_scalarMult_new(bump, shadRay.vec, 0.0001);
  vec3f_add_new(shadRay.startPos, pixHit.loc, bump);
  
  // Check if the shadow ray intersects with any object
  // If so, return a value of .2 to indicate the pixel is in shadow
  RayHit shadHit;
  shadHit.t = -1;
  for(int i=0;i<numSpheres;i++) {
    shadHit = sphereIntersect(&shadRay, &spheres[i]);
    if(shadHit.t>0)
      return 0.2;
  }

  for(int i=0;i<numTriangles;i++) {
    shadHit = triangleIntersect(&shadRay, &triangles[i]);
    if(shadHit.t>0)
      return 0.2;
  }
 
  // If the pixel is unobstructed from the light calculate
  // the diffuse shading value
  return fmax(0.2,vec3f_dot(pixHit.normal, shadRay.vec));    
}

/* 
   Calculate the ultimate color of a reflective surface 
   using relfection rays
   
   @return a RayHit whose color is the value returned by the reflected
           rays. If more then 10 reflections occur, the color is black.
*/
RayHit reflect(RayHit pixHit, int depth){
  RayHit reflHit;
  if(depth == 9) {
    vec3f_set(reflHit.color, 0, 0, 0);
    return reflHit;
  }
  Ray r; // Reflect ray
  float ddotn = vec3f_dot(pixHit.inRay->vec, pixHit.normal);
  // Calculate r
  float n[3];
  vec3f_scalarMult_new(n, pixHit.normal, 2*ddotn);
  vec3f_sub_new(r.vec, pixHit.inRay->vec, n);
  vec3f_normalize_new(r.vec, r.vec);
  
  // Bump starting location
  float bump[3];
  vec3f_scalarMult_new(bump, r.vec, 0.0001);
  vec3f_add_new(r.startPos, pixHit.loc, bump);
  
  // Find nearest reflHit
  RayHit checkHit;
  reflHit.t=100000;
  reflHit.isRefl = -1;
  
  // Check for intersections with spheres
  for(int i=0; i<numSpheres; i++){
    checkHit = sphereIntersect(&r, &spheres[i]);
    if(checkHit.t > 0 && checkHit.t < reflHit.t){
      reflHit = checkHit;
    }
  }
  
  // Check for intersections with triangles
  for(int i=0; i<numTriangles; i++){
    checkHit = triangleIntersect(&r, &triangles[i]);
    if(checkHit.t > 0 && checkHit.t < reflHit.t) {
      reflHit = checkHit;
    }
  }
  
  if(reflHit.isRefl==1)
    return reflect(reflHit, depth+1);
  else if(reflHit.isRefl==0) {
    float shading = diffuseShading(reflHit);
    vec3f_scalarMult_new(reflHit.color, reflHit.color, shading);
    return reflHit;
  } else {
    vec3f_set(reflHit.color, 0, 0, 0);
    return reflHit;
  }
  
}

/* Gets a ray from the camera to a pixel */
Ray getRay(Perspective p, int screenCoord[2]) {
  float ratio  = (float) p.worldWidth/p.screenWidth;
  float worldX = -1 * (p.worldWidth/2) + ratio*screenCoord[0]+ratio/2;
  float worldY = (p.worldWidth/2) - ratio*screenCoord[1] - ratio/2;

  // 3D Pixel position
  float pixPos[3];
  vec3f_set(pixPos, worldX , worldY,-2);
  float vector[3];
  vec3f_sub_new(vector, pixPos, p.camPos);
  
  Ray r;
  vec3f_normalize_new(r.vec, vector);
  vec3f_copy(r.startPos, p.camPos);
  return r;
}

int main(int argc, char** argv) {

  // contains RGB values for each pixel in output, in this case 512*512*3.
  unsigned char arrayContainingImage[786432];
  char * filename;
  
  if(argc == 1) {
    printf("ERROR: This program reqire input argument 'reference' or 'custom' to select the scene to render\n");
    exit(1);
  }
  
  // Determine which scene to draw
  if(strcmp(argv[1], "reference")==0)
    setSceneRef();
  else if(strcmp(argv[1], "custom")==0)
    setSceneCustom();
  else {
    printf("ERROR: This program reqire input argument 'reference' or 'custom' to select the scene to render\n");
    exit(1);
  }
  
  // Create perpective
  Perspective camera;
  vec3f_set(camera.camPos, 0,0,0);
  camera.screenDist  = -2;
  camera.worldWidth  = 2;
  camera.screenWidth = 512;


  // append the file extension
  filename = strcat(argv[1], ".png");

  Ray ray;
  RayHit nearestHit, checkHit;
  // Find color for each pixel on the screen
  for(int x=0; x<512; x++) {
    for(int y=0; y<512; y++) {
      // Get Ray
      int pixel[]={x,y};
      ray = getRay(camera, pixel);
      
      // Init nearest hit
      nearestHit.t=100000;
      nearestHit.isRefl = -1;
     
      // Check for intersections with spheres
      for(int i=0; i<numSpheres; i++){
	checkHit = sphereIntersect(&ray, &spheres[i]);
	if(checkHit.t > 0 && checkHit.t < nearestHit.t){
	  nearestHit = checkHit;
	}
      }

      // Check for intersections with triangles
      for(int i=0; i<numTriangles; i++){
	checkHit = triangleIntersect(&ray, &triangles[i]);
	if(checkHit.t > 0 && checkHit.t < nearestHit.t) {
	  nearestHit = checkHit;
	}
      }

      float shading;
      // Find shadow and reflection values
      if(nearestHit.isRefl==0){
	shading = diffuseShading(nearestHit);
      } else if(nearestHit.isRefl==1){
	nearestHit = reflect(nearestHit, 0);
      } else {
	nearestHit.color[0]=0;
	nearestHit.color[1]=0;
	nearestHit.color[2]=0;
      }
      // Set color of pixel
      arrayContainingImage[getArrayIndex(y,x,0)] = shading * nearestHit.color[0]*255;
      arrayContainingImage[getArrayIndex(y,x,1)] = shading * nearestHit.color[1]*255;
      arrayContainingImage[getArrayIndex(y,x,2)] = shading * nearestHit.color[2]*255;
    }
  }
  // write an image out to a PNG file
  stbi_write_png(filename, camera.screenWidth, camera.screenWidth, 3, arrayContainingImage, camera.screenWidth*3);

  free(materials);
  free(spheres);
  free(triangles);

}
