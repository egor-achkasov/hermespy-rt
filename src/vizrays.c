#include "../inc/vec3.h" /* for Vec3 */
#include "../inc/scene.h" /* for Scene, Mesh, scene_load */
#include "../inc/ray.h" /* for Ray */

#include "../test/test.h" /* for g_numPaths, g_numBounces, g_numTx, g_numRx */

#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <ctype.h> /* for tolower */ 

/**
 * Structs
 */

typedef struct {
  Vec3 o, t;
  float d;
  float lastX, lastY;
  float yaw, pitch, roll;
} Camera;

/**
 * Globals
 */

Camera g_cam = {
  {0.f, 0.f, 5.f},
  {0.f, 0.f, 0.f},
  5.f,
  0.f, 0.f,
  0.f, 0.f, 0.f
};
uint8_t g_mouseDown = 0;

Ray* g_rays = NULL;
uint8_t *g_active = NULL;

Scene g_scene = {0};

uint32_t g_bounce_cur = 0;

/**
 * Loading functions
 */

void loadRays(const char* filename) {
  FILE* file = fopen(filename, "rb");
  if (!file) {
    printf("Failed to open file\n");
    exit(8);
  }
  
  uint32_t fileSize = (g_numBounces + 1) * g_numTx * g_numPaths * 2 * 3;
  g_rays = (Ray*)malloc(fileSize * sizeof(float));
  if (fread(g_rays, sizeof(float), fileSize, file) != fileSize) {
    printf("Failed to read rays\n");
    exit(8);
  }
  fclose(file);

  /* active.bin */
  file = fopen("active.bin", "rb");
  if (!file) {
    printf("Failed to open file\n");
    exit(8);
  }
  uint32_t activeSize = (g_numBounces + 1) * (g_numTx * g_numPaths / 8 + 1);
  g_active = (uint8_t*)malloc(activeSize * sizeof(uint8_t));
  if (fread(g_active, sizeof(uint8_t), activeSize, file) != activeSize) {
    printf("Failed to read active\n");
    exit(8);
  }
  fclose(file);
}

/**
 * Draw functions
 */

void drawScene() {
  glBegin(GL_TRIANGLES);
  for (uint32_t i = 0; i < g_scene.num_meshes; i++) {
    Mesh* mesh = &g_scene.meshes[i];
    float mesh_color = (float)i / g_scene.num_meshes;
    glColor3f(mesh_color, 1.f - mesh_color, 0.f);

    for (uint32_t j = 0; j < mesh->num_triangles; j++) {
      uint32_t idx = mesh->is[j * 3];
      Vec3 v1 = mesh->vs[idx];
      idx = mesh->is[j * 3 + 1];
      Vec3 v2 = mesh->vs[idx];
      idx = mesh->is[j * 3 + 2];
      Vec3 v3 = mesh->vs[idx];
      glVertex3f(v1.x, v1.y, v1.z);
      glVertex3f(v2.x, v2.y, v2.z);
      glVertex3f(v3.x, v3.y, v3.z);
    }
  }
  glEnd();
}

void drawRays() {
  float color = (float)g_bounce_cur / g_numBounces;
  glColor3f(color, 0.0f, 1.0f - color);
  
  uint8_t active_bit = 1;
  uint32_t num_skipped = 0;
  for (uint32_t tx = 0; tx < g_numTx; ++tx) {
    glBegin(GL_LINES);
 
    uint32_t active_byte = (tx * (g_numBounces + 1) + g_bounce_cur) * (g_numPaths / 8 + 1);
    uint32_t ray_ind_base = (tx * (g_numBounces + 1) + g_bounce_cur) * g_numPaths;

    for (uint32_t path = 0; path < g_numPaths; ++path, active_bit <<= 1) {
      if (!active_bit) {
        active_byte++;
        active_bit = 1;
      }
      if (!(g_active[active_byte] & active_bit))
      {
        num_skipped++;
        continue;
      }
      
      uint32_t idx = ray_ind_base + path;
      GLfloat x = g_rays[idx].o.x;
      GLfloat y = g_rays[idx].o.y;
      GLfloat z = g_rays[idx].o.z;
      /* Origin */
      glVertex3f(x, y, z);
      /* Destination */
      glVertex3f(
        x + g_rays[idx].d.x,
        y + g_rays[idx].d.y,
        z + g_rays[idx].d.z
      );
    }
    glEnd();
    
    /* Draw points */
    glPointSize(5.0f);
    glBegin(GL_POINTS);
    active_byte = g_bounce_cur * (g_numTx * g_numPaths / 8 + 1);
    active_bit = 1;
    for (uint32_t path = 0; path < g_numPaths; ++path) {
      if (!active_bit) {
        active_byte++;
        active_bit = 1;
      }
      if (!(g_active[active_byte] & active_bit))
        continue;
      uint32_t idx = ray_ind_base + path;
      glVertex3f(g_rays[idx].o.x, g_rays[idx].o.y, g_rays[idx].o.z);
    }
    glEnd();
  }
  
  printf("\rSkipped %d paths", num_skipped);
  fflush(stdout);
}

/**
 * Callbacks
 */

void reshape(int w, int h) {
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0f, (float)w/h, 0.1f, 1000.0f);
  glMatrixMode(GL_MODELVIEW);
}

void display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  
  /* Look at target from camera position */
  gluLookAt(
      g_cam.o.x, g_cam.o.y, g_cam.o.z,
      g_cam.t.x, g_cam.t.y, g_cam.t.z,
      0.0f, 1.0f, 0.0f
  );

  drawScene();
  drawRays();
  
  glutSwapBuffers();
}

void updateCamera() {
  float x = g_cam.d * cosf(g_cam.pitch) * cosf(g_cam.yaw);
  float y = g_cam.d * sinf(g_cam.pitch);
  float z = g_cam.d * cosf(g_cam.pitch) * sinf(g_cam.yaw);
  
  g_cam.o.x = g_cam.t.x + x;
  g_cam.o.y = g_cam.t.y + y;
  g_cam.o.z = g_cam.t.z + z;

  gluLookAt(g_cam.o.x, g_cam.o.y, g_cam.o.z,
            g_cam.t.x, g_cam.t.y, g_cam.t.z,
            0.f, 1.f, 0.f);
}

void keyboard(unsigned char key, int x, int y) {
  float moveSpeed = 0.1f;
  float tiltSpeed = 0.15f;
  /* Speed up if SHIFT is down */
  if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
    moveSpeed *= 5;
    key = tolower(key);
  }
  
  /* Calculate forward direction */
  Vec3 dir = {
    g_cam.t.x - g_cam.o.x,
    g_cam.t.y - g_cam.o.y,
    g_cam.t.z - g_cam.o.z
  };
  float norm = sqrtf(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
  Vec3 forward = {dir.x / norm, dir.y / norm, dir.z / norm};
  
  /* Calculate right vector (cross product of forward and up (0,0,1)) */
  norm = sqrtf(forward.z * forward.z + forward.x * forward.x);
  Vec3 right = {forward.z / norm, 0.0f, -forward.x / norm};
  
  switch (key) {
    case 'w': /* Move target forward */
      g_cam.t.x += forward.x * moveSpeed;
      g_cam.t.z += forward.z * moveSpeed;
      break;
    case 's': /* Move target backward */
      g_cam.t.x -= forward.x * moveSpeed;
      g_cam.t.z -= forward.z * moveSpeed;
      break;
    case 'a': /* Move target left */
      g_cam.t.x += right.x * moveSpeed;
      g_cam.t.z += right.z * moveSpeed;
      break;
    case 'd': /* Move target right */
      g_cam.t.x -= right.x * moveSpeed;
      g_cam.t.z -= right.z * moveSpeed;
      break;

    case 'q': /* Roll left */
      g_cam.roll += tiltSpeed;
      break;
    case 'e': /* Roll right */
      g_cam.roll -= tiltSpeed;
      break;

    case 'x': /* Increase g_curBounce */
      if (g_bounce_cur < g_numBounces) g_bounce_cur++;
      break;
    case 'z': /* Decrease g_curBounce */
      if (g_bounce_cur > 0) g_bounce_cur--;
      break;
  }
  
  updateCamera();
  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {
  if (button == GLUT_LEFT_BUTTON) {
    g_mouseDown = (state == GLUT_DOWN);
    g_cam.lastX = x;
    g_cam.lastY = y;
  }
  /* Handle scroll wheel */
  else if (button == 3) { /* Scroll up */
    g_cam.d = fmaxf(0.1f, g_cam.d - 0.5f); /* Prevent going too close */
    updateCamera();
    glutPostRedisplay();
  }
  else if (button == 4) { /* Scroll down */
    g_cam.d += 0.5f;
    updateCamera();
    glutPostRedisplay();
  }
}

void mouseMotion(int x, int y) {
  if (!g_mouseDown) return;

  float sensitivity = 0.005f;
  float dx = (float)(x - g_cam.lastX);
  float dy = (float)(y - g_cam.lastY);
  
  g_cam.yaw -= dx * sensitivity;
  g_cam.pitch -= dy * sensitivity;
  
  /* Clamp pitch to prevent flipping */
  if (g_cam.pitch > 1.5f) g_cam.pitch = 1.5f;
  if (g_cam.pitch < -1.5f) g_cam.pitch = -1.5f;
  
  g_cam.lastX = x;
  g_cam.lastY = y;
  
  updateCamera();
  glutPostRedisplay();
}

/**
 * Main
 */

int main(int argc, char** argv) {
  const char* rays_filename = "rays.bin";
  
  /* Load ray data */
  loadRays(rays_filename);
  if (argc > 1)
    g_scene = scene_load(argv[1]);
  
  /* Initialize GLUT */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(800, 600);
  glutCreateWindow("Ray Visualization");
  
  /* Set up OpenGL */
  glEnable(GL_DEPTH_TEST);
  glClearColor(0.f, 0.f, 0.f, 1.0f);
  
  /* Register callbacks */
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(mouseMotion);
  
  updateCamera();
  glutMainLoop();
  
  return 0;
}
