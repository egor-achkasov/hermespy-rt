#include "test.h" /* for g_numPaths, g_numBounces, g_numTx, g_numRx */

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
  float x, y, z;
} Vec3;

typedef struct {
  Vec3 o, t;
  float d;
  float lastX, lastY;
  float yaw, pitch;
} Camera;

typedef struct {
  uint32_t num_vertices;
  Vec3* vs;
  uint32_t num_triangles;
  uint32_t* is;
} Mesh;

typedef struct {
  uint32_t num_meshes;
  Mesh* meshes;
} Scene;

/**
 * Globals
 */

Camera g_cam = {
  {0.f, 0.f, 5.f},
  {0.f, 0.f, 0.f},
  5.f,
  0.f, 0.f,
  0.f, 0.f
};
uint8_t g_mouseDown = 0;

float* g_rays = NULL;
uint8_t *g_active = NULL;

Scene g_scene = {0};

int g_bounce_cur = 0;

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
  g_rays = (float*)malloc(fileSize * sizeof(float));
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

void loadScene(const char* filename) {
  FILE* file = fopen(filename, "rb");
  if (!file) {
    printf("Failed to open file\n");
    exit(8);
  }

  /* Magic */
  char magic[3];
  if (fread(magic, 1, 3, file) != 3) {
    printf("Failed to read magic\n");
    exit(8);
  }
  if (magic[0] != 'H' || magic[1] != 'R' || magic[2] != 'T') {
    printf("Invalid file format\n");
    exit(8);
  }

  /* Scene */
  if (fread(&g_scene.num_meshes, sizeof(uint32_t), 1, file) != 1) {
    printf("Failed to read num_meshes\n");
    exit(8);
  }
  g_scene.meshes = (Mesh*)malloc(g_scene.num_meshes * sizeof(Mesh));

  /* Meshes */
  for (Mesh* mesh = g_scene.meshes;
       mesh != g_scene.meshes + g_scene.num_meshes;
       ++mesh)
  {
    if (fread(&mesh->num_vertices, sizeof(uint32_t), 1, file) != 1) {
      printf("Failed to read num_vertices\n");
      exit(8);
    }
    mesh->vs = (Vec3*)malloc(mesh->num_vertices * sizeof(Vec3));
    if (fread(mesh->vs, sizeof(Vec3), mesh->num_vertices, file) != mesh->num_vertices) {
      printf("Failed to read vertices\n");
      exit(8);
    }
    if (fread(&mesh->num_triangles, sizeof(uint32_t), 1, file) != 1) {
      printf("Failed to read num_triangles\n");
      exit(8);
    }
    mesh->is = (uint32_t*)malloc(mesh->num_triangles * 3 * sizeof(uint32_t));
    if (fread(mesh->is, sizeof(uint32_t), mesh->num_triangles * 3, file) != mesh->num_triangles * 3) {
      printf("Failed to read indices\n");
      exit(8);
    }
    fseek(file, sizeof(uint32_t), SEEK_CUR); /* Skip material index */
  }
  fclose(file);
}

/**
 * Draw functions
 */

void drawScene() {
  glColor3f(1.0f, 1.0f, 1.0f);
  glBegin(GL_TRIANGLES);
  for (uint32_t i = 0; i < g_scene.num_meshes; i++) {
    Mesh* mesh = &g_scene.meshes[i];
    float mesh_color = (float)i / g_scene.num_meshes;
    glColor3f(0.f, mesh_color, 1.f - mesh_color);

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
  
  glBegin(GL_LINES);
  int idx_base = g_bounce_cur * g_numTx * g_numPaths * 2 * 3;
  uint32_t active_byte = g_bounce_cur * (g_numTx * g_numPaths / 8 + 1);
  uint8_t active_bit = 1;
  uint32_t num_skipped = 0;
  for (int path = 0; path < g_numPaths; ++path, active_bit <<= 1) {
    int idx = idx_base + path * 2 * 3;

    if (!active_bit) {
      active_byte++;
      active_bit = 1;
    }
    if (!(g_active[active_byte] & active_bit))
    {
      num_skipped++;
      continue;
    }
    
    /* Origin */
    float x1 = g_rays[idx];
    float y1 = g_rays[idx + 1];
    float z1 = g_rays[idx + 2];
    
    /* Destination */
    float x2 = x1 + g_rays[idx + 3];
    float y2 = y1 + g_rays[idx + 4];
    float z2 = z1 + g_rays[idx + 5];
    
    glVertex3f(x1, y1, z1);
    glVertex3f(x2, y2, z2);
  }
  glEnd();
  printf("\rSkipped %d paths", num_skipped);
  fflush(stdout);
  
  /* Draw points */
  glPointSize(5.0f);
  glBegin(GL_POINTS);
  active_byte = g_bounce_cur * (g_numTx * g_numPaths / 8 + 1);
  active_bit = 1;
  for (int path = 0; path < g_numPaths; ++path) {
    int idx = idx_base + path * 2 * 3;
    if (!active_bit) {
      active_byte++;
      active_bit = 1;
    }
    if (!(g_active[active_byte] & active_bit))
      continue;
    glVertex3f(g_rays[idx], g_rays[idx + 1], g_rays[idx + 2]);
  }
  glEnd();
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
            0.0f, 1.0f, 0.0f);
}

void keyboard(unsigned char key, int x, int y) {
  float speed = 0.1f;
  /* Speed up if SHIFT is down */
  if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
    speed *= 5;
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
      g_cam.t.x += forward.x * speed;
      g_cam.t.z += forward.z * speed;
      break;
    case 's': /* Move target backward */
      g_cam.t.x -= forward.x * speed;
      g_cam.t.z -= forward.z * speed;
      break;
    case 'a': /* Move target left */
      g_cam.t.x += right.x * speed;
      g_cam.t.z += right.z * speed;
      break;
    case 'd': /* Move target right */
      g_cam.t.x -= right.x * speed;
      g_cam.t.z -= right.z * speed;
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
    loadScene(argv[1]);
  
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
