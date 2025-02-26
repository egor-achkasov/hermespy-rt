#include "test.h" /* for g_numPaths, g_numBounces, g_numTx, g_numRx */

#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

/**
 * Structs
 */

typedef struct {
  float x, y, z;
} Vec3;

typedef struct {
  Vec3 o;
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

Camera g_cam = {{0.f, 0.f, 5.f}, 0.f, 0.f, 0.f, 0.f};
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
  fread(g_rays, sizeof(float), fileSize, file);
  fclose(file);

  /* active.bin */
  file = fopen("active.bin", "rb");
  if (!file) {
    printf("Failed to open file\n");
    exit(8);
  }
  uint32_t activeSize = (g_numBounces + 1) * (g_numTx * g_numPaths / 8 + 1);
  g_active = (uint8_t*)malloc(activeSize * sizeof(uint8_t));
  fread(g_active, sizeof(uint8_t), activeSize, file);
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
  fread(magic, 1, 3, file);
  if (magic[0] != 'H' || magic[1] != 'R' || magic[2] != 'T') {
    printf("Invalid file format\n");
    exit(8);
  }

  /* Scene */
  fread(&g_scene.num_meshes, sizeof(uint32_t), 1, file);
  g_scene.meshes = (Mesh*)malloc(g_scene.num_meshes * sizeof(Mesh));

  /* Meshes */
  for (Mesh* mesh = g_scene.meshes;
       mesh != g_scene.meshes + g_scene.num_meshes;
       ++mesh)
  {
    fread(&mesh->num_vertices, sizeof(uint32_t), 1, file);
    mesh->vs = (Vec3*)malloc(mesh->num_vertices * sizeof(Vec3));
    fread(mesh->vs, sizeof(Vec3), mesh->num_vertices, file);
    fread(&mesh->num_triangles, sizeof(uint32_t), 1, file);
    mesh->is = (uint32_t*)malloc(mesh->num_triangles * 3 * sizeof(uint32_t));
    fread(mesh->is, sizeof(uint32_t), mesh->num_triangles * 3, file);
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
  printf("Skipped %d paths\n", num_skipped);
  
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

void display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  
  /* Apply camera transformation */
  glRotatef(g_cam.pitch, 1.0f, 0.0f, 0.0f);
  glRotatef(g_cam.yaw, 0.0f, 1.0f, 0.0f);
  glTranslatef(-g_cam.o.x, -g_cam.o.y, -g_cam.o.z);

  drawScene();
  drawRays();
  
  glutSwapBuffers();
}

void reshape(int w, int h) {
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0f, (float)w/h, 0.1f, 1000.0f);
  glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y) {
  float speed = 0.5f;
  switch (key) {
    case 'w': // Forward
      g_cam.o.x += sin(g_cam.yaw * M_PI/180.0f) * speed;
      g_cam.o.z -= cos(g_cam.yaw * M_PI/180.0f) * speed;
      break;
    case 's': // Backward
      g_cam.o.x -= sin(g_cam.yaw * M_PI/180.0f) * speed;
      g_cam.o.z += cos(g_cam.yaw * M_PI/180.0f) * speed;
      break;
    case 'a': // Left
      g_cam.o.x -= cos(g_cam.yaw * M_PI/180.0f) * speed;
      g_cam.o.z -= sin(g_cam.yaw * M_PI/180.0f) * speed;
      break;
    case 'd': // Right
      g_cam.o.x += cos(g_cam.yaw * M_PI/180.0f) * speed;
      g_cam.o.z += sin(g_cam.yaw * M_PI/180.0f) * speed;
      break;
    case 'q': // Up
      g_cam.o.y += speed;
      break;
    case 'e': // Down
      g_cam.o.y -= speed;
      break;
    case 'x': // Next bounce
      if (g_bounce_cur != g_numBounces)
        g_bounce_cur++;
      break;
    case 'z': // Previous bounce
      if (g_bounce_cur != 0)
        g_bounce_cur--;
      break;
    case 27: // Escape key
      free(g_rays);
      for (uint32_t i = 0; i < g_scene.num_meshes; i++) {
        free(g_scene.meshes[i].vs);
        free(g_scene.meshes[i].is);
      }
      free(g_scene.meshes);
      exit(0);
  }
  glutPostRedisplay();
}

void motion(int x, int y) {
  if (!g_mouseDown) return;
  g_cam.yaw += (x - g_cam.lastX) * 0.2f;
  g_cam.pitch += (y - g_cam.lastY) * 0.2f;
  g_cam.pitch = fmax(-89.0f, fmin(89.0f, g_cam.pitch));
  g_cam.lastX = x;
  g_cam.lastY = y;
  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {
  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_DOWN) {
      g_mouseDown = 1;
      g_cam.lastX = x;
      g_cam.lastY = y;
    } else {
      g_mouseDown = 0;
    }
  }
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
  glutMotionFunc(motion);
  
  glutMainLoop();
  
  return 0;
}
