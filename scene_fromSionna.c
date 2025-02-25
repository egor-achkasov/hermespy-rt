/* vim: set tabstop=2:softtabstop=2:shiftwidth=2:noexpandtab */

#include "scene.h" /* for HRT_Scene, HRT_Mesh, HRT_Material */

#include <stdio.h> /* for FILE, fopen, fclose, fseek, fgets, sscanf, printf */
#include <stdlib.h> /* for malloc */
#include <stdint.h> /* for uint8_t, uint32_t */

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("Usage: %s <scene.xml>\n", argv[0]);
    return 1;
  }

  /* TODO: parse the XML file */
  HRT_Mesh mesh_box;
  mesh_box.num_vertices = 8;
  float vs_temp[] = {
    5.f, 5.f, 0.f,
    -5.f, 5.f, 0.f,
    -5.f, -5.f, 0.f,
    5.f, -5.f, 0.f,
    5.f, 5.f, 5.f,
    -5.f, 5.f, 5.f,
    -5.f, -5.f, 5.f,
    5.f, -5.f, 5.f,
  };
  mesh_box.vs = (float*)vs_temp;
  mesh_box.num_triangles = 12;
  uint32_t is_temp[] = {
    0, 1, 2,
    0, 2, 3,
    0, 4, 5,
    0, 5, 1,
    1, 5, 6,
    1, 6, 2,
    2, 6, 7,
    2, 7, 3,
    3, 7, 4,
    3, 4, 0,
    4, 7, 6,
    4, 6, 5,
  };
  mesh_box.is = (uint32_t*)is_temp;
  mesh_box.material_index = MATERIAL_CONCRETE;

  HRT_Mesh mesh_simple_reflector;
  mesh_simple_reflector.num_vertices = 4;
  float vs_temp2[] = {
    -0.5f, -0.5f, 0.f,
    0.5f, -0.5f, 0.f,
    0.5f, 0.5f, 0.f,
    -0.5f, 0.5f, 0.f,
  };
  mesh_simple_reflector.vs = (float*)vs_temp2;
  mesh_simple_reflector.num_triangles = 2;
  uint32_t is_temp2[] = {
    0, 1, 2,
    0, 2, 3,
  };
  mesh_simple_reflector.is = (uint32_t*)is_temp2;
  mesh_simple_reflector.material_index = MATERIAL_CONCRETE;

  HRT_Scene scene;
  scene.num_meshes = 1;
  scene.meshes = &mesh_simple_reflector;

  FILE *fp = fopen("scene.hrt", "wb");
  if (fp == NULL) {
    printf("Error: cannot open file\n");
    return 1;
  }

  /* MAGIC */
  fwrite("HRT", 1, 3, fp);
  /* SCENE */
  /* num_meshes */
  fwrite(&scene.num_meshes, sizeof(uint32_t), 1, fp);
  /* meshes */
  for (uint32_t i = 0; i != scene.num_meshes; ++i) {
    /* MESH */
    fwrite(&scene.meshes[i].num_vertices, sizeof(uint32_t), 1, fp);
    fwrite(scene.meshes[i].vs, sizeof(float), 3 * scene.meshes[i].num_vertices, fp);
    fwrite(&scene.meshes[i].num_triangles, sizeof(uint32_t), 1, fp);
    fwrite(scene.meshes[i].is, sizeof(uint32_t), 3 * scene.meshes[i].num_triangles, fp);
    fwrite(&scene.meshes[i].material_index, sizeof(uint32_t), 1, fp);
  }

  fclose(fp);

  return 0;
}
