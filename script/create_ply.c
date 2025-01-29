/* vim: set tabstop=2:softtabstop=2:shiftwidth=2:noexpandtab */

#include "material.h" /* for Material */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct {
  unsigned int num_vertices;
  float *vertices;
  unsigned int num_faces;
  int *faces;
} Mesh;
  
int main(int argc, char *argv[]) {
  if (argc < 3) {
    printf("Usage: %s <materials.bin> <path_to_mesh1.ply> [path_to_mesh2.ply] ...\n", argv[0]);
    return 1;
  }

  /* materials */
  FILE *fp = fopen(argv[1], "rb");
  if (!fp) {
    printf("Cannot open file %s\n", argv[2]);
    return 1;
  }
  int num_materials;
  int num_indices;
  fread(&num_materials, sizeof(int), 1, fp);
  fread(&num_indices, sizeof(int), 1, fp);
  Material *materials = (Material *)malloc(sizeof(Material) * num_materials);
  int *material_indices = (int *)malloc(sizeof(int) * num_indices);
  fread(material_indices, sizeof(int), num_indices, fp);
  for (int i = 0; i < num_materials; ++i) {
    fread(&materials[i].name_size, sizeof(int), 1, fp);
    materials[i].name = (char *)malloc(sizeof(char) * materials[i].name_size);
    fread(materials[i].name, sizeof(char), materials[i].name_size, fp);
    fread(&materials[i].a, sizeof(float), 1, fp);
    fread(&materials[i].b, sizeof(float), 1, fp);
    fread(&materials[i].c, sizeof(float), 1, fp);
    fread(&materials[i].d, sizeof(float), 1, fp);
    fread(&materials[i].s, sizeof(float), 1, fp);
    fread(&materials[i].alpha, sizeof(uint8_t), 1, fp);
  }
  fclose(fp);

  /* meshes */
  char line[256];
  unsigned int num_vertices;
  unsigned int num_faces;
  unsigned int num_meshes = argc - 2;
  Mesh *meshes = (Mesh *)malloc(sizeof(Mesh) * num_meshes);
  for (int i = 2; i < argc; ++i) {
    Mesh *mesh = &meshes[i - 2];
    if ((fp = fopen(argv[i], "rb")) == NULL) {
      printf("Cannot open file %s\n", argv[1]);
      return 1;
    }
    fseek(fp, 36, SEEK_SET);
    fgets(line, 32, fp);
    sscanf(line, "element vertex %d\n", &num_vertices);
    fseek(fp, 85, SEEK_CUR);
    fgets(line, 32, fp);
    sscanf(line, "element face %d\n", &num_faces);
    fseek(fp, 39, SEEK_CUR);
    fgets(line, 12, fp);
    sscanf(line, "end_header\n");
    mesh->vertices = (float *)malloc(12 * num_vertices);
    mesh->faces = (int *)malloc(sizeof(int) * 3 * num_faces);
    mesh->num_vertices = num_vertices;
    mesh->num_faces = num_faces;
    for (int i = 0; i < num_vertices; ++i) {
      fread(&mesh->vertices[3 * i], 4, 3, fp);
      fseek(fp, 8, SEEK_CUR);
    }
    for (int i = 0; i < num_faces; ++i) {
      fseek(fp, 1, SEEK_CUR);
      fread(&mesh->faces[3 * i], 4, 3, fp);
    }
    fclose(fp);
  }

  /* unite meshes */
  num_vertices = 0;
  num_faces = 0;
  for (int i = 0; i < num_meshes; ++i) {
    num_vertices += meshes[i].num_vertices;
    num_faces += meshes[i].num_faces;
  }
  unsigned int sum = meshes[0].num_vertices;
  for (int i = 1; i < num_meshes; ++i) {
    for (unsigned int j = 0; j != meshes[i].num_faces; ++j)
      for (int k = 0; k < 3; ++k)
        meshes[i].faces[3 * j + k] += sum;
    sum += meshes[i].num_vertices;
  }

  fp = fopen("mesh.ply", "wb");
  fprintf(fp, "ply\n");
  fprintf(fp, "format binary_little_endian 1.0\n");
  fprintf(fp, "element vertex %d\n", num_vertices);
  fprintf(fp, "property float x\n");
  fprintf(fp, "property float y\n");
  fprintf(fp, "property float z\n");
  fprintf(fp, "element radio_material %d\n", num_materials);
  fprintf(fp, "property list uint char name\n");
  fprintf(fp, "property float a\n");
  fprintf(fp, "property float b\n");
  fprintf(fp, "property float c\n");
  fprintf(fp, "property float d\n");
  fprintf(fp, "property float s\n");
  fprintf(fp, "property char alpha\n");
  fprintf(fp, "element face %d\n", num_faces);
  fprintf(fp, "property list uchar uint vertex_index\n");
  fprintf(fp, "property uint material_index\n");
  fprintf(fp, "end_header\n");
  for (unsigned int j = 0; j != num_meshes; ++j)
    fwrite(meshes[j].vertices, 4, 3 * meshes[j].num_vertices, fp);
  for (int i = 0; i != num_materials; ++i) {
    fwrite(&materials[i].name_size, sizeof(int), 1, fp);
    fwrite(materials[i].name, sizeof(char), materials[i].name_size, fp);
    fwrite(&materials[i].a, sizeof(float), 1, fp);
    fwrite(&materials[i].b, sizeof(float), 1, fp);
    fwrite(&materials[i].c, sizeof(float), 1, fp);
    fwrite(&materials[i].d, sizeof(float), 1, fp);
    fwrite(&materials[i].s, sizeof(float), 1, fp);
    fwrite(&materials[i].alpha, sizeof(uint8_t), 1, fp);
  }
  char n = 3;
  for (unsigned int j = 0; j != num_meshes; ++j)
    for (unsigned int k = 0; k != meshes[j].num_faces; ++k) {
      fwrite(&n, 1, 1, fp);
      fwrite(&meshes[j].faces[3 * k], 4, 3, fp);
      fwrite(&material_indices[j], 4, 1, fp);
    }
  fclose(fp);
}
