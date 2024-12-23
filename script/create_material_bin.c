#include <stdio.h>

#define NUM_MATERIALS 1
#define NUM_INDICES 6

typedef struct {
  int name_size;
  char *name;
  float a, b, c, d;
} Material;

int main() {
  int material_indices[NUM_INDICES] = {0, 0, 0, 0, 0, 0};
  Material materials[NUM_MATERIALS] = {
    {5, "metal", 1., 0., 10000000., 0.}
  };

  FILE *fp = fopen("materials.bin", "wb");
  int num_materials = NUM_MATERIALS;
  int num_indices = NUM_INDICES;
  fwrite(&num_materials, sizeof(int), 1, fp);
  fwrite(&num_indices, sizeof(int), 1, fp);
  fwrite(material_indices, sizeof(int), NUM_INDICES, fp);
  for (int i = 0; i < NUM_MATERIALS; ++i) {
    fwrite(&materials[i].name_size, sizeof(int), 1, fp);
    fwrite(materials[i].name, sizeof(char), materials[i].name_size, fp);
    fwrite(&materials[i].a, sizeof(float), 1, fp);
    fwrite(&materials[i].b, sizeof(float), 1, fp);
    fwrite(&materials[i].c, sizeof(float), 1, fp);
    fwrite(&materials[i].d, sizeof(float), 1, fp);
  }
  fclose(fp);

  return 0;
}
