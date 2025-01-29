#include "material.h" /* for Material */

#include <stdio.h>
#include <stdint.h> /* for uint8_t */

Material mats[2] = {
  {5, "metal", 1.f, 0.f, 10000000.f, 0.f, 0.f, 0},
  {8, "concrete", 5.24f, 0.f, 0.0462f, 0.7822f, 0.5f, 4}
};

/* Configure your scene's materials here */
#define NUM_MATERIALS 1 /* number of materials in the scene */
#define NUM_INDICES 6 /* number of material indices in the scene (= number of meshes) */
int material_indices[NUM_INDICES] = {1, 1, 1, 1, 1, 1};

int main() {
  FILE *fp = fopen("materials.bin", "wb");
  int num_materials = NUM_MATERIALS;
  int num_indices = NUM_INDICES;
  fwrite(&num_materials, sizeof(int), 1, fp);
  fwrite(&num_indices, sizeof(int), 1, fp);
  fwrite(material_indices, sizeof(int), NUM_INDICES, fp);
  for (int i = 0; i < NUM_MATERIALS; ++i) {
    fwrite(&mats[i].name_size, sizeof(int), 1, fp);
    fwrite(mats[i].name, sizeof(char), mats[i].name_size, fp);
    fwrite(&mats[i].a, sizeof(float), 1, fp);
    fwrite(&mats[i].b, sizeof(float), 1, fp);
    fwrite(&mats[i].c, sizeof(float), 1, fp);
    fwrite(&mats[i].d, sizeof(float), 1, fp);
    fwrite(&mats[i].s, sizeof(float), 1, fp);
    fwrite(&mats[i].alpha, sizeof(uint8_t), 1, fp);
  }
  fclose(fp);

  return 0;
}
