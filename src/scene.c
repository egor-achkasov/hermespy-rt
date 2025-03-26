#include "../inc/scene.h" /* for Scene, Mesh, Material */
#include "../inc/materials.h" /* for g_hrt_materials */
#include "../inc/common.h" /* for PERROT_CLEANUP_EXIT */

#include <stdio.h>  /* for fopen, FILE, fclose */

void scene_save(
  IN Scene* scene,
  IN const char* filepath
)
{
  FILE *fp = fopen(filepath, "wb");
  if (!fp)
    PERROR_CLEANUP_EXIT("Error: cannot open file", 8);

  /* MAGIC */
  fwrite("HRT", 1, 3, fp);
  /* SCENE */
  /* num_meshes */
  fwrite(&scene->num_meshes, sizeof(uint32_t), 1, fp);
  /* meshes */
  for (uint32_t i = 0; i != scene->num_meshes; ++i) {
    /* MESH */
    Mesh *mesh = &scene->meshes[i];
    fwrite(&mesh->num_vertices, sizeof(uint32_t), 1, fp);
    fwrite(mesh->vs, sizeof(Vec3), mesh->num_vertices, fp);
    fwrite(&mesh->num_triangles, sizeof(uint32_t), 1, fp);
    fwrite(mesh->is, sizeof(uint32_t), 3 * mesh->num_triangles, fp);
    fwrite(&mesh->material_index, sizeof(uint32_t), 1, fp);
    fwrite(&mesh->velocity, sizeof(Vec3), 1, fp);
  }

  fclose(fp);
}

Scene scene_load(IN const char *filepath)
{
  /* Open the HRT file */
  FILE *f = fopen(filepath, "rb");
  if (!f)
    PERROR_CLEANUP_EXIT("Could not open scene file", 8);

  /* Parse the file */
  /* MAGIC */
  char magic[3];
  if (fread(magic, 1, 3, f) != 3) exit(8);
  if (strncmp(magic, "HRT", 3)) exit(8);
  /* SCENE */
  Scene scene;
  /* num_meshes */
  if (fread(&scene.num_meshes, sizeof(uint32_t), 1, f) != 1) exit(8);
  if (scene.num_meshes == 0)
    PERROR_CLEANUP_EXIT("Scene has no meshes", 8);
  if (scene.num_meshes > 1000)
    PERROR_CLEANUP_EXIT("Scene has too many meshes", 8);
  /* meshes */
  scene.meshes = (Mesh*)malloc(scene.num_meshes * sizeof(Mesh));
  for (uint32_t i = 0; i != scene.num_meshes; ++i) {
    /* MESH */
    Mesh *mesh = &scene.meshes[i];
    /* num_vertices */
    if(fread(&mesh->num_vertices, sizeof(uint32_t), 1, f) != 1) exit(8);
    /* vertices */
    mesh->vs = (Vec3*)malloc(mesh->num_vertices * sizeof(Vec3));
    if (fread(mesh->vs, sizeof(float) * 3, mesh->num_vertices, f) != mesh->num_vertices)
      exit(8);
    /* num_triangles */
    if (fread(&mesh->num_triangles, sizeof(uint32_t), 1, f) != 1) exit(8);
    /* triangles */
    mesh->is = (uint32_t*)malloc(mesh->num_triangles * 3 * sizeof(uint32_t));
    if (fread(mesh->is, sizeof(uint32_t), mesh->num_triangles * 3, f) != mesh->num_triangles * 3)
      exit(8);
    /* material_index */
    if (fread(&mesh->material_index, sizeof(uint32_t), 1, f) != 1) exit(8);
    /* velocity */
    if (fread(&mesh->velocity, sizeof(Vec3), 1, f) != 1) exit(8);
    /* normals */
    mesh->ns = NULL;
  }

  fclose(f);
  return scene;
}

