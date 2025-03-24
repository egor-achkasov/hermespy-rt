/* vim: set tabstop=2:softtabstop=2:shiftwidth=2:noexpandtab */

#include "scene.h" /* for HRT_Scene, HRT_Mesh, HRT_Material */

#include <stdio.h> /* for FILE, fopen, fclose, fseek, fgets, sscanf, printf */
#include <stdlib.h> /* for malloc */
#include <stdint.h> /* for uint8_t, uint32_t */
#include <string.h> /* for strncmp, strrchr, strcmp */
#include <dirent.h> /* for opendir, readdir, closedir, DT_REG */


/**
 * Hardcoded scenes
 */

/* Box */
float mesh_box_vs[] = {
  5.f, 5.f, 0.f,
  -5.f, 5.f, 0.f,
  -5.f, -5.f, 0.f,
  5.f, -5.f, 0.f,
  5.f, 5.f, 5.f,
  -5.f, 5.f, 5.f,
  -5.f, -5.f, 5.f,
  5.f, -5.f, 5.f,
};
uint32_t mesh_box_is[] = {
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
HRT_Mesh mesh_box = {
  .num_vertices = 8,
  .vs = (float*)mesh_box_vs,
  .num_triangles = 12,
  .is = (uint32_t*)mesh_box_is,
  .material_index = MATERIAL_CONCRETE,
  .velocity = {0.f, 0.f, 0.f},
};
HRT_Scene scene_box = {
  .num_meshes = 1,
  .meshes = &mesh_box,
};
const char* scene_box_filename = "box.xml";

/* Simple reflector */
float mesh_simpleReflector_vs[] = {
  -0.5f, -0.5f, 0.f,
  0.5f, -0.5f, 0.f,
  0.5f, 0.5f, 0.f,
  -0.5f, 0.5f, 0.f,
};
uint32_t mesh_simpleReflector_is[] = {
  0, 1, 2,
  0, 2, 3,
};
HRT_Mesh mesh_simpleReflector = {
  .num_vertices = 4,
  .vs = (float*)mesh_simpleReflector_vs,
  .num_triangles = 2,
  .is = (uint32_t*)mesh_simpleReflector_is,
  .material_index = MATERIAL_CONCRETE,
  .velocity = {0.f, 0.f, 0.f},
};
HRT_Scene scene_simpleReflector = {
  .num_meshes = 1,
  .meshes = &mesh_simpleReflector,
};
const char* scene_simpleReflector_filename = "simple_reflector.xml";

/**
 * Read a Sionna scene
 */

/* Read a binary PLY file
 *
 * The PLY file have the following header:
 * ply
 * format binary_little_endian 1.0
 * element vertex <num_vertices>
 * property float x
 * property float y
 * property float z
 * property float s
 * property float t
 * element face <num_triangles>
 * property list uchar int vertex_index
 * end_header
 * 
 * \param filepath Path to the PLY file
 * \return The mesh
 */
HRT_Mesh readPly(const char* filepath) {
  FILE* f = fopen(filepath, "rb");
  if (f == NULL) {
    perror("Error: cannot open file");
    exit(8);
  }

  char buf[256];
  HRT_Mesh mesh = {
    .num_vertices = 0,
    .vs = NULL,
    .num_triangles = 0,
    .is = NULL,
    .material_index = 0,
    .velocity = {0.f, 0.f, 0.f},
  };

  /* read header */
  while(fgets(buf, 256, f)) {
    if (!strncmp(buf, "end_header", 10))
      break;
    if (!strncmp(buf, "element vertex ", 15)) {
      if (sscanf(buf, "element vertex %u", &mesh.num_vertices) != 1) {
        perror("Error: cannot read number of vertices");
        exit(8);
      }
    }
    else if (!strncmp(buf, "element face ", 13))
      if (sscanf(buf, "element face %u", &mesh.num_triangles) != 1) {
        perror("Error: cannot read number of triangles");
        exit(8);
      }
  }

  /* check if vertex and face elements are found */
  if (mesh.num_vertices == 0 || mesh.num_triangles == 0) {
    perror("Error: PLY element vertex or element face not found");
    exit(8);
  /* check if the numbers are not too big */
  } else if (mesh.num_vertices > 1000000 || mesh.num_triangles > 1000000) {
    perror("Error: PLY element vertex or element face too big");
    exit(8);
  }

  /* allocate memory for vertices and faces */
  mesh.vs = (float*)malloc(3 * mesh.num_vertices * sizeof(float));
  mesh.is = (uint32_t*)malloc(3 * mesh.num_triangles * sizeof(uint32_t));

  /* read vertices */
  for (uint32_t i = 0; i != mesh.num_vertices * 3; i += 3) {
    if (fread(&mesh.vs[i], sizeof(float), 3, f) != 3) {
      perror("Error: cannot read vertex");
      goto free_and_error;
    }
    /* skip s and t */
    if (fseek(f, 8, SEEK_CUR) != 0) {
      perror("Error: cannot skip s and t");
      goto free_and_error;
    }
  }

  /* read faces */
  for (uint32_t i = 0; i != mesh.num_triangles * 3; i += 3) {
    uint8_t n;
    if (fread(&n, sizeof(uint8_t), 1, f) != 1) {
      perror("Error: cannot read number of vertices in face");
      goto free_and_error;
    }
    if (n != 3) {
      perror("Error: face is not a triangle");
      goto free_and_error;
    }
    if (fread(&mesh.is[i], sizeof(uint32_t), 3, f) != 3) {
      perror("Error: cannot read face");
      goto free_and_error;
    }
  }

  fclose(f);
  return mesh;

free_and_error:
  free(mesh.vs);
  free(mesh.is);
  fclose(f);
  exit(8);
}

/** Read a scene config CSV file.
 * 
 * The CSV file must be a text file with the following format:
 * name,material_index,velocity_x,velocity_y,velocity_z
 * 
 * The file can have any number of lines except the first one, which is the header.
 * If a mesh is not in the CSV file,
 * it will have the default (0) material index and velocity.
 *
 * \param filepath Path to the CSV file
 * \param num_meshes Output number of meshes in the CSV file
 * \param mesh_filenames Output array of mesh filenames
 * \param material_indicies Output array of material indicies
 * \param velocity Output array of velocities. Size [num_meshes * 3]
 */
void readCsv(
  const char* filepath,
  uint32_t* num_meshes,
  char** mesh_filenames,
  uint32_t* material_indicies,
  float* velocity
)
{
  FILE* f = fopen(filepath, "r");
  if (f == NULL) {
    perror("Error: cannot open file");
    exit(8);
  }

  /* Check header */
  char buf[256];
  if (fgets(buf, 256, f) == NULL) {
    perror("Error: cannot read header");
    exit(8);
  }
  if (strncmp(buf, "name,material_index,velocity_x,velocity_y,velocity_z\n", 48)) {
    perror("Error: invalid header");
    exit(8);
  }

  /* Count number of lines */
  *num_meshes = 0;
  while (fgets(buf, 256, f))
    ++(*num_meshes);
  /* Return to the beginning of the file */
  rewind(f);
  /* Skip header */
  if (fgets(buf, 256, f) == NULL) {
    perror("Error: cannot read header");
    exit(8);
  }

  /* Allocate memory for the arrays */
  mesh_filenames = (char**)malloc(*num_meshes * sizeof(char*));
  material_indicies = (uint32_t*)malloc(*num_meshes * sizeof(uint32_t));
  velocity = (float*)malloc(*num_meshes * 3 * sizeof(float));

  /* Read lines */
  for (uint32_t i = 0; i != *num_meshes; ++i) {
    /* Read line */
    if (fgets(buf, 256, f) == NULL) {
      perror("Error: cannot read line");
      exit(8);
    }
    /* Parse line */
    char name[50];
    int rc = sscanf(
      buf,
      "%49[^,],%u,%f,%f,%f\n",
      name,
      &material_indicies[i],
      &velocity[i * 3],
      &velocity[i * 3 + 1],
      &velocity[i * 3 + 2]
    );
    if (rc != 4) {
      perror("Error: cannot parse line");
      exit(8);
    }
    /* Copy name */
    mesh_filenames[i] = (char*)malloc(strlen(name) + 1);
    if (mesh_filenames[i] == NULL) {
      perror("Error: cannot allocate memory for mesh filename");
      exit(8);
    }
    strcpy(mesh_filenames[i], name);
  }
}


/* Read a Sionna .xml scene. Assumes meshes are in a "meshes" directory.
 *
 * \param filepath Path to the .xml scene file
 * \param scene The scene to fill
 */
void readScene(const char* filepath, HRT_Scene* scene) {
  /* Check if filepath ends with ".xml" */
  size_t filepath_len = strlen(filepath);
  if (filepath_len > 4 && strcmp(filepath + filepath_len - 4, ".xml") != 0) {
    perror("Error: scene file must end with .xml");
    exit(8);
  }

  /* Read the scene config CSV file */
  char* csv_path = (char*)malloc(filepath_len + 1);
  if (!csv_path) {
    perror("Error: failed to allocate csv_path");
    exit(8);
  }
  strcpy(csv_path, filepath);
  strcpy(csv_path + filepath_len - 4, ".csv");
  uint32_t csv_num_meshes;
  char** csv_mesh_filenames;
  uint32_t* csv_material_indicies;
  float* csv_velocity;
  readCsv(
    csv_path,
    &csv_num_meshes,
    csv_mesh_filenames,
    csv_material_indicies,
    csv_velocity
  );

  /* Get the meshes directory path */
  const char* last_slash = strrchr(filepath, '/');
  size_t base_len = (last_slash ? (last_slash - filepath + 1) : 0);
  char* dir_path = (char*)malloc(base_len + strlen("meshes/") + 1);
  if (!dir_path) {
    perror("Error: failed to allocate dir_path");
    exit(8);
  }
  if (base_len > 0) {
    strncpy(dir_path, filepath, base_len);
    dir_path[base_len] = '\0';
  } else
    dir_path[0] = '\0'; /* No directory, start empty */
  strcat(dir_path, "meshes/");

  /* Scan the directory for .ply files */
  DIR* dir = opendir(dir_path);
  if (!dir) {
    perror("Failed to open meshes directory");
    free(dir_path);
    exit(8);
  }

  /* Count the number of .ply files */
  uint32_t ply_count = 0;
  struct dirent* entry;
  while ((entry = readdir(dir))) {
    if (entry->d_type != DT_REG) continue;
    size_t name_len = strlen(entry->d_name);
    if (name_len > 4 && strncmp(entry->d_name + name_len - 4, ".ply", 4) == 0)
      ply_count++;
  }

  /* Allocate memory for the meshes array */
  scene->num_meshes = ply_count;
  scene->meshes = (HRT_Mesh*)malloc(ply_count * sizeof(HRT_Mesh));
  if (!scene->meshes && ply_count > 0) {
    perror("Error: failed to allocate meshes array");
    free(dir_path);
    closedir(dir);
    exit(8);
  }

  /* Read each .ply file */
  rewinddir(dir); // Reset directory pointer
  uint32_t mesh_idx = 0;
  while ((entry = readdir(dir)) && mesh_idx < ply_count) {
    /* Continue if not a regular file or not a .ply file */
    size_t name_len = strlen(entry->d_name);
    if ((entry->d_type != DT_REG)
    || (name_len <= 4 || strcmp(entry->d_name + name_len - 4, ".ply")))
      continue;

    /* Construct full filepath for the .ply file */
    char* ply_path = (char*)malloc(strlen(dir_path) + name_len + 1);
    if (!ply_path) {
      perror("Error: failed to allocate ply_path");
      free(dir_path);
      closedir(dir);
      free(scene->meshes);
      exit(8);
    }
    strcpy(ply_path, dir_path);
    strcat(ply_path, entry->d_name);

    /* Read mesh */
    scene->meshes[mesh_idx] = readPly(ply_path);

    /* Get material index and velocity from the CSV file */
    for (uint32_t i = 0; i != csv_num_meshes; ++i)
      if (strcmp(entry->d_name, csv_mesh_filenames[i]) == 0) {
        scene->meshes[mesh_idx].material_index = csv_material_indicies[i];
        scene->meshes[mesh_idx].velocity[0] = csv_velocity[i * 3];
        scene->meshes[mesh_idx].velocity[1] = csv_velocity[i * 3 + 1];
        scene->meshes[mesh_idx].velocity[2] = csv_velocity[i * 3 + 2];
        break;
      }

    ++mesh_idx;
    free(ply_path);
  }

  // Cleanup
  closedir(dir);
  free(dir_path);
}


/**
 * Main
 */

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("Usage: %s <scene.xml>\n", argv[0]);
    return 1;
  }
  const char* scene_filepath = argv[1];
  const char* scene_filename = strrchr(scene_filepath, '/');
  if (scene_filename == NULL)
    exit(8);
  else
    ++scene_filename;

  HRT_Scene* scene_ptr;

  /* read scene */
  /* check if the scene is a hardcoded scene */
  if (!strcmp(scene_filename, scene_box_filename))
    scene_ptr = &scene_box;
  else if (!strcmp(scene_filename, scene_simpleReflector_filename))
    scene_ptr = &scene_simpleReflector;
  else {
    scene_ptr = (HRT_Scene*)malloc(sizeof(HRT_Scene));
    readScene(scene_filepath, scene_ptr);
  }

  FILE *fp = fopen("scene.hrt", "wb");
  if (fp == NULL) {
    perror("Error: cannot open file");
    return 1;
  }

  /* MAGIC */
  fwrite("HRT", 1, 3, fp);
  /* SCENE */
  /* num_meshes */
  fwrite(&scene_ptr->num_meshes, sizeof(uint32_t), 1, fp);
  /* meshes */
  for (uint32_t i = 0; i != scene_ptr->num_meshes; ++i) {
    /* MESH */
    fwrite(&scene_ptr->meshes[i].num_vertices, sizeof(uint32_t), 1, fp);
    fwrite(scene_ptr->meshes[i].vs, sizeof(float), 3 * scene_ptr->meshes[i].num_vertices, fp);
    fwrite(&scene_ptr->meshes[i].num_triangles, sizeof(uint32_t), 1, fp);
    fwrite(scene_ptr->meshes[i].is, sizeof(uint32_t), 3 * scene_ptr->meshes[i].num_triangles, fp);
    fwrite(&scene_ptr->meshes[i].material_index, sizeof(uint32_t), 1, fp);
    fwrite(&scene_ptr->meshes[i].velocity, sizeof(float), 3, fp);
  }

  fclose(fp);

  return 0;
}
