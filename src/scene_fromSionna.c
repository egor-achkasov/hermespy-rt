/* vim: set tabstop=2:softtabstop=2:shiftwidth=2:noexpandtab */

#include "../inc/vec3.h" /* for Vec3 */
#include "../inc/common.h" /* for IN, OUT, PERROR_CLEANUP_EXIT */
#include "../inc/scene.h" /* for Scene, Mesh, Material, scene_save */
#include "../inc/materials.h" /* for g_hrt_materials, MATERIAL_CONCRETE */

#include <stdio.h> /* for FILE, fopen, fclose, fseek, fgets, sscanf, printf */
#include <stdlib.h> /* for malloc */
#include <stdint.h> /* for uint8_t, uint32_t */
#include <string.h> /* for strncmp, strrchr, strcmp, snprintf, strcpy, strlen, strncpy, strchr, strstr */

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
Mesh mesh_box = {
  .num_vertices = 8,
  .vs = (Vec3*)mesh_box_vs,
  .num_triangles = 12,
  .is = (uint32_t*)mesh_box_is,
  .material_index = MATERIAL_CONCRETE,
  .velocity = {0.f, 0.f, 0.f},
};
Scene scene_box = {
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
Mesh mesh_simpleReflector = {
  .num_vertices = 4,
  .vs = (Vec3*)mesh_simpleReflector_vs,
  .num_triangles = 2,
  .is = (uint32_t*)mesh_simpleReflector_is,
  .material_index = MATERIAL_CONCRETE,
  .velocity = {0.f, 0.f, 0.f},
};
Scene scene_simpleReflector = {
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
Mesh readPly(const char* filepath) {
  FILE* f = fopen(filepath, "rb");
  if (f == NULL)
    PERROR_CLEANUP_EXIT("Error: cannot open file", 8);

  char buf[256];
  Mesh mesh = {
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
      if (sscanf(buf, "element vertex %u", &mesh.num_vertices) != 1)
        PERROR_CLEANUP_EXIT("Error: cannot read number of vertices", 8);
    }
    else if (!strncmp(buf, "element face ", 13))
      if (sscanf(buf, "element face %u", &mesh.num_triangles) != 1)
        PERROR_CLEANUP_EXIT("Error: cannot read number of triangles", 8);
  }

  /* check if vertex and face elements are found */
  if (mesh.num_vertices == 0 || mesh.num_triangles == 0)
    PERROR_CLEANUP_EXIT("Error: PLY element vertex or element face not found", 8);
  /* check if the numbers are not too big */
  else if (mesh.num_vertices > 1000000 || mesh.num_triangles > 1000000)
    PERROR_CLEANUP_EXIT("Error: PLY element vertex or element face too big", 8);

  /* allocate memory for vertices and faces */
  mesh.vs = (Vec3*)malloc(mesh.num_vertices * sizeof(Vec3));
  mesh.is = (uint32_t*)malloc(3 * mesh.num_triangles * sizeof(uint32_t));

  /* read vertices */
  for (uint32_t i = 0; i != mesh.num_vertices; ++i) {
    if (fread(&mesh.vs[i], sizeof(Vec3), 1, f) != 1)
      PERROR_CLEANUP_EXIT("Error: cannot read vertex", 8, mesh.vs, mesh.is);
    /* skip s and t */
    if (fseek(f, 8, SEEK_CUR) != 0)
      PERROR_CLEANUP_EXIT("Error: cannot skip s and t", 8, mesh.vs, mesh.is);
  }

  /* read faces */
  for (uint32_t i = 0; i != mesh.num_triangles * 3; i += 3) {
    uint8_t n;
    if (fread(&n, sizeof(uint8_t), 1, f) != 1)
      PERROR_CLEANUP_EXIT("Error: cannot read number of vertices in face", 8, mesh.vs, mesh.is);
    if (n != 3)
      PERROR_CLEANUP_EXIT("Error: face is not a triangle", 8, mesh.vs, mesh.is);
    if (fread(&mesh.is[i], sizeof(uint32_t), 3, f) != 3)
      PERROR_CLEANUP_EXIT("Error: cannot read face", 8, mesh.vs, mesh.is);
  }

  fclose(f);
  return mesh;
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
 * \param velocity Output array of velocity vectors. Size [num_meshes]
 */
void readCsv(
  IN const char* filepath,
  OUT uint32_t* csv_num_meshes,
  OUT char*** csv_mesh_filenames,
  OUT uint32_t** csv_material_indicies,
  OUT Vec3** csv_velocities
)
{
  FILE* f = fopen(filepath, "r");
  if (!f)
    PERROR_CLEANUP_EXIT("Error: cannot open file", 8);

  /* Check header */
  char buf[256];
  if (!fgets(buf, 256, f))
    PERROR_CLEANUP_EXIT("Error: cannot read header", 8);
  if (strncmp(buf, "name,material_index,velocity_x,velocity_y,velocity_z\n", 48))
    PERROR_CLEANUP_EXIT("Error: invalid header", 8);

  /* Count number of lines */
  *csv_num_meshes = 0;
  while (fgets(buf, 256, f))
    ++(*csv_num_meshes);
  if (*csv_num_meshes == 0) {
    fclose(f);
    return;
  }
  /* Return to the beginning of the file */
  rewind(f);
  /* Skip header */
  if (!fgets(buf, 256, f))
    PERROR_CLEANUP_EXIT("Error: cannot read header", 8);

  /* Allocate memory for the arrays */
  *csv_mesh_filenames = (char**)malloc(*csv_num_meshes * sizeof(char*));
  *csv_material_indicies = (uint32_t*)malloc(*csv_num_meshes * sizeof(uint32_t));
  *csv_velocities = (Vec3*)malloc(*csv_num_meshes * sizeof(Vec3));

  /* Read lines */
  for (uint32_t i = 0; i != *csv_num_meshes; ++i) {
    /* Read line */
    if (!fgets(buf, 256, f))
      PERROR_CLEANUP_EXIT("Error: cannot read line", 8);
    /* Parse line */
    char name[50];
    int rc = sscanf(
      buf,
      "%49[^,],%u,%f,%f,%f\n",
      name,
      &(*csv_material_indicies)[i],
      &(*csv_velocities)[i].x,
      &(*csv_velocities)[i].y,
      &(*csv_velocities)[i].z
    );
    if (rc != 4)
      PERROR_CLEANUP_EXIT("Error: cannot parse line", 8);
    /* Copy name */
    (*csv_mesh_filenames)[i] = (char*)malloc(strlen(name) + 1);
    if (!csv_mesh_filenames[i])
      PERROR_CLEANUP_EXIT("Error: cannot allocate memory for mesh filename", 8);
    strcpy((*csv_mesh_filenames)[i], name);
  }
}

/** Read a Sionna scene XML file.
 * 
 * Extracts each shape's name, file path and material.
 * 
 * \param filepath Path to the XML file
 * \param xml_num_meshes Output number of meshes in the XML file
 * \param xml_mesh_filepaths Output array of mesh file paths (null-terminated strings)
 * \param xml_mesh_names Output array of mesh names (null-terminated strings)
 * \param xml_mesh_materials Output array of mesh materials (null-terminated strings)
 */
void readXml(
  IN const char* filepath,
  OUT uint32_t* xml_num_meshes,
  OUT char*** xml_mesh_filepaths,
  OUT char*** xml_mesh_names,
  OUT char*** xml_mesh_materials
)
{
  FILE* f = fopen(filepath, "r");
  if (f == NULL) 
    PERROR_CLEANUP_EXIT("Error: cannot open the xml file", 8);

  fseek(f, 0, SEEK_END);
  size_t fsize = (size_t)ftell(f);
  fseek(f, 0, SEEK_SET);

  char* buf = (char*)malloc(fsize + 1);
  if (!buf)
    PERROR_CLEANUP_EXIT("Error: cannot allocate memory for the xml file", 8);  
  if (fread(buf, 1, fsize, f) != fsize)
    PERROR_CLEANUP_EXIT("Error: cannot read the xml file", 8, buf);
  buf[fsize] = '\0';
  fclose(f);

  /* Count the amount of "<shape" entries */
  char** shape_positions = NULL;
  *xml_num_meshes = 0;
  char* pos = buf, *end;
  while ((pos = strstr(pos, "<shape")) != NULL) {
    ++*xml_num_meshes;
    shape_positions = (char**)realloc(shape_positions, *xml_num_meshes * sizeof(char*));
    if (!shape_positions)
      PERROR_CLEANUP_EXIT("Error: cannot reallocate memory for shape positions", 8, buf);
    shape_positions[*xml_num_meshes - 1] = pos;
    pos += 6;
  }
  if (*xml_num_meshes == 0)
    PERROR_CLEANUP_EXIT("Error: no shapes found in the xml file", 8, buf);
  if (!shape_positions)
    PERROR_CLEANUP_EXIT("Error: cannot allocate memory for shape positions", 8, buf);

  /* Allocate memory for the mesh data */
  *xml_mesh_filepaths = (char**)malloc(*xml_num_meshes * sizeof(char*));
  *xml_mesh_names = (char**)malloc(*xml_num_meshes * sizeof(char*));
  *xml_mesh_materials = (char**)malloc(*xml_num_meshes * sizeof(char*));
  if (!*xml_mesh_filepaths || !*xml_mesh_names || !*xml_mesh_materials)
    PERROR_CLEANUP_EXIT("Error: cannot allocate memory for mesh data", 8, buf);

  /* For each "<shape" entry, extract the name, file path and material */
  char *e_name, *e_filepath, *e_material;
  for (uint32_t i = 0; i != *xml_num_meshes; ++i) {
    pos = shape_positions[i] + 6;

    /* Find mesh name */
    if ((pos = strstr(pos, "name=\"")) == NULL)
      PERROR_CLEANUP_EXIT("Error: cannot find mesh name", 8,
        buf, shape_positions, *xml_mesh_filepaths, *xml_mesh_names, *xml_mesh_materials);
    pos += 6;
    end = strchr(pos, '\"');
    if (!end)
      PERROR_CLEANUP_EXIT("Error: cannot find mesh name end", 8,
        buf, shape_positions, *xml_mesh_filepaths, *xml_mesh_names, *xml_mesh_materials);
    e_name = (char*)malloc(end - pos + 1);
    if (!e_name)
      PERROR_CLEANUP_EXIT("Error: cannot allocate memory for mesh name", 8,
        buf, shape_positions, *xml_mesh_filepaths, *xml_mesh_names, *xml_mesh_materials);
    strncpy(e_name, pos, end - pos);
    e_name[end - pos] = '\0';

    /* Find mesh file path */
    if (!(pos = strstr(pos, "<string name=\"filename\""))) 
      PERROR_CLEANUP_EXIT("Error: cannot find mesh file path", 8, buf, e_name,
        shape_positions, *xml_mesh_filepaths, *xml_mesh_names, *xml_mesh_materials);
    if (!(pos = strstr(pos, "value=\"")))
      PERROR_CLEANUP_EXIT("Error: cannot find mesh file path value", 8, buf, e_name,
        shape_positions, *xml_mesh_filepaths, *xml_mesh_names, *xml_mesh_materials);
    pos += 7;
    end = strchr(pos, '\"');
    if (!end)
      PERROR_CLEANUP_EXIT("Error: cannot find mesh file path end", 8, buf, e_name,
        shape_positions, *xml_mesh_filepaths, *xml_mesh_names, *xml_mesh_materials);
    e_filepath = (char*)malloc(end - pos + 1);
    if (!e_filepath)
      PERROR_CLEANUP_EXIT("Error: cannot allocate memory for mesh file path", 8, buf, e_name,
        shape_positions, *xml_mesh_filepaths, *xml_mesh_names, *xml_mesh_materials);
    strncpy(e_filepath, pos, end - pos);
    e_filepath[end - pos] = '\0';

    /* Find mesh material */
    if (!(pos = strstr(pos, "id=\"mat-itu_")))
      PERROR_CLEANUP_EXIT("Error: cannot find mesh material", 8, buf, e_name, e_filepath,
        shape_positions, *xml_mesh_filepaths, *xml_mesh_names, *xml_mesh_materials);
    pos += 12;
    end = strchr(pos, '\"');
    if (!end)
      PERROR_CLEANUP_EXIT("Error: cannot find mesh material end", 8, buf, e_name, e_filepath,
        shape_positions, *xml_mesh_filepaths, *xml_mesh_names, *xml_mesh_materials);
    e_material = (char*)malloc(end - pos + 1);
    if (!e_material)
      PERROR_CLEANUP_EXIT("Error: cannot allocate memory for mesh material", 8, buf, e_name, e_filepath,
        shape_positions, *xml_mesh_filepaths, *xml_mesh_names, *xml_mesh_materials);
    strncpy(e_material, pos, end - pos);
    e_material[end - pos] = '\0';

    /* Save mesh data */
    (*xml_mesh_filepaths)[i] = e_filepath;
    (*xml_mesh_names)[i] = e_name;
    (*xml_mesh_materials)[i] = e_material;
  }

  free(shape_positions);
  free(buf);
}

/* Read a Sionna .xml scene. Assumes meshes are in a "meshes" directory.
 *
 * \param filepath Path to the .xml scene file
 * \param scene The scene to fill
 */
void readScene(const char* filepath, Scene* scene) {
  /* Check if filepath ends with ".xml" */
  size_t filepath_len = strlen(filepath);
  if (filepath_len > 4 && strcmp(filepath + filepath_len - 4, ".xml") != 0)
    PERROR_CLEANUP_EXIT("Error: scene file must end with .xml", 8);

  /* Read the xml file */
  uint32_t xml_num_meshes;
  char** xml_mesh_filepaths;
  char** xml_mesh_names;
  char** xml_mesh_materials;
  readXml(
    filepath,
    &xml_num_meshes,
    &xml_mesh_filepaths,
    &xml_mesh_names,
    &xml_mesh_materials
  );

  /* Read the scene config CSV file */
  char* csv_path = (char*)malloc(filepath_len + 1);
  if (!csv_path)
    PERROR_CLEANUP_EXIT("Error: cannot allocate memory for csv_path", 8);
  strcpy(csv_path, filepath);
  strcpy(csv_path + filepath_len - 4, ".csv");
  uint32_t csv_num_meshes;
  char** csv_mesh_filenames;
  uint32_t* csv_material_indicies;
  Vec3* csv_velocities;
  readCsv(
    csv_path,
    &csv_num_meshes,
    &csv_mesh_filenames,
    &csv_material_indicies,
    &csv_velocities
  );
  free(csv_path);

  /* Read the meshes */

  /* Get the scene directory path */
  size_t scene_dir_len = strrchr(filepath, '/') - filepath + 1;
  /* if filepath is "/path/to/scene.xml", scene_dir is "/path/to/" */
  char* scene_dir = (char*)malloc(scene_dir_len + 7);
  if (!scene_dir)
    PERROR_CLEANUP_EXIT("Error: failed to allocate scene_dir", 8);
  strncpy(scene_dir, filepath, scene_dir_len);

  /* For each ply file from the xml file */
  scene->num_meshes = xml_num_meshes;
  scene->meshes = (Mesh*)malloc(xml_num_meshes * sizeof(Mesh));
  if (!scene->meshes && xml_num_meshes > 0)
    PERROR_CLEANUP_EXIT("Error: failed to allocate meshes array", 8);
  for (uint32_t i = 0; i != xml_num_meshes; ++i) {
    /* Get the mesh file path */
    size_t mesh_filepath_len = scene_dir_len + strlen(xml_mesh_filepaths[i]) + 1;
    char* mesh_filepath = (char*)malloc(mesh_filepath_len);
    snprintf(
      mesh_filepath,
      mesh_filepath_len,
      "%.*s%s",
      (int)scene_dir_len,
      scene_dir,
      xml_mesh_filepaths[i]
    );
    /* Read the mesh */
    scene->meshes[i] = readPly(mesh_filepath);
    free(mesh_filepath);
    /* Set the material index from the material name from xml */
    scene->meshes[i].material_index = get_material_index(xml_mesh_materials[i]);
    /* Set the material index and velocity from CSV if present */
    for (uint32_t j = 0; j != csv_num_meshes; ++j)
      if (strcmp(xml_mesh_names[i], csv_mesh_filenames[j]) == 0) {
        scene->meshes[i].material_index = csv_material_indicies[j];
        scene->meshes[i].velocity = csv_velocities[j];
        break;
      }
  }

  free(scene_dir);
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
  if (scene_filename == NULL) exit(8);
  else ++scene_filename;

  Scene* scene_ptr;

  /* read scene */
  /* check if the scene is a hardcoded scene */
  if (!strcmp(scene_filename, scene_box_filename))
    scene_ptr = &scene_box;
  else if (!strcmp(scene_filename, scene_simpleReflector_filename))
    scene_ptr = &scene_simpleReflector;
  else {
    scene_ptr = (Scene*)malloc(sizeof(Scene));
    readScene(scene_filepath, scene_ptr);
  }

  /* Save the scene as HRT */
  scene_save(scene_ptr, "scene.hrt");

  return 0;
}
