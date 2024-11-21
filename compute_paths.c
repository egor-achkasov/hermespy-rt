#include "compute_paths.h"

#include <libxml/parser.h>
#include <libxml/tree.h>

#include <stddef.h> /* for size_t */
#include <stdlib.h> /* for exit, malloc, free */
#include <libgen.h> /* for dirname */
#include <string.h> /* for strlen, sprintf */
#include <stdio.h>  /* for fopen, FILE, fclose */
#include <math.h>   /* for sin, cos, sqrt */

#define PI 3.14159265358979323846
#define EPS 1e-6

typedef struct {
    float x, y, z;
} Vec3;

typedef struct {
    Vec3 o, d;
} Ray;

typedef struct {
    char *name;
    Vec3 reflectance;
} Material;

typedef struct {
    Vec3 *vertices;
    size_t num_vertices;
    int32_t *indices;
    size_t num_indices;
    Vec3 *normals;
    Material *material;
} Mesh;

typedef struct {
    Mesh *meshes[256]; /* TODO: dynamic allocation, fix possible overflow */
    size_t num_meshes;
    Material *materials[256]; /* TODO: dynamic allocation, fix possible overflow */
    size_t num_materials;
} Scene;

/* ==== VECTOR OPERATIONS ==== */

Vec3 vec3_sub(const Vec3 *a, const Vec3 *b)
{
    return (Vec3){a->x - b->x, a->y - b->y, a->z - b->z};
}
Vec3 vec3_add(const Vec3 *a, const Vec3 *b)
{
    return (Vec3){a->x + b->x, a->y + b->y, a->z + b->z};
}
Vec3 vec3_cross(const Vec3 *a, const Vec3 *b)
{
    return (Vec3){
        a->y * b->z - a->z * b->y,
        a->z * b->x - a->x * b->z,
        a->x * b->y - a->y * b->x
    };
}
float vec3_dot(const Vec3 *a, const Vec3 *b)
{
    return a->x * b->x + a->y * b->y + a->z * b->z;
}
Vec3 vec3_scale(const Vec3 *a, float s)
{
    return (Vec3){a->x * s, a->y * s, a->z * s};
}

/* ==== MESH LOADING ==== */

/** Parse a Mitsuba XML "bsdf" node with a diffuse material.
 * Extracts "rgb" property with "reflectance" name.
 * 
 * \param node the xml node
 * \param scene the scene to populate. Parsed material will be added to scene->materials.
 */
void parse_xmlnode_bsdf(IN xmlNodePtr node, OUT Scene *scene)
{
    char *material_name;
    xmlChar *prop;
    xmlNodePtr child, grandchild;

    /* <bsdf type="diffuse" name="bsdf"> */
    child = node->children;
    while (child && child->type != XML_ELEMENT_NODE)
        child = child->next;
    if (!child) return;
    if (xmlStrcmp(child->name, (const xmlChar *)"bsdf")
    || xmlStrcmp(xmlGetProp(child, (const xmlChar *)"type"), (const xmlChar *)"diffuse"))
        return;

    /* <rgb value="0.603827 0.090842 0.049707" name="reflectance"/> */
    grandchild = child->children;
    while (grandchild && grandchild->type != XML_ELEMENT_NODE)
        grandchild = grandchild->next;
    if (!grandchild) return;
    if (xmlStrcmp(grandchild->name, (const xmlChar *)"rgb")
    || xmlStrcmp(xmlGetProp(grandchild, (const xmlChar *)"name"), (const xmlChar *)"reflectance"))
        return;
    prop = xmlGetProp(grandchild, (const xmlChar *)"value");
    if (!prop) return;

    /* Get the material name from the original node */
    material_name = (char*)xmlGetProp(node, (const xmlChar *)"id");
    if (!material_name) return;

    /* Create the material object */
    scene->materials[scene->num_materials] = (Material*)malloc(sizeof(Material));
    scene->materials[scene->num_materials]->name = material_name;
    sscanf((char*)prop, "%f %f %f",
        &scene->materials[scene->num_materials]->reflectance.x,
        &scene->materials[scene->num_materials]->reflectance.y,
        &scene->materials[scene->num_materials]->reflectance.z);    
    xmlFree(prop);
    ++scene->num_materials;
}

/** Load a mesh from a PLY file.
 * 
 * Supports only format binary_little_endian 1.0.
 * 
 * All faces must be triangles, meaning
 * indices set as "property list uchar int"
 * with each uchar being 3.
 * 
 * Vertices must contain float x, y, z, u and v coordinates.
 * 
 * \param mesh_filepath path to the PLY file
 * \return the loaded mesh
 */
Mesh* load_mesh_ply(const char *mesh_filepath)
{
    FILE *f = fopen(mesh_filepath, "rb");
    if (!f) {
        fprintf(stderr, "Could not open file %s\n", mesh_filepath);
        exit(8);
    }

    char buff[256];

    /* HEADER */
    /* ply */
    fread(buff, 1, 4, f);
    if (strncmp(buff, "ply\n", 4)) exit(8);
    /* format binary_little_endian 1.0 */
    fread(buff, 1, 32, f);
    if (strncmp(buff, "format binary_little_endian 1.0\n", 32)) exit(8);

    Mesh *mesh = (Mesh*)malloc(sizeof(Mesh));
    while (fgets(buff, sizeof(buff), f)) {
        if (strncmp(buff, "element vertex", 14) == 0)
            sscanf(buff, "element vertex %zu", &mesh->num_vertices);
        else if (strncmp(buff, "element face", 12) == 0)
            sscanf(buff, "element face %zu", &mesh->num_indices);
        else if (strncmp(buff, "end_header", 10) == 0)
            break;
    }
    mesh->num_indices *= 3;
    mesh->vertices = (Vec3*)malloc(mesh->num_vertices * sizeof(Vec3));
    mesh->indices = (int32_t*)malloc(mesh->num_indices * sizeof(int32_t));
    mesh->normals = (Vec3*)malloc(mesh->num_indices / 3 * sizeof(Vec3));

    /* VERTICES */
    for (size_t i = 0; i < mesh->num_vertices; ++i) {
        fread(&mesh->vertices[i], sizeof(Vec3), 1, f);
        fseek(f, 8, SEEK_CUR); /* skip u and v */
    }

    /* INDICES */
    /* also calculate normals */
    uint8_t n;
    Vec3 *v1, *v2, *v3, u, v;
    for (size_t i = 0; i < mesh->num_indices; i += 3) {
        fread(&n, 1, 1, f);
        if (n != 3) {
            fprintf(stderr, "Only trianglar faces are supported\n");
            exit(8);
        }
        fread(&mesh->indices[i], sizeof(int32_t), 3, f);
        /* calculate normal */
        v1 = &mesh->vertices[mesh->indices[i]];
        v2 = &mesh->vertices[mesh->indices[i + 1]];
        v3 = &mesh->vertices[mesh->indices[i + 2]];
        u = vec3_sub(v2, v1);
        v = vec3_sub(v3, v1);
        mesh->normals[i / 3] = vec3_cross(&u, &v);
    }

    fclose(f);
    return mesh;
}

/** Parse a Mitsuba XML "shape" node with a mesh.
 * Assumes all the materials are already loaded in the scene.
 * 
 * \param node the xml node
 * \param dir the directory of the scene file
 * \param scene_filepath the path to the scene file
 * \param scene the scene to populate. Parsed mesh will be added to scene->meshes.
 * \param material_name the name of the material
 */
void parse_xmlnode_shape(
    IN xmlNodePtr node,
    IN char* dir,
    IN const char *scene_filepath,
    OUT Scene *scene)
{
    xmlChar *prop;
    xmlNodePtr child;
    char *mesh_filepath = NULL;
    char *material_name = NULL;

    prop = xmlGetProp(node, (const xmlChar *)"type");
    if (!prop || xmlStrcmp(prop, (const xmlChar *)"ply")) {
        xmlFree(prop);
        return NULL;
    }
    xmlFree(prop);
    /* <shape type="ply" id="... found, now get the mesh */
    for (child = node->children;
        child && !(material_name && mesh_filepath);
        child = child->next)
    {
        if (child->type != XML_ELEMENT_NODE)
            continue;
        /* get material name */
        if (!xmlStrcmp(child->name, (const xmlChar *)"ref")) {
            material_name = (char*)xmlGetProp(child, (const xmlChar *)"id");
            continue;
        }
        /* get mesh filepath */
        if (xmlStrcmp(child->name, (const xmlChar *)"string"))
            continue;
        prop = xmlGetProp(child, (const xmlChar *)"name");
        if (!prop || xmlStrcmp(prop, (const xmlChar *)"filename"))
            continue;
        xmlFree(prop);
        prop = xmlGetProp(child, (const xmlChar *)"value");
        if (!prop) continue;
        mesh_filepath = (char*)malloc(strlen(dir) + xmlStrlen(prop) + 2);
        sprintf(mesh_filepath, "%s/%s", dir, prop);
    }
    if (!mesh_filepath || !material_name) return;

    /* load the mesh */
    scene->meshes[scene->num_meshes] = load_mesh_ply(mesh_filepath);
    free(mesh_filepath);
    xmlFree(prop);

    /* Find the material */
    for (size_t i = 0; i < scene->num_materials; ++i)
        if (!strcmp(scene->materials[i]->name, material_name))
            scene->meshes[scene->num_meshes]->material = scene->materials[i];
    if (!scene->meshes[scene->num_meshes]->material) {
        fprintf(stderr, "Material %s not found\n", material_name);
        exit(8);
    }
    ++scene->num_meshes;
}

/** Load a scene from a Mitsuba XML file.
 * 
 * \param scene_filepath path to the Mitsuba XML file
 * \return the loaded mesh
 */
Scene load_scene(const char *scene_filepath)
{
    Scene scene = {0};
    xmlNodePtr root;
    xmlNodePtr cur = NULL;
    scene.num_meshes = 0;

    /* Open the xml file */
    xmlDocPtr doc = xmlReadFile(scene_filepath, NULL, XML_PARSE_NOBLANKS);
    if (!doc) {
        fprintf(stderr, "Could not parse file %s\n", scene_filepath);
        exit(8);
    }
    root = xmlDocGetRootElement(doc);
    if (!root) {
        fprintf(stderr, "Could not get root element of %s\n", scene_filepath);
        exit(8);
    }

    /* Get materials */
    for (cur = root->children; cur; cur = cur->next) {
        if (cur->type != XML_ELEMENT_NODE)
            continue;
        if (!xmlStrcmp(cur->name, (const xmlChar *)"bsdf"))
            parse_xmlnode_bsdf(cur, &scene);
    }

    /* Get shapes */
    char* dir = dirname((char*)scene_filepath);
    for (cur = root->children; cur; cur = cur->next) {
        if (cur->type != XML_ELEMENT_NODE)
            continue;
        if (!xmlStrcmp(cur->name, (const xmlChar *)"shape"))
            parse_xmlnode_shape(cur, dir, scene_filepath, &scene);
    }

    return scene;
}

/* ==== RT ==== */

/** Compute Moeleer-Trumbore intersection algorithm.
 * 
 * \param ray the ray
 * \param mesh the mesh
 * \param t output distance to the hit point. If no hit, t is not modified.
 * \param ind output index of the hit triangle. If no hit, i is not modified.
 */
void moeller_trumbore(Ray *ray, Mesh *mesh, float *t, size_t *ind)
{
    /* TODO BVH */
    /* for each triangle */
    Vec3 v1, v2, v3, e1, e2, re2_cross, s, se1_cross;
    float d, u, v;
    float dist = 1e9;
    float dist_tmp;
    for (size_t i = 0; i < mesh->num_indices; i += 3) {
        v1 = mesh->vertices[mesh->indices[i]];
        v2 = mesh->vertices[mesh->indices[i + 1]];
        v3 = mesh->vertices[mesh->indices[i + 2]];
        e1 = vec3_sub(&v2, &v1);
        e2 = vec3_sub(&v3, &v1);
        re2_cross = vec3_cross(&ray->d, &e2);
        d = vec3_dot(&e1, &re2_cross);
        if (d > -EPS && d < EPS) continue;
        s = vec3_sub(&ray->o, &v1);
        u = vec3_dot(&s, &re2_cross) / d;
        if ((u < 0. && fabs(u) > EPS)
        || (u > 1. && fabs(u - 1.) > EPS))
            continue;
        se1_cross = vec3_cross(&s, &e1);
        v = vec3_dot(&ray->d, &se1_cross) / d;
        if ((v < 0. && fabs(v) > EPS)
        || (u + v > 1. && fabs(u + v - 1.) > EPS))
            continue;
        dist_tmp = vec3_dot(&e2, &se1_cross) / d;
        if (dist_tmp > EPS && dist_tmp < dist) {
            dist = dist_tmp;
            *t = dist;
            *ind = i;
        }
    }
}

/* ==== MAIN FUNCTION ==== */

void compute_paths(
    IN const char *scene_filepath, /* path to the scene file */
    IN const float *rx_positions,  /* shape (num_rx, 3) */
    IN const float *tx_positions,  /* shape (num_tx, 3) */
    IN const float *rx_velocities, /* shape (num_rx, 3) */
    IN const float *tx_velocities, /* shape (num_tx, 3) */
    IN float carrier_frequency,    /* > 0.0 */
    IN int num_rx,                 /* number of receivers */
    IN int num_tx,                 /* number of transmitters */
    IN int num_paths,              /* number of paths */
    IN int num_bounces,            /* number of bounces */
    OUT float *a,                  /* output array of gains (num_paths,) */
    OUT int32_t *tau               /* output array of delays (num_paths,) */
)
{
    /* Load the scene */
    Scene scene = load_scene(scene_filepath);

    /* Calculate fibonacci sphere */
    Vec3 *ray_directions = (Vec3*)malloc(num_paths * sizeof(Vec3));
    float k, phi, theta;
    for (size_t i = 0; i < num_paths; ++i) {
        k = (float)i + .5;
        phi = 1. - 2. * k / num_paths;
        theta = PI * (1. + sqrt(5.)) * k;
        ray_directions[i] = (Vec3){
            cos(theta) * sin(phi),
            sin(theta) * sin(phi),
            cos(phi)
        };
    }

    /* Create num_path rays for each tx */
    /* Shape (num_tx, num_paths) */
    Ray **rays = (Ray**)malloc(num_tx * sizeof(Ray*));
    Vec3 tx_position;
    for (size_t i = 0; i < num_tx; ++i) {
        rays[i] = (Ray*)malloc(num_paths * sizeof(Ray));
        tx_position = (Vec3){tx_positions[i * 3], tx_positions[i * 3 + 1], tx_positions[i * 3 + 2]};
        for (size_t j = 0; j < num_paths; ++j) {
            rays[i][j].o = tx_position;
            rays[i][j].d = ray_directions[j];
        }
    }
    free(ray_directions);

    /* Calculate the paths */
    /* TODO calculate a and tau in moeller_trumbore */
    /* shape (num_bounces, num_tx, num_paths) */
    Vec3 ***hits = (Vec3***)malloc(num_bounces * sizeof(Vec3**));
    size_t ***hit_indices = (size_t***)malloc(num_bounces * sizeof(size_t**));
    float t;
    size_t ind;
    for (size_t i = 0; i < num_bounces; ++i) {
        hits[i] = (Vec3**)malloc(num_tx * sizeof(Vec3*));
        hit_indices[i] = (size_t**)malloc(num_tx * sizeof(size_t*));
        for (size_t j = 0; j < num_tx; ++j) {
            hits[i][j] = (Vec3*)malloc(num_paths * sizeof(Vec3));
            hit_indices[i][j] = (size_t*)malloc(num_paths * sizeof(size_t));
            for (size_t k = 0; k < num_paths; ++k)
                for (size_t l = 0; l < scene.num_meshes; ++l) {
                    t = -1.;
                    ind = -1;
                    moeller_trumbore(&rays[j][k], scene.meshes[l], &t, &ind);
                    hits[i][j][k] = vec3_scale(&rays[j][k].d, t);
                    hits[i][j][k] = vec3_add(&hits[i][j][k], &rays[j][k].o);
                    rays[j][k].o = hits[i][j][k];
                    hit_indices[i][j][k] = ind;
                }
        }
    }

    /* TODO Remove */
    a = hits;
    tau = hit_indices;

    /* Free */
    /* Free the hits */
    /*
    for (size_t i = 0; i < num_bounces; ++i) {
        for (size_t j = 0; j < num_tx; ++j) {
            free(hits[i][j]);
            free(hit_indices[i][j]);
        }
        free(hits[i]);
        free(hit_indices[i]);
    }
    free(hits);
    free(hit_indices);
    */
    /* Free the rays */
    for (size_t i = 0; i < num_tx; ++i)
        free(rays[i]);
    free(rays);
    /* Free the scene */
    for (size_t i = 0; i < scene.num_meshes; ++i) {
        free(scene.meshes[i]->vertices);
        free(scene.meshes[i]->indices);
        free(scene.meshes[i]->normals);
        free(scene.meshes[i]);
    }
}
