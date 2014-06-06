/**
 * \file    src/plex.c
 * \author  Scott Wales <scott.wales@unimelb.edu.au>
 *
 * Copyright 2014 ARC Centre of Excellence for Climate Systems Science
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <mpi.h>
#include <petsc.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef struct {
    size_t   vertex_count;
    double * vertex_coords; // Co-ordinates of each vector (length vertex_count*dimension)
    size_t   face_count;
    int    * face_vertices; // Surrounding vertices of each face (length face_count*corners)

    size_t   dimension;     // Dimension of each vertex
    size_t   corners;       // Number of corners in a cell
} meshdata;

// Write to stderr with MPI rank & function name
#define log(fmt, ...) \
    do { \
        int _msize, _mrank; \
        MPI_Comm_rank(MPI_COMM_WORLD,&_mrank); \
        MPI_Comm_size(MPI_COMM_WORLD,&_msize); \
        fprintf(stderr, "[%d/%d] %s:\t" fmt, _mrank, _msize, __FUNCTION__, __VA_ARGS__); \
    } while (0)

// Create a single square
meshdata BasicMesh() {
    meshdata data;

    data.vertex_count  = 4;
    data.dimension     = 2;
    data.vertex_coords = calloc(data.dimension*data.vertex_count,
                                sizeof(*(data.vertex_coords)));

    double coords[] = {
        0,0,  0,1,
        1,0,  1,1,
    };
    memcpy(data.vertex_coords, coords, 
           data.vertex_count*data.dimension*sizeof(*(data.vertex_coords)));

    data.face_count    = 1;
    data.corners       = 4;
    data.face_vertices = calloc(data.corners*data.vertex_count,
                                sizeof(*(data.face_vertices)));

    for (size_t i=0; i<4;++i) {
        data.face_vertices[i] = i;
    }

    return data;
}

// Read Mike's data
meshdata FileMesh() {
    FILE * meshfile = fopen("CubedSphere.mesh","r");

    meshdata data;
    fscanf(meshfile," %zu ",&(data.vertex_count));
    assert(data.vertex_count == 1572866);
    data.dimension     = 3;
    data.vertex_coords = calloc(data.dimension*data.vertex_count,
                                sizeof(*(data.vertex_coords)));

    for (size_t i=0; i<data.dimension*data.vertex_count;++i) {
        int count = fscanf(meshfile, " %lf ", &(data.vertex_coords[i]));
        if (count != 1) {
            perror("Reading vertices");
        }
    }

    fscanf(meshfile," %zu ",&(data.face_count));
    assert(data.face_count == 1572864);
    data.corners       = 4;
    data.face_vertices = calloc(data.corners*data.face_count,
                                sizeof(*(data.face_vertices)));

    for (size_t i=0; i<data.corners*data.face_count;++i) {
        int count = fscanf(meshfile, " %d ", &(data.face_vertices[i]));
        if (count != 1) {
            perror("Reading faces");
        }
        --(data.face_vertices[i]); // File is 1-based, we want 0-based
    }

    fclose(meshfile);
    
    return data;
}

// Distribute the mesh to all processors
void ScatterMesh(meshdata * data) {
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Bcast(&(data->vertex_count),
              1,                     // Number of items
              MPI_UINT64_T,          // Type
              0,                     // Root
              MPI_COMM_WORLD);

    MPI_Bcast(&(data->face_count) , 1 , MPI_UINT64_T , 0 , comm);
    MPI_Bcast(&(data->dimension)  , 1 , MPI_UINT64_T , 0 , comm);
    MPI_Bcast(&(data->corners)    , 1 , MPI_UINT64_T , 0 , comm);

    // Allocate arrays on rank != 0
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank != 0) {
        data->vertex_coords = calloc(data->dimension*data->vertex_count,
                                    sizeof(*(data->vertex_coords)));
        data->face_vertices = calloc(data->corners*data->face_count,
                                    sizeof(*(data->face_vertices)));
    }
    MPI_Bcast(data->vertex_coords , data->dimension*data->vertex_count,
              MPI_DOUBLE , 0 , comm);
    MPI_Bcast(data->face_vertices , data->corners*data->face_count,
              MPI_INT    , 0 , comm);
}

// Manually create section
// (is this handled by the library somewhere?)
void SetupSection(DM mesh, size_t dof) {
    PetscSection section;
    PetscSectionCreate(MPI_COMM_WORLD, &section);

    int chartStart, chartEnd;
    DMPlexGetChart(mesh, &chartStart, &chartEnd);
    PetscSectionSetChart(section, chartStart, chartEnd);

    int depth;
    DMPlexGetDepth(mesh, &depth);
    for (int s=0; s<depth; ++s){
        int start, end;
        DMPlexGetDepthStratum(mesh, s, &start, &end);
        for (int i=start; i<end; ++i){
            PetscSectionSetDof(section,i,dof);
        }
    }
    PetscSectionSetUp(section);
    DMSetDefaultSection(mesh, section);
}

int main(int argc, char ** argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    int mpirank;
    MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);

    meshdata data;
    if (mpirank == 0) {
        // Read in the mesh
        data = FileMesh();
    }
    ScatterMesh(&data);

    log("Vertex count %zu, Face count %zu\n",data.vertex_count, data.face_count);

    // Create a mesh from the list of vertices and faces (undocumented fun fact:
    // this should only be run on one node, else you get N copies of the points)
    DM mesh;
    if (mpirank == 0) {
        DMPlexCreateFromCellList(MPI_COMM_WORLD, 
                                 2,                  // Topological dimension
                                 data.face_count,    // Number of faces
                                 data.vertex_count,  // Number of vertices
                                 data.corners,       // Corners per face
                                 PETSC_FALSE,        // Interoplate faces and edges
                                 data.face_vertices, // Vertices for each cell
                                 data.dimension,     // Vertex spatial dimension
                                 data.vertex_coords, // Vertex co-ordinates
                                 &mesh);
    } else {
        DMPlexCreateFromCellList(MPI_COMM_WORLD,2,0,0,data.corners,PETSC_FALSE,
                                 NULL,data.dimension,NULL,&mesh);
    }

    // Mesh is only on rank 0, split amongst all ranks
    DM newmesh;
    DMPlexDistribute(mesh,NULL,0,&newmesh);
    if (newmesh) {
        mesh = newmesh;
    }   

    size_t dof = 1;
    SetupSection(mesh, dof);

    // Create a field on the mesh
    Vec temperature;
    DMCreateGlobalVector(mesh, &temperature);

    // Create a local field for this rank
    Vec local_temperature;
    DMCreateLocalVector(mesh, &local_temperature);

    // Get and print size info
    int global_temperature_size;
    VecGetSize(temperature, &global_temperature_size);
    log("Global temperature size %d\n", global_temperature_size);
    int local_temperature_size;
    VecGetSize(local_temperature, &local_temperature_size);
    log("Local temperature size %d\n", local_temperature_size);

    // Set local values
    VecSet(local_temperature, mpirank);

    // Scatter
    DMLocalToGlobalBegin(mesh, local_temperature, INSERT_VALUES, temperature);
    DMLocalToGlobalEnd(mesh, local_temperature, INSERT_VALUES, temperature);

    return PetscFinalize();
}
