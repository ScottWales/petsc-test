/**
 * \file    src/benchmark.c
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

#include <petsc.h>
#include <assert.h>
#include <stdlib.h>

// Write to stderr with MPI rank & function name
#define log(fmt, ...) \
    do { \
        int _msize, _mrank; \
        MPI_Comm_rank(MPI_COMM_WORLD,&_mrank); \
        MPI_Comm_size(MPI_COMM_WORLD,&_msize); \
        fprintf(stderr, "[%d/%d] %s:\t" fmt, _mrank, _msize, __FUNCTION__, __VA_ARGS__); \
    } while (0)

// Read in a mesh from a file
DM NewMesh(MPI_Comm comm, unsigned overlap) {
    FILE * input = fopen("CubedSphere.mesh","r");

    // Read in vertex co-ordinates
    size_t dim = 3;
    size_t vertex_count;
    fscanf(input, " %zu ", &vertex_count);
    assert(vertex_count == 1572866);
    double * vertices = calloc(vertex_count*dim, sizeof(*vertices));
    for (size_t vertex=0; vertex<vertex_count; ++vertex) {
        for (int d=0; d<dim; ++d) {
            int read = fscanf(input, " %lf ", &vertices[vertex*dim+d]);
            assert(read == 1);
        }
    }
    assert(vertices[0] == -1.00000000000000000);
    assert(vertices[3] == -0.99999529380957619);

    // Read in each face's vertex id
    size_t corners = 4;
    size_t face_count;
    fscanf(input, " %zu ", &face_count);
    assert(face_count == 1572864);
    int * faces = calloc(face_count*corners, sizeof(*faces));
    for (size_t face=0; face<face_count; ++face) {
        for (int c=0; c<corners; ++c) {
            int read = fscanf(input, " %d ", &faces[corners*face+c]);
            assert(read == 1);

            // Vertex ids are 0-based in petsc
            --faces[corners*face+c];
        }
    }
    assert(faces[0] == 1232838-1);
    assert(faces[4] == 1233870-1);
    fclose(input);

    int rank;
    MPI_Comm_rank(comm,&rank);

    DM mesh = NULL;
    if (rank == 0) {
        DMPlexCreateFromCellList(comm,         // MPI Communicator
                                 2,            // Spatial dimension of graph
                                 face_count,   // Number of faces
                                 vertex_count, // Number of vertices
                                 corners,      // Number of vertices per face
                                 PETSC_FALSE,  // Create edges
                                 faces,        // Vertex IDs for each face
                                 dim,          // Vertex dimension
                                 vertices,     // Co-ordinates of each vertex
                                 &mesh);       // Output
    } else {
        DMPlexCreateFromCellList(comm, 2, 0, 0, corners, PETSC_FALSE,
                                 NULL, dim, NULL, &mesh);
    }

    DM distmesh = NULL;
    DMPlexDistribute(mesh, NULL, overlap, &distmesh);
    if (distmesh) {
        mesh = distmesh;
    }

    free(faces);
    free(vertices);

    return mesh;
}

void test(DM mesh) {
    const int dim = 2;
    const int field_count = 1;
    int components[field_count];
    int dof[field_count*(dim+1)];
    
    // Field 0 - temperature
    components[0]      =  1;
    dof[0*(dim+1) + 0] =  1; // vertices
    dof[0*(dim+1) + 1] =  0; // faces

    // Field 1 - pressure
    components[1]      =  1;
    dof[1*(dim+1) + 0] = 64; // vertices
    dof[1*(dim+1) + 1] =  0; // faces

    PetscSection section;
    DMPlexCreateSection(mesh,
                        dim,         // Spatial dimension
                        field_count, // Number of fields
                        components,  // Component count for each field
                        dof,         // DoF for each field component
                        0,           // Number of boundary conditions
                        NULL,        // BC Field
                        NULL,        // BC Points
                        &section);
    DMSetDefaultSection(mesh, section);

    // Temperature
    // Create a global vector
    Vec temperature = NULL;
    DMCreateGlobalVector(mesh, &temperature);

    // Get a temporary local vector & a pointer to its data
    Vec localTemp = NULL;
    DMGetLocalVector(mesh, &localTemp);
    PetscScalar * data = NULL;
    VecGetArray(localTemp, &data);

    // Loop over points in the section
    int start, end;
    DMGetDefaultSection(mesh, &section);
    PetscSectionGetChart(section, &start, &end);
    for (int point=start; point<end; ++point) {
        // Get the size & offset of field 0
        int dof, offset;
        PetscSectionGetFieldDof(section, point, 0, &dof);
        PetscSectionGetFieldOffset(section, point, 0, &offset); 

        // Set values
        for (int d=0; d<dof; ++d){
            data[offset+d] = point * 100 + d;            
        }
    }
    VecRestoreArray(localTemp, &data);

    // Broadcast values
    DMLocalToGlobalBegin(mesh, localTemp, INSERT_VALUES, temperature);
    DMLocalToGlobalEnd(mesh, localTemp, INSERT_VALUES, temperature);
    // Now each node should have all temperature values

    // Free the temp vector
    DMRestoreLocalVector(mesh,&localTemp);
}

// Create a section where each vertex has the given degrees of freedom
PetscSection NewSection(MPI_Comm comm, DM mesh, uint32_t dof) {
    PetscSection section;
    PetscSectionCreate(comm, &section);

    // Set the section size from the mesh
    int chartStart, chartEnd;
    DMPlexGetChart(mesh, &chartStart, &chartEnd);
    PetscSectionSetChart(section, chartStart, chartEnd);

    // Set the dof for each vertex
    // Can also do faces/edges depending on the integer argument to
    // GetDepthStratum
    int vertexStart, vertexEnd;
    DMPlexGetDepthStratum(mesh, 0, &vertexStart, &vertexEnd);
    for (int vertex = vertexStart; vertex<vertexEnd; ++vertex) {
        PetscSectionSetDof(section, vertex, dof);
    }

    // Finalise setup
    PetscSectionSetUp(section);
    return section;
}

// Create a field on the vertices with given degrees of freedom
Vec NewField(MPI_Comm comm, DM mesh, uint32_t dof) {
    Vec field;
    VecCreate(comm, &field);

    // Create a new section
    PetscSection section = NewSection(comm, mesh, dof);

    // Set the vector size
    int size;
    PetscSectionGetStorageSize(section, &size);
    VecSetSizes(field, size, PETSC_DETERMINE);

    // Finalise setup
    VecSetFromOptions(field);
    return field;
}

void DescribeVec(Vec v) {
    int low,high,size,total;
    VecGetLocalSize(v, &size);
    VecGetSize(v, &total);
    VecGetOwnershipRange(v, &low, &high);
    log("Local size %d/%d\n",size,total);
    log("Ownership range %d - %d\n",low,high);
}


int main(int argc, char ** argv) {
    // Setup library (similar to MPIInit)
    PetscInitialize(&argc,&argv,
                    NULL,              // 'database' (config?) file
                    "UKMO Benchmark"); // Help message

    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);

    // Mesh setup
    unsigned overlap = 1;
    DM mesh = NewMesh(comm, overlap);

    // The co-ordinate field
    Vec coords;
    DMGetCoordinatesLocal(mesh, &coords);

    test(mesh);

    // Create a 1-Dof field
    Vec temperature = NewField(comm, mesh, 1);
    DescribeVec(temperature);

    // Set values
    PetscScalar * data = NULL;
    VecGetArray(temperature, &data);

    int low,high;
    VecGetOwnershipRange(temperature, &low, &high);
    for (int vertex = low; vertex<high; ++vertex) {
        data[vertex] = vertex*100 + rank;
    }
    VecRestoreArray(temperature, &data);

    // Update vector
    VecAssemblyBegin(temperature);
    VecAssemblyEnd(temperature);

    // Create a 64-Dof field
    Vec pressure = NewField(comm, mesh, 64);

    // Exchange

    // Cleanup
    return PetscFinalize();
}
