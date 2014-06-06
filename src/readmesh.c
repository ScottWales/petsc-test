/**
 * \file    src/readmesh.c
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

#include "readmesh.h"

// Read in a mesh from a file
DM ReadMesh(MPI_Comm comm, const char * filename, unsigned overlap) {
    FILE * input = fopen(filename,"r");

    // Read in vertex co-ordinates
    size_t dim = 3;
    size_t vertex_count;
    fscanf(input, " %zu ", &vertex_count);
    double * vertices = calloc(vertex_count*dim, sizeof(*vertices));
    for (size_t vertex=0; vertex<vertex_count; ++vertex) {
        for (int d=0; d<dim; ++d) {
            int read = fscanf(input, " %lf ", &vertices[vertex*dim+d]);
        }
    }

    // Read in each face's vertex id
    size_t corners = 4;
    size_t face_count;
    fscanf(input, " %zu ", &face_count);
    int * faces = calloc(face_count*corners, sizeof(*faces));
    for (size_t face=0; face<face_count; ++face) {
        for (int c=0; c<corners; ++c) {
            int read = fscanf(input, " %d ", &faces[corners*face+c]);

            // Vertex ids are 0-based in petsc
            --faces[corners*face+c];
        }
    }
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
