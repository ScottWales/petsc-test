/**
 * \file    src/benchmark2.c
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
#include <petsc.h>
#include <stdbool.h>
#include <assert.h>

typedef struct {
    DM mesh;
    PetscSection section;
    Vec data; // Global data array
    int vertexStart;
    int vertexEnd;
} Field;

// Create a new section & vector
Field NewField(DM mesh, unsigned dof) {
    Field f;
    f.mesh = mesh;

    const int dim = 2;
    const int field_count = 1;
    int field_components[field_count];
    int field_dofs[field_count*(dim+1)];

    // Setup each field
    for (int f=0; f<field_count; ++f) {

        // Set dof for vertices only (d=0)
        field_dofs[f*(dim+1)+0] = dof;
        for (int d=1; d<= dim; ++d) {
            field_dofs[f*(dim+1)+d] = 0;
        }

        // Only 1 component per field at the moment
        field_components[f] = 1;
    }

    // Create the section
    DMPlexCreateSection(mesh,
                        dim,
                        field_count,
                        field_components,
                        field_dofs,
                        0,
                        NULL,
                        NULL,
                        &(f.section));
    DMSetDefaultSection(mesh, f.section);

    // Pre-fill vertex counts
    DMPlexGetDepthStratum(mesh, 0, &(f.vertexStart), &(f.vertexEnd));

    // Create a distributed vector for the data
    DMGetGlobalVector(mesh, &(f.data));
    VecScale(f.data, NAN);

    return f;
}

// Print field stats
void DescribeField(Field field) {
    // Mesh diagnostics
    int faceStart, faceEnd;
    DMPlexGetDepthStratum(field.mesh, 1, &faceStart, &faceEnd);
    fprintf(stderr, "Mesh faces: %d - %d (%d total)\n",
            faceStart, faceEnd, faceEnd-faceStart);

    int vertexStart, vertexEnd;
    DMPlexGetDepthStratum(field.mesh, 0, &vertexStart, &vertexEnd);
    fprintf(stderr, "Mesh vertices: %d - %d (%d total)\n",
            vertexStart, vertexEnd, vertexEnd-vertexStart);

    // Section diagnostics
    int graphStart, graphEnd;
    PetscSectionGetChart(field.section, &graphStart, &graphEnd);
    fprintf(stderr, "Section graph vertices: %d - %d (%d total)\n",
            graphStart, graphEnd, graphEnd-graphStart);
    int offsetStart, offsetEnd;
    PetscSectionGetOffsetRange(field.section, &offsetStart, &offsetEnd);
    fprintf(stderr, "Section offset range: %d - %d (%d total)\n",
            offsetStart, offsetEnd, offsetEnd-offsetStart);

    // Vector diagnostics
    int localSize, globalSize;
    VecGetLocalSize(field.data, &localSize);
    VecGetSize(field.data, &globalSize);
    fprintf(stderr, "Data vector local size: %d\n",localSize);
    fprintf(stderr, "Data vector global size: %d\n",globalSize);

}

int main(int argc, char ** argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);

    bool diagnostics = true;

    int overlap = 0;
    DM mesh = ReadMesh(comm,"CubedSphere.mesh",overlap);

    int dof = 1;
    Field temp = NewField(mesh, dof);

    if (diagnostics){
        DescribeField(temp);
    }

    // Set the local vector values
    double * buffer = NULL;
    int vecsize;
    VecGetArray(temp.data,&buffer);
    VecGetLocalSize(temp.data,&vecsize);
    for (int v = temp.vertexStart; v < temp.vertexEnd; ++v) {
        int dof = 0;
        int offset = 0;
        PetscSectionGetDof(temp.section, v, &dof);
        PetscSectionGetFieldOffset(temp.section, v, 0, &offset);
        
        if (offset >= vecsize) {
            fprintf(stderr,"ERROR at vertex %d: offset %d greater than size %d\n",
                    v,offset,vecsize);
//        assert(offset < vecsize);
        }
        buffer[offset] = v;
    }
    VecRestoreArray(temp.data,&buffer);
//    VecAssemblyBegin(temp.data);
//    VecAssemblyEnd(temp.data);

    if (diagnostics){
        double * buffer = NULL;
        VecGetArray(temp.data,&buffer);

        for (int v = temp.vertexStart; v < temp.vertexEnd; ++v) {
            int dof = 0;
            int offset = 0;
            PetscSectionGetDof(temp.section, v, &dof);
            PetscSectionGetOffset(temp.section, v, &offset);

            fprintf(stdout,"% 9d\t% 9d\t%f\n",
                    offset,v,buffer[offset]);
        }

        VecRestoreArray(temp.data,&buffer);
    }

    return PetscFinalize();
}
