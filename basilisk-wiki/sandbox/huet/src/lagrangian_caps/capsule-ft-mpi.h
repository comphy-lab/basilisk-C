/**
# MPI synchronization of the Lagrangian mesh

In the current implementation, each processor has its own copy of the Lagrangian
mesh. In order to update the position of all the local copies of the Lagrangian
mesh, the processors need to communicate to each other the contributions of
their flow field to the velocities of each Lagrangian nodes.
*/

void reduce_lagVel(lagMesh* mesh) {
  int li = dimension; // length of one item
  double* send_data = (double*)malloc(mesh->nln*li*sizeof(double));
  double* recv_data = (double*)malloc(mpi_npe*mesh->nln*li*sizeof(double));
  for(int i=0; i<mesh->nln; i++) {
    send_data[li*i] = mesh->nodes[i].lagVel.x;
    send_data[li*i+1] = mesh->nodes[i].lagVel.y;
    #if dimension > 2
    send_data[li*i+2] = mesh->nodes[i].lagVel.z;
    #endif
  }
  MPI_Allgather(send_data, mesh->nln*li, MPI_DOUBLE, recv_data,
    mesh->nln*li, MPI_DOUBLE, MPI_COMM_WORLD);

  for(int i=0; i<mesh->nln; i++) {
    mesh->nodes[i].lagVel.x = 0.;
    mesh->nodes[i].lagVel.y = 0.;
    #if dimension > 2
    mesh->nodes[i].lagVel.z = 0.;
    #endif
  }

  for(int k=0; k<mpi_npe; k++) {
    for(int i=0; i<mesh->nln; i++) {
      mesh->nodes[i].lagVel.x += recv_data[k*li*mesh->nln+li*i];
      mesh->nodes[i].lagVel.y += recv_data[k*li*mesh->nln+li*i+1];
      #if dimension > 2
      mesh->nodes[i].lagVel.z += recv_data[k*li*mesh->nln+li*i+2];
      #endif
    }
  }
  free(send_data);
  free(recv_data);
}
