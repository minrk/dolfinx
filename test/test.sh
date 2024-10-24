set -ex
env | grep PATH
env | grep PREFIX
# copy activation from impi_rt activate.bat
export I_MPI_ROOT="$CONDA_PREFIX\Library"
export PYTHONUNBUFFERED=1
export "PATH=$I_MPI_ROOT\bin\libfabric\utils:$I_MPI_ROOT\bin\libfabric:$PATH"
echo PATH=$PATH
export I_MPI_DEBUG=1000
which mpiexec
which -a python
mpiexec --version
mpiexec -n 1 python -c 'from mpi4py import MPI; print(MPI.COMM_WORLD.rank)'
mpiexec -n 2 python -c 'from mpi4py import MPI; print(MPI.COMM_WORLD.rank)'
python -c 'from mpi4py import MPI; print(MPI.COMM_WORLD.rank)'
mpiexec -n 1 python test/test_mesh.py
python test/test_dolfinx.py
mpiexec -v -np 1 python -c "from mpi4py import MPI; import dolfinx; dolfinx.mesh.create_rectangle(comm=MPI.COMM_WORLD, points=((0, 0), (2, 1)), n=(32, 16))"
mpiexec -v -np 2 python -c "from mpi4py import MPI; import dolfinx; dolfinx.mesh.create_rectangle(comm=MPI.COMM_WORLD, points=((0, 0), (2, 1)), n=(32, 16))"
