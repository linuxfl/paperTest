# paperTest

1. install some package

sudo apt-get install -y python-pip python-numpy mpich python-dev
sudo pip install mpi4py

2. how to run the application

mpirun -f hostfile -np $numofprocess python application_name
