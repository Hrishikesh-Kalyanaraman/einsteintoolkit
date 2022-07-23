FROM fedora

RUN dnf install -y which vim findutils git gfortran make tar gcc-c++ python3 python3-pip gsl-devel openblas-devel lapack-devel hdf5-devel
RUN pip3 install termcolor scons

WORKDIR /
RUN git clone https://peterdiener@bitbucket.org/peterdiener/selfforce-1d.git SelfForce-1D
WORKDIR /SelfForce-1D
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 0
RUN scons -f SConstruct.fedora
