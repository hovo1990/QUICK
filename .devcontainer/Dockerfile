############################
### Base MPI CUDA 12.0.1 ###
### -- ! For QUICK compatibility
############################
FROM nvidia/cuda:12.0.1-devel-ubuntu22.04 AS base-mpi-cuda-12.0.1

# FROM  nvidia/cuda:12.6.0-devel-ubuntu24.04 AS base-mpi-cuda-12.6.0

RUN apt-get update -y \
 && apt-get install -y \
    gfortran \
    gdb \
    valgrind \
    binutils \
    # -- ! cmake this version does not have vs code cmake debugging capabilities
    g++ \
    git \
    openmpi-bin \
    openmpi-common \
    libopenmpi-dev \
    bash \
    python3 \
    python3-pip \
    python3-pydot \ 
    ca-certificates \
    gpg \
    wget \
    build-essential \
    unzip \
    graphviz \
    # -- ! https://forums.developer.nvidia.com/t/vscode-debug-matrixmul-example-fails-with-cuda-driver-has-hit-an-internal-error/239437/2
    libcudadebugger1


# -- * Install newest cmake from official kitware apt repository
RUN test -f /usr/share/doc/kitware-archive-keyring/copyright || \
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /usr/share/keyrings/kitware-archive-keyring.gpg >/dev/null

RUN echo 'deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ jammy main' | tee /etc/apt/sources.list.d/kitware.list >/dev/null 
RUN apt-get update


RUN test -f /usr/share/doc/kitware-archive-keyring/copyright || \
    rm /usr/share/keyrings/kitware-archive-keyring.gpg
    

RUN  apt-get update && apt-get install kitware-archive-keyring

RUN apt-get -y  install cmake

# -- * For ubuntu 24.04
# RUN rm /usr/lib/python3.12/EXTERNALLY-MANAGED



RUN pip3 install fortls fprettify


#-- * ubuntu 22.04 cmake version is cmake version 3.22.1
#-- * ubuntu 24.04 cmake version is cmake version 3.28.3


# -- * Create a new user
RUN useradd -ms /bin/bash vscode

# -- * Switch to the new user
USER vscode


# # Manually run steps from quick.rc
# ENV QUICK_INSTALL /src/install
# ENV QUICK_BASIS $QUICK_INSTALL/basis
# ENV PATH $PATH:$QUICK_INSTALL/bin

