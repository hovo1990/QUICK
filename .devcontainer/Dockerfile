############################
### Base MPI CUDA 12.0.1 ###
############################
FROM nvidia/cuda:12.0.1-devel-ubuntu22.04 AS base-mpi-cuda-12.0.1

RUN apt-get update -y \
 && apt-get install -y \
    gfortran \
    cmake \
    g++ \
    git \
    openmpi-bin \
    openmpi-common \
    libopenmpi-dev \
    bash \
    python3 \
    python3-pip

RUN pip3 install fortls


# Create a new user
RUN useradd -ms /bin/bash vscode

# Switch to the new user
USER vscode

