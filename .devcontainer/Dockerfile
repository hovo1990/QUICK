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
    bash


