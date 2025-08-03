# Build with: docker build . -t sapphirepp/sapphirepp:v1.2.0
# Run with: docker run --rm -ti -v $(pwd):/results sapphirepp/sapphirepp:v1.2.0 bash
FROM dealii/dealii:v9.6.0-jammy
USER root
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    python3 \
    python3-pip \
    libopenmpi-dev \
    libgsl-dev \
    libboost-all-dev \
    libhdf5-openmpi-dev \
    liblapack-dev \
    libblas-dev \
    libgmp-dev \
    libmpfr-dev \
    libmuparser-dev \
    libtbb-dev \
    libassimp-dev \
    curl

# Install nodejs-v18
RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - && \
    apt-get -y install nodejs 

# Install doxygen-v1.10.0
RUN sudo apt-get install -y cmake graphviz texlive-bibtex-extra && \
    wget https://github.com/doxygen/doxygen/releases/download/Release_1_10_0/doxygen-1.10.0.linux.bin.tar.gz && \
    tar -xzf doxygen-1.10.0.linux.bin.tar.gz && \
    mv doxygen-1.10.0/bin/doxygen /usr/bin/doxygen && \
    rm -rf doxygen-1.10.0 doxygen-1.10.0.linux.bin.tar.gz

# Add LLVM APT repo and install clang-tidy-19
RUN wget https://apt.llvm.org/llvm.sh && \
    chmod +x llvm.sh && \
    ./llvm.sh 19 && \
    rm llvm.sh && \
    apt-get install -y clang-tidy-19

# Install sapphirepp
USER dealii
# RUN git clone https://github.com/sapphirepp/sapphirepp.git && \
#   cd sapphirepp && \
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DDOC=ON -DCMAKE_CXX_CLANG_TIDY="clang-tidy-19;--use-color" && \
#   make -C build -j 2 && \
#   make -C build doc
