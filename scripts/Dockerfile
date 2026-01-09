# Build with:
#   docker build . -t sapphirepp
# If you need to configure a Proxy use: 
#   docker build --build-arg https_proxy=http://$https_proxy . -t sapphirepp
# Run with: 
#   docker run --rm -ti -u root -v $(pwd)/results:/home/dealii/sapphirepp/results sapphirepp bash
FROM dealii/dealii:v9.7.1-noble
USER root
# Allow run mpi as root (workaround for host filesystem owner matching problem)
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ARG NJOBS=2      # Jobs used for building
# Install editors, nodejs (for actions/checkout), clang-tidy and clang-format
RUN apt-get update && apt-get install -y \
    nano vim \
    nodejs \
    clang-tidy-19 \
    clang-format-19

# Install doxygen-v1.15.0 (compile from source needed for arm64)
RUN apt-get install -y graphviz texlive-bibtex-extra flex bison
RUN wget https://github.com/doxygen/doxygen/releases/download/Release_1_15_0/doxygen-1.15.0.src.tar.gz && \
    tar -xzf doxygen-1.15.0.src.tar.gz && \
    cmake -S doxygen-1.15.0 -B doxygen-build && \
    make -C doxygen-build -j $NJOBS install && \
    rm -rf doxygen-build doxygen-1.15.0 doxygen-1.15.0.src.tar.gz

# Install sapphirepp
USER dealii
ARG SAPPHIREPP_VERSION=main
RUN git clone https://github.com/sapphirepp/sapphirepp.git sapphirepp && \
    cd sapphirepp && \
    git checkout $SAPPHIREPP_VERSION
WORKDIR /home/dealii/sapphirepp
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DEXAMPLES=ON -DTESTS=ON -DDOC=ON && \
    make -C build -j $NJOBS && \
    make -C build doc

# Install miniforge and sapphireppplot
ENV PYTHONPATH=
RUN cd .. && \
    wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash Miniforge3.sh -b -p $HOME/miniforge3 && \
    rm Miniforge3.sh
ENV PATH=$HOME/miniforge3/bin:$PATH
RUN conda config --set auto_activate false && \
    conda init bash

ARG SAPPHIREPPPLOT_VERSION=main
WORKDIR /home/dealii
RUN git clone https://github.com/sapphirepp/sapphireppplot.git sapphireppplot && \
    cd sapphireppplot && \
    git checkout $SAPPHIREPPPLOT_VERSION
WORKDIR /home/dealii/sapphireppplot
RUN conda env create -f environment.yml
RUN export PATH=$HOME/miniforge3/envs/sapplot/bin:$PATH && \
    pip install -e '.[dev,docs]'
WORKDIR /home/dealii/sapphirepp
