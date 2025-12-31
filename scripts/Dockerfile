# Build with:
#   docker build . -t sapphirepp
# If you need to configure a Proxy use: 
#   docker build --build-arg https_proxy=http://$https_proxy . -t sapphirepp
# Run with: 
#   docker run --rm -ti -u root -v $(pwd)/results:/home/dealii/sapphirepp/results sapphirepp bash
FROM dealii/dealii:v9.7.1-noble
USER root
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
