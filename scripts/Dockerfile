# Build with: docker build . -t sapphirepp
#   if you need to configure a Proxy use: 
#   docker build --build-arg https_proxy=http://www-cache.mpi-hd.mpg.de:3128 . -t sapphirepp
# Run with: docker run --rm -ti -v $(pwd):/results sapphirepp bash
FROM dealii/dealii:v9.6.2-noble
USER root
# Install editors, nodejs (for actions/checkout), clang-tidy and clang-format
RUN apt-get update && apt-get install -y \
    nano vim \
    nodejs \
    clang-tidy-19 \
    clang-format-19

# Install doxygen-v1.10.0
RUN apt-get install -y graphviz texlive-bibtex-extra && \
    wget https://github.com/doxygen/doxygen/releases/download/Release_1_10_0/doxygen-1.10.0.linux.bin.tar.gz && \
    tar -xzf doxygen-1.10.0.linux.bin.tar.gz && \
    mv doxygen-1.10.0/bin/doxygen /usr/bin/doxygen && \
    rm -rf doxygen-1.10.0 doxygen-1.10.0.linux.bin.tar.gz

# Install sapphirepp
USER dealii
ARG SAPPHIREPP_VERSION=main
RUN git clone https://github.com/sapphirepp/sapphirepp.git sapphirepp && \
    cd sapphirepp && \
    git checkout $SAPPHIREPP_VERSION
WORKDIR /home/dealii/sapphirepp
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DEXAMPLES=ON -DTESTS=ON -DDOC=ON && \
    make -C build -j 2 && \
    make -C build doc
