FROM registry.gitlab.com/nsf-muses/common:dev AS build
############################################################################################
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -yq \
    build-essential cmake git pkg-config \
    openmpi-bin libopenmpi-dev \
    libhdf5-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Copy all the current dir 
ENV KOKKOSSRC=/tmp/kokkos
ENV CABANASRC=/tmp/Cabana
ENV PREFIX=/usr/local
WORKDIR /tmp
COPY . /tmp/

# Install Kokkos and Cabana
RUN git clone https://github.com/kokkos/kokkos.git

WORKDIR "$KOKKOSSRC"
RUN pwd
RUN git checkout 4.1.00 && \
    mkdir -p "$KOKKOSSRC/build"
    
WORKDIR "$KOKKOSSRC/build"
RUN pwd

RUN cmake -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DCMAKE_CXX_COMPILER=g++ \
      -DCMAKE_CXX_STANDARD=17 \
      -DCMAKE_CXX_EXTENSIONS=Off \
      -DCMAKE_BUILD_TYPE="Release" \
      -DKokkos_ENABLE_COMPILER_WARNINGS=ON \
      -DKokkos_ENABLE_CUDA=Off \
      -DKokkos_ENABLE_CUDA_LAMBDA=Off \
      -DKokkos_ENABLE_OPENMP=On \
      -DKokkos_ENABLE_SERIAL=On \
      -DKokkos_ENABLE_TESTS=Off \
      -DKokkos_ARCH_AMPERE80=Off "$KOKKOSSRC"

RUN cmake --build . --target install -j

WORKDIR /tmp
RUN git clone https://github.com/ECP-copa/Cabana.git "$CABANASRC"
WORKDIR "$CABANASRC" 
RUN git checkout 2b642f8dbba760e2008f8b81a4196e6adf100f42 && \
    mkdir -p "$CABANASRC/build"

WORKDIR "$CABANASRC/build"
RUN cmake -DCMAKE_BUILD_TYPE="Release" \
      -DCMAKE_CXX_STANDARD=17 \
      -DCMAKE_PREFIX_PATH="/usr/local;/usr" \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DCMAKE_CXX_COMPILER=g++ \
      -DCabana_REQUIRE_CUDA=OFF \
      -DCabana_ENABLE_TESTING=OFF \
      -DCabana_ENABLE_EXAMPLES=OFF "$CABANASRC"

RUN cmake --build . --target install -j

WORKDIR "$CABANASRC/build"
WORKDIR /tmp
RUN cmake -S . -B build -DCMAKE_PREFIX_PATH="$PREFIX"
WORKDIR /tmp/build
RUN make


FROM registry.gitlab.com/nsf-muses/common:dev AS deps
#############################
# Install git for pip install
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    python3 python3-pip python3-venv git \
 && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
# COPY requirements.txt /opt/
# RUN pip install --no-cache-dir -r /opt/requirements.txt


FROM registry.gitlab.com/nsf-muses/common:dev AS final
#####################
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -yq \
        make python3 python3-pip \
        libhdf5-103-1 libhdf5-cpp-103-1 \
        openmpi-bin libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

ARG USERNAME=CCAKE
ARG UID=1000
## Create user with desired UID and set "home" to /opt
RUN useradd --uid $UID --non-unique --no-log-init --no-create-home --home-dir /opt --shell /bin/bash $USERNAME
RUN chown -R $UID:$UID /opt
USER $USERNAME
## Install compiled executable and add to PATH
# COPY --from=deps /opt/venv /opt/venv
# ENV PATH="/opt/venv/bin:/opt/src:/opt/.local/bin:${PATH}"
## Install Python dependencies# Install source code with compiled output
WORKDIR /opt
COPY --from=build --chown=$UID:$UID /usr/local /usr/local
COPY --from=build --chown=$UID:$UID /tmp/ /opt/ 
# COPY --from=build --chown=$UID:$UID /tmp/dependencies dependencies
# COPY --chown=$UID:$UID api api
# COPY --chown=$UID:$UID input input
# RUN rm -rf src/*.cpp src/*.hpp
## remove build from container image
# RUN rm -rf src/build
## Copy manifest, Dockerfile, Readme.md, LICENSE, and *.sh files
# COPY --chown=$UID:$UID manifest.yaml ./
COPY --chown=$UID:$UID Dockerfile ./
# COPY --chown=$UID:$UID LICENSE ./
# COPY --chown=$UID:$UID Readme.md ./
## create input and output directories
# RUN mkdir -p /opt/input /opt/output
# change to source code directory
WORKDIR /opt/
