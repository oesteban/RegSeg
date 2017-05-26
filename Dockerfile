# Use Ubuntu Trusty LTS
FROM ubuntu:trusty-20170330

# Pre-cache neurodebian key
COPY .docker/neurodebian.gpg /root/.neurodebian.gpg

# Prepare environment
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
                    curl \
                    bzip2 \
                    ca-certificates \
                    xvfb \
                    pkg-config && \
    curl -sSL http://neuro.debian.net/lists/trusty.us-ca.full >> /etc/apt/sources.list.d/neurodebian.sources.list && \
    apt-key add /root/.neurodebian.gpg && \
    (apt-key adv --refresh-keys --keyserver hkp://ha.pool.sks-keyservers.net 0xA5D32F012649A5A9 || true) && \
    apt-get update

# Installing Ubuntu packages and cleaning up
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
                    libinsighttoolkit4-dev=4.7.0-1~nd14.04+1 \
                    cmake=2.8.12.2-0ubuntu3 \
                    g++ \
                    build-essential \
                    libjsoncpp-dev \
                    libvtk6-dev \
                    libvtkgdcm2-dev \
                    libboost-filesystem-dev \
                    libboost-system-dev \
                    libboost-program-options-dev \
                    libfftw3-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY ./Code /root/regseg/src
RUN mkdir /root/regseg/Release && \
    cd /root/regseg/Release && \
    cmake ../src/ -G"Unix Makefiles"  -DCMAKE_BUILD_TYPE=Release -DITK_DIR=/usr/local/lib/cmake/ITK-4.7/ && \
    make -j$( grep -c ^processor /proc/cpuinfo ) && \
    make install

ENTRYPOINT ["/usr/local/bin/regseg"]

# Store metadata
ARG BUILD_DATE
ARG VCS_REF
ARG VERSION
LABEL org.label-schema.build-date=$BUILD_DATE \
      org.label-schema.name="RegSeg" \
      org.label-schema.description="RegSeg -" \
      org.label-schema.url="https://github.com/oesteban/RegSeg" \
      org.label-schema.vcs-ref=$VCS_REF \
      org.label-schema.vcs-url="https://github.com/oesteban/RegSeg" \
      org.label-schema.version=$VERSION \
      org.label-schema.schema-version="1.0"

# Installing and setting up miniconda
# RUN curl -sSLO https://repo.continuum.io/miniconda/Miniconda3-4.3.11-Linux-x86_64.sh && \
#     bash Miniconda3-4.3.11-Linux-x86_64.sh -b -p /usr/local/miniconda && \
#     rm Miniconda3-4.3.11-Linux-x86_64.sh

# ENV PATH=/usr/local/miniconda/bin:$PATH \
#     LANG=C.UTF-8 \
#     LC_ALL=C.UTF-8

# # Installing precomputed python packages
# RUN conda install -c conda-forge -y openblas=0.2.19; \
#     sync && \
#     conda install -c conda-forge -y \
#                      numpy=1.12.0 \
#                      scipy=0.19.0 \
#                      scikit-learn=0.18.1 \
#                      matplotlib=2.0.0 \
#                      pandas=0.19.2 \
#                      libxml2=2.9.4 \
#                      libxslt=1.1.29 \
#                      sympy=1.0 \
#                      statsmodels=0.8.0 \
#                      dipy=0.11.0 \
#                      traits=4.6.0 \
#                      psutil=5.2.2 \
#                      sphinx=1.5.4; \
#     sync &&  \
#     chmod +x /usr/local/miniconda/bin/* && \
#     conda clean --all -y; sync && \
#     python -c "from matplotlib import font_manager"
