FROM basilisk-minimal

RUN rm -f /etc/apt/apt.conf.d/docker-clean

RUN sed -i \
      's/^Components: main$/Components: main contrib non-free non-free-firmware/' \
      /etc/apt/sources.list.d/debian.sources
      
RUN apt-get update && \
    apt-get install -y \
      -o APT::Keep-Downloaded-Packages=true \
      libglfw3-dev \
      nvidia-cuda-dev \
      hipcc \
      ocl-icd-opencl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /basilisk
