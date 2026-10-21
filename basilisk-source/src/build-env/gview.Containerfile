FROM basilisk-minimal

RUN rm -f /etc/apt/apt.conf.d/docker-clean

RUN apt-get update && \
    apt-get install -y \
      -o APT::Keep-Downloaded-Packages=true \
      g++ libglfw3-dev libimgui-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /basilisk
