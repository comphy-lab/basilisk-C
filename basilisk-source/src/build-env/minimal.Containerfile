FROM docker.io/library/debian:stable

RUN rm -f /etc/apt/apt.conf.d/docker-clean

RUN apt-get update && \
    apt-get install -y \
      -o APT::Keep-Downloaded-Packages=true \
      gcc make darcs \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /basilisk
