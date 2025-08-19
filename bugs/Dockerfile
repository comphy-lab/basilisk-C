# Docker container for Basilisk compilation on Linux
FROM ubuntu:22.04

# Set non-interactive mode for apt
ENV DEBIAN_FRONTEND=noninteractive

# Install essential build tools and Basilisk dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    make \
    gawk \
    gnuplot \
    imagemagick \
    ffmpeg \
    graphviz \
    valgrind \
    git \
    darcs \
    binutils \
    vim \
    curl \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Create a working directory
WORKDIR /workspace

# Set up environment variables that Basilisk might need
ENV BASILISK=/workspace/basilisk-source/src
ENV PATH=$PATH:$BASILISK

# Create a bootstrap script that will set up Basilisk properly
# This script will run when the container starts and ensure video tools are available
RUN echo '#!/bin/bash\n\
echo "🔧 Setting up Basilisk environment..."\n\
\n\
# Check if Basilisk source is available (should be mounted)\n\
if [ -d "/workspace/basilisk-source/src" ]; then\n\
    cd /workspace/basilisk-source/src\n\
    \n\
    # Set up config file for build\n\
    if [ ! -f "config" ]; then\n\
        if [ -f "config.gcc" ]; then\n\
            echo "🔗 Linking config.gcc -> config"\n\
            ln -sf config.gcc config\n\
        fi\n\
    fi\n\
    \n\
    # Build Basilisk if not already built\n\
    if [ ! -f "qcc" ] || [ ! -x "ppm2mp4" ]; then\n\
        echo "📦 Building Basilisk (this may take a moment)..."\n\
        make -j$(nproc) qcc 2>/dev/null || {\n\
            echo "⚠️  Initial qcc build failed, trying basic make..."\n\
            make qcc 2>/dev/null || true\n\
        }\n\
        \n\
        # Ensure video tools are executable\n\
        make ppm2mp4 ppm2mpeg ppm2ogv ppm2gif 2>/dev/null || {\n\
            echo "🎬 Making video tools executable manually..."\n\
            chmod +x ppm2mp4 ppm2mpeg ppm2ogv ppm2gif 2>/dev/null || true\n\
        }\n\
    fi\n\
    \n\
    # Create symlinks in /usr/local/bin for global access\n\
    echo "🔗 Creating tool symlinks..."\n\
    for tool in ppm2mp4 ppm2mpeg ppm2ogv ppm2gif; do\n\
        if [ -f "/workspace/basilisk-source/src/$tool" ]; then\n\
            ln -sf "/workspace/basilisk-source/src/$tool" "/usr/local/bin/$tool" 2>/dev/null || true\n\
            echo "  ✓ $tool -> /usr/local/bin/$tool"\n\
        else\n\
            echo "  ⚠ $tool not found"\n\
        fi\n\
    done\n\
    \n\
    # Special handling for qcc\n\
    if [ -f "/workspace/basilisk-source/src/qcc" ]; then\n\
        ln -sf "/workspace/basilisk-source/src/qcc" "/usr/local/bin/qcc" 2>/dev/null || true\n\
        echo "  ✓ qcc -> /usr/local/bin/qcc"\n\
    else\n\
        echo "  ⚠ qcc not found - Basilisk build may have failed"\n\
        echo "  ℹ This is expected if basilisk-source is not mounted or needs building"\n\
    fi\n\
    \n\
    echo "✅ Basilisk environment setup complete"\n\
    echo "🎬 Video tools: $(ls /usr/local/bin/ppm2* 2>/dev/null | wc -l) available"\n\
    echo "🎞️  ffmpeg: $(which ffmpeg)"\n\
    echo "🔧 qcc: $(which qcc 2>/dev/null || echo \"not found - will be built during compilation\")"\n\
else\n\
    echo "⚠️  Basilisk source not found at /workspace/basilisk-source/src"\n\
    echo "   Make sure volume mounts are properly configured"\n\
fi\n\
\n\
# Execute the passed command\n\
exec "$@"\n' > /usr/local/bin/basilisk-init.sh && \
    chmod +x /usr/local/bin/basilisk-init.sh

# Set the init script as the entrypoint
ENTRYPOINT ["/usr/local/bin/basilisk-init.sh"]

# Default command
CMD ["/bin/bash"]