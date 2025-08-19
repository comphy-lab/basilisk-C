#!/bin/bash

# Helper script to build and run Docker container for Basilisk compilation

echo "🔨 Building Docker container for Basilisk..."
docker build -t basilisk-dev .

if [ $? -ne 0 ]; then
    echo "❌ Docker build failed"
    exit 1
fi

echo "🚀 Starting Docker container..."
echo "📁 Your current directory will be mounted at /host in the container"
echo "📂 Parent directory will be mounted at /project (for basilisk-source access)"
echo "🔧 Basilisk environment will be set up automatically on container start"

# Get the absolute paths
BUGS_DIR="$(pwd)"
PROJECT_DIR="$(dirname "$(pwd)")"
BASILISK_SOURCE_DIR="$PROJECT_DIR/basilisk-source"

# Check if basilisk-source exists
if [ ! -d "$BASILISK_SOURCE_DIR" ]; then
    echo "⚠️  Warning: $BASILISK_SOURCE_DIR not found"
    echo "    Make sure you're in the basilisk-C/bugs directory and basilisk-source exists"
fi

# Run container with proper mounts
docker run -it --rm \
    -v "$BUGS_DIR":/host \
    -v "$PROJECT_DIR":/project \
    -v "$BASILISK_SOURCE_DIR":/workspace/basilisk-source \
    --name basilisk-work \
    basilisk-dev