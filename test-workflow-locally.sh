#!/bin/bash
# Test the Darcs sync workflow locally using Docker
# This simulates the GitHub Actions workflow to verify HTTPS connectivity
#
# First run: builds a Docker image with Darcs pre-installed (~10-15 min)
# Subsequent runs: uses cached image (seconds)
#
# Usage:
#   ./test-workflow-locally.sh          # Use cached image or build if needed
#   ./test-workflow-locally.sh --rebuild # Force rebuild the image

set -e

IMAGE_NAME="darcs-test"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=============================================="
echo "Testing Darcs Sync Workflow Locally (Docker)"
echo "=============================================="
echo ""

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo "Error: Docker is not running. Please start Docker and try again."
    exit 1
fi

echo "Docker is running."

# Check if we need to build the image
BUILD_IMAGE=false
if [[ "$1" == "--rebuild" ]]; then
    BUILD_IMAGE=true
    echo "Forcing image rebuild..."
elif ! docker image inspect "$IMAGE_NAME" > /dev/null 2>&1; then
    BUILD_IMAGE=true
    echo "Image '$IMAGE_NAME' not found. Building (this takes ~10-15 min first time)..."
fi

if [[ "$BUILD_IMAGE" == true ]]; then
    echo ""
    echo "Building Docker image with Darcs pre-installed..."
    echo "This only needs to happen once."
    echo ""
    docker build -f "$SCRIPT_DIR/Dockerfile.darcs-test" -t "$IMAGE_NAME" "$SCRIPT_DIR"
    echo ""
    echo "Image built successfully! Future runs will be fast."
    echo ""
fi

echo ""
echo "Running Darcs sync test..."
echo ""

# Run the test using the pre-built image
docker run --rm "$IMAGE_NAME"

echo ""
echo "Local test completed successfully!"
exit 0
