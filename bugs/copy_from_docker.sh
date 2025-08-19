#!/bin/bash

# Script to copy files from Basilisk Docker container to local

echo "📁 Copying files from Basilisk Docker container..."

# Create results directory if it doesn't exist
mkdir -p docker-results

case "$1" in
    "executable" | "exe")
        echo "📥 Copying executable..."
        docker cp basilisk-work:/host/master-mac ./docker-results/
        echo "✅ master-mac copied to docker-results/"
        ;;
    "movies" | "mp4")
        echo "🎬 Copying movie files..."
        docker cp basilisk-work:/host/vort-master.mp4 ./docker-results/ 2>/dev/null || echo "⚠️  vort-master.mp4 not found"
        docker cp basilisk-work:/host/vort-slave.mp4 ./docker-results/ 2>/dev/null || echo "⚠️  vort-slave.mp4 not found"
        # Copy any PPM files
        docker exec basilisk-work find /workspace/basilisk-C/bugs -name "*.ppm" -exec cp {} /host/docker-results/ \; 2>/dev/null
        echo "✅ Movie files copied (if they exist)"
        ;;
    "output" | "data")
        echo "📊 Copying output data files..."
        docker exec basilisk-work find /workspace/basilisk-C/bugs -name "*.dat" -exec cp {} /host/docker-results/ \; 2>/dev/null
        docker exec basilisk-work find /workspace/basilisk-C/bugs -name "*.log" -exec cp {} /host/docker-results/ \; 2>/dev/null
        docker exec basilisk-work find /workspace/basilisk-C/bugs -name "*.txt" -exec cp {} /host/docker-results/ \; 2>/dev/null
        echo "✅ Data files copied (if they exist)"
        ;;
    "all" | "")
        echo "📦 Copying all simulation results..."
        
        # Copy executable
        docker cp basilisk-work:/host/master-mac ./docker-results/ 2>/dev/null || echo "⚠️  master-mac not found"
        
        # Copy object files
        docker cp basilisk-work:/host/slave.o ./docker-results/ 2>/dev/null || echo "⚠️  slave.o not found"
        
        # Copy any output files from /host (mounted directory)
        docker exec basilisk-work bash -c "cd /host && find . -maxdepth 1 -type f \( -name '*.mp4' -o -name '*.ppm' -o -name '*.dat' -o -name '*.log' -o -name '*.txt' -o -name '*.plt' \) -exec cp {} /host/docker-results/ \;" 2>/dev/null
        
        echo "✅ All available files copied"
        ;;
    "list")
        echo "📋 Files available in container:"
        docker exec basilisk-work ls -la /host/
        ;;
    *)
        echo "📁 Basilisk Docker File Copy Helper"
        echo ""
        echo "Usage: $0 [type]"
        echo ""
        echo "Types:"
        echo "  all         - Copy all simulation results (default)"
        echo "  executable  - Copy just the master-mac executable"  
        echo "  movies      - Copy movie/animation files (*.mp4, *.ppm)"
        echo "  data        - Copy data files (*.dat, *.log, *.txt)"
        echo "  list        - List available files in container"
        echo ""
        echo "Files are copied to: ./docker-results/"
        echo ""
        echo "Examples:"
        echo "  $0              # Copy everything"
        echo "  $0 all          # Copy everything"
        echo "  $0 executable   # Copy just master-mac"
        echo "  $0 movies       # Copy animations"
        echo "  $0 list         # Show what's available"
        ;;
esac

# Show what was copied
if [ "$1" != "list" ] && [ "$1" != "--help" ]; then
    echo ""
    echo "📂 Files in docker-results/:"
    ls -la docker-results/ 2>/dev/null || echo "No files copied or directory doesn't exist"
fi