#!/bin/bash

# Basilisk Master-Slave Coupling Compilation Script
# This script compiles the master-slave coupling example with proper environment setup

echo "🔧 Basilisk Master-Slave Coupling Compiler"
echo "==========================================="

# Set up Basilisk environment
export BASILISK=/workspace/basilisk-C/basilisk-source/src
export PATH=$PATH:$BASILISK:/usr/local/bin

# Ensure we're in the correct directory context
# Check if we're in the mounted /host directory (preferred) or container copy
if [ -d "/host" ] && [ -f "/host/master-mac.c" ]; then
    cd /host
    echo "📍 Using mounted directory: /host"
elif [ -d "/workspace/basilisk-C/bugs" ] && [ -f "/workspace/basilisk-C/bugs/master-mac.c" ]; then
    cd /workspace/basilisk-C/bugs
    echo "📍 Using container directory: /workspace/basilisk-C/bugs"
else
    echo "❌ Error: Cannot find source files in expected locations"
    exit 1
fi

echo "📍 Working directory: $(pwd)"
echo "🛠️  BASILISK path: $BASILISK"
echo "🔍 qcc location: $(which qcc)"
echo "🎬 ppm2mp4 location: $(which ppm2mp4 2>/dev/null || echo \"not found\")"
echo "🎞️  ffmpeg location: $(which ffmpeg)"
echo ""

# Check if source files exist
if [ ! -f "master-mac.c" ]; then
    echo "❌ Error: master-mac.c not found!"
    exit 1
fi

if [ ! -f "slave-mac.c" ]; then
    echo "❌ Error: slave-mac.c not found!"
    exit 1
fi

echo "✅ Source files found"
echo ""

# Ensure Basilisk tools are built and available
if [ ! -f "$BASILISK/qcc" ]; then
    echo "🔨 qcc not found, building Basilisk first..."
    cd "$BASILISK"
    
    # Set up config if needed
    if [ ! -f "config" ] && [ -f "config.gcc" ]; then
        echo "🔗 Setting up config.gcc -> config"
        ln -sf config.gcc config
    fi
    
    # Build dependencies first
    echo "📦 Building AST library..."
    make -C ast libast.a || {
        echo "❌ Failed to build AST library"
        exit 1
    }
    
    echo "📦 Building include and postproc..."
    make include.o postproc.o || {
        echo "❌ Failed to build include.o postproc.o"
        exit 1
    }
    
    # Now build qcc
    echo "📦 Building qcc compiler..."
    make qcc || {
        echo "❌ Failed to build qcc"
        exit 1
    }
    
    # Make video tools executable
    chmod +x ppm2mp4 ppm2mpeg ppm2ogv ppm2gif 2>/dev/null || true
    
    # Return to work directory
    cd - > /dev/null
    
    echo "✅ Basilisk tools built successfully"
fi

# Clean previous build
echo "🧹 Cleaning previous build..."
rm -f slave.o master-mac

# Step 1: Compile slave to object file
echo "📝 Step 1: Compiling slave-mac.c to object file..."
qcc -O2 -fno-common -D_OBJECT -c slave-mac.c -o slave.o

if [ $? -ne 0 ]; then
    echo "❌ Failed to compile slave-mac.c"
    exit 1
fi

echo "✅ slave.o created"

# Step 2: Filter symbols with objcopy (this is the key step that works on Linux)
echo "🔧 Step 2: Filtering symbols with GNU objcopy..."
objcopy -G slave_step -G slave_stop -G slave_interpolate slave.o

if [ $? -ne 0 ]; then
    echo "❌ Failed to filter symbols with objcopy"
    exit 1
fi

echo "✅ Symbols filtered (kept: slave_step, slave_stop, slave_interpolate)"

# Step 3: Compile master and link with filtered slave object
echo "🔗 Step 3: Linking master-mac.c with filtered slave.o..."
qcc -O2 master-mac.c slave.o -o master-mac -lm

if [ $? -ne 0 ]; then
    echo "❌ Failed to link master-mac"
    exit 1
fi

echo "✅ master-mac executable created"
echo ""

# Verify the executable
if [ -f "master-mac" ]; then
    echo "🎉 SUCCESS! Compilation completed successfully."
    echo ""
    echo "📊 File information:"
    ls -la master-mac slave.o
    echo ""
    echo "🚀 To run the simulation:"
    echo "   ./master-mac"
    echo ""
    echo "⏱️  To run for limited time:"
    echo "   timeout 30s ./master-mac"
else
    echo "❌ ERROR: master-mac executable not created"
    exit 1
fi
