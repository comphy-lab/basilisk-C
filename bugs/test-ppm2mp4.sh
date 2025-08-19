#!/bin/bash

# Test script to verify ppm2mp4 functionality in Docker
# This creates a simple PPM file and converts it to MP4

echo "🧪 Testing ppm2mp4 functionality..."

# Create a simple test PPM file (100x100 red square)
cat > test.ppm << 'EOF'
P3
100 100
255
EOF

# Add red pixels (100x100 = 10000 pixels)
for i in {1..10000}; do
    echo "255 0 0" >> test.ppm
done

echo "✅ Created test.ppm (100x100 red square)"

# Test ppm2mp4 directly
echo "🎬 Testing ppm2mp4 conversion..."

# Create a series of PPM files for video
for i in {1..5}; do
    cp test.ppm "frame-$(printf %04d $i).ppm"
done

echo "✅ Created 5 frame files"

# Test ppm2mp4 conversion  
echo "🎞️  Converting to MP4..."
if command -v ppm2mp4 >/dev/null; then
    # Concatenate all PPM frames and pipe to ppm2mp4 (this is how Basilisk uses it)
    cat frame-*.ppm | ppm2mp4 test-output.mp4
    if [ -f "test-output.mp4" ]; then
        echo "✅ ppm2mp4 conversion successful!"
        ls -la test-output.mp4
    else
        echo "❌ ppm2mp4 conversion failed - no output file"
    fi
else
    echo "❌ ppm2mp4 not found in PATH"
fi

# Test ffmpeg directly
echo "🎥 Testing ffmpeg directly..."
if command -v ffmpeg >/dev/null; then
    ffmpeg -f image2 -i frame-%04d.ppm -r 25 -c:v libx264 -crf 18 direct-ffmpeg.mp4 -y 2>/dev/null
    if [ -f "direct-ffmpeg.mp4" ]; then
        echo "✅ Direct ffmpeg conversion successful!"
        ls -la direct-ffmpeg.mp4
    else
        echo "❌ Direct ffmpeg conversion failed"
    fi
else
    echo "❌ ffmpeg not found in PATH"
fi

echo ""
echo "🔍 Tool locations:"
echo "   ppm2mp4: $(which ppm2mp4 2>/dev/null || echo 'not found')"
echo "   ffmpeg: $(which ffmpeg 2>/dev/null || echo 'not found')"
echo "   avconv: $(which avconv 2>/dev/null || echo 'not found')"

# Clean up
rm -f test.ppm frame-*.ppm

echo ""
echo "🎯 Test complete!"