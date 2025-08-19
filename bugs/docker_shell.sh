#!/bin/bash

# Helper script to interact with the Basilisk Docker container

case "$1" in
    "shell" | "")
        echo "🐳 Opening interactive shell in Basilisk container..."
        echo "💡 Starting in /host (your local files)"
        docker exec -it basilisk-work /bin/bash -c "cd /host && exec /bin/bash"
        ;;
    "run")
        echo "🚀 Running master-mac simulation..."
        docker exec basilisk-work /bin/bash -c "cd /host && export BASILISK=/workspace/basilisk-C/basilisk-source/src && export PATH=\$PATH:\$BASILISK && ./master-mac"
        ;;
    "compile")
        echo "🔨 Recompiling master-slave coupling..."
        docker exec basilisk-work /bin/bash -c "cd /host && export BASILISK=/workspace/basilisk-C/basilisk-source/src && export PATH=\$PATH:\$BASILISK && ./compile.sh"
        ;;
    "status")
        echo "📊 Container status:"
        docker ps | grep basilisk-work
        echo ""
        echo "📁 Files in /host directory (your local files):"
        docker exec basilisk-work ls -la /host/
        echo ""
        echo "📁 Files in /workspace directory (container copy):"
        docker exec basilisk-work ls -la /workspace/basilisk-C/bugs/ 2>/dev/null || echo "Container copy not available"
        ;;
    "stop")
        echo "🛑 Stopping container..."
        docker stop basilisk-work
        ;;
    "start")
        echo "▶️ Starting container..."
        docker start basilisk-work
        ;;
    "logs")
        echo "📜 Container logs:"
        docker logs basilisk-work
        ;;
    *)
        echo "🐳 Basilisk Docker Helper"
        echo ""
        echo "Usage: $0 [command]"
        echo ""
        echo "Commands:"
        echo "  shell     - Open interactive bash shell (default)"
        echo "  run       - Run the master-mac simulation"
        echo "  compile   - Recompile the master-slave coupling"
        echo "  status    - Show container and file status"
        echo "  stop      - Stop the container"
        echo "  start     - Start the container"
        echo "  logs      - Show container logs"
        echo ""
        echo "Examples:"
        echo "  $0                  # Open shell"
        echo "  $0 shell           # Open shell"
        echo "  $0 run             # Run simulation"
        echo "  $0 compile         # Recompile"
        ;;
esac