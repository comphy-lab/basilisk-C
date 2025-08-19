# Updated Docker Workflow (Host-Based Development)

## 🎯 **Live Sync Development Model**

**Key Principle**: Edit on Mac, compile in Docker, outputs sync automatically

### **Directory Structure**

| Location | Purpose | Sync Status |
|----------|---------|-------------|
| **`/host`** | Your Mac `bugs/` directory | ✅ Live bidirectional sync |
| **`/workspace`** | Container Basilisk copy | ❌ Static (for compilation only) |
| **`/project`** | Your Mac `basilisk-C/` parent | ✅ Live bidirectional sync |

### **Recommended Workflow**

#### 1. **Edit Files Anywhere**
```bash
# On your Mac (any editor)
vim master-mac.c
vim slave-mac.c
```

#### 2. **Compile in Docker**  
```bash
./docker_shell.sh compile    # Uses /host directory automatically
```

#### 3. **Run Simulations**
```bash
./docker_shell.sh run        # Runs from /host directory
```

#### 4. **Access Results Immediately**
```bash
ls -la *.mp4                 # Available on Mac immediately!
ls -la master-mac slave.o    # Compiled files visible locally
```

## 🔄 **No Copy Step Needed!**

**OLD workflow:**
```bash
# Edit locally → Compile in Docker → Copy results back ❌
./copy_from_docker.sh
```

**NEW workflow:**
```bash
# Edit locally → Compile in Docker → Results already local! ✅
ls -la *.mp4   # Videos are already in your local directory
```

## 🐚 **Interactive Shell**

```bash
./docker_shell.sh           # Opens shell in /host directory
cd /host                     # Already here!
pwd                          # Shows: /host
ls                          # Shows your local Mac files
./compile.sh                 # Works directly
./master-mac                 # Generates videos in local directory
exit                        # Videos remain on your Mac
```

## 📊 **Status Monitoring**

```bash
./docker_shell.sh status    # Shows both /host and /workspace
```

Output example:
```
📁 Files in /host directory (your local files):
-rwxr-xr-x 1 root root 525720 master-mac
-rw-r--r-- 1 root root 235582 vort-master.mp4
-rw-r--r-- 1 root root 327522 vort-slave.mp4

📁 Files in /workspace directory (container copy):
-rw-r--r-- 1 root root 3292 master-mac.c
-rw-r--r-- 1 root root 1405 slave-mac.c
```

## ⚡ **Key Benefits**

1. **Real-time sync**: Edit in any Mac editor, see changes instantly in Docker
2. **No file management**: No copying back and forth
3. **All outputs local**: Videos, executables, data files appear on Mac automatically  
4. **Fast iteration**: Edit → Compile → Run → Results (no copy step)
5. **Version control ready**: All changes tracked in git on Mac side

## 🎬 **Video Generation**

Videos are generated directly in your Mac directory:
```bash
./docker_shell.sh run       # Run simulation
ls -la vort-*.mp4           # Videos appear locally immediately
open vort-master.mp4        # View with Mac video player
```

## 🔧 **Development Commands**

| Command | Action | Location |
|---------|--------|----------|
| `./docker_shell.sh` | Interactive shell | `/host` (your Mac files) |
| `./docker_shell.sh compile` | Recompile | `/host` using `compile.sh` |
| `./docker_shell.sh run` | Run simulation | `/host` |
| `./docker_shell.sh status` | Show file status | Both locations |

## ✅ **Migration Complete**

- ✅ All scripts updated to use `/host` directory
- ✅ Video generation working without warnings  
- ✅ Live sync enabled for seamless development
- ✅ No manual file copying required
- ✅ Local git tracking of all changes

**Result**: True hybrid Mac/Linux development with zero friction! 🚀