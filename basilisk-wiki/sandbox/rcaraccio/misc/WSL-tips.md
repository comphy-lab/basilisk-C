# WSL tips

This file collects a handful of tips I've found while running this sandbox under WSL2
(Ubuntu on Windows). Some are general WSL gotchas — disk usage, core allocation — and
some are specific to Basilisk workflows like the `jview` interactive viewer or MPI
builds. It's an informal, growing scratchpad rather than a polished guide, and I'll try
to keep adding to it as I hit new issues.

## "My Ubuntu install takes 140 GB but `du` only shows 16 GB"

Two independent effects stack up here.

### 1. `du -hs *` massively undercounts

Run in your home directory, `du -hs *` has three blind spots:

- **Home directory only** — it never sees `/usr`, `/var`, `/tmp`, `/opt`, etc.
- **The `*` glob skips dotfiles** — `.local`, `.cache`, `.vscode-server`, ... are invisible.
- Even within home, subtotals hide behind the glob.

Get the real numbers instead:

```bash
df -h /                                   # actual usage of the whole distro
du -hx --max-depth=1 / 2>/dev/null | sort -hr | head -20   # per top-level dir
du -hs ~/.[!.]* ~/* 2>/dev/null | sort -hr | head          # home incl. dotfiles
```

### 2. The WSL2 virtual disk never shrinks on its own

WSL2 stores the whole distro in an `ext4.vhdx` file on Windows. That file **grows on
demand but never shrinks automatically** — once usage peaks (big builds, downloads,
Docker layers) the vhdx stays that large even after you delete the files. So Windows
reports the peak size, not current usage. The gap = reclaimable slack.

## Reclaiming space

### Step 1 — free real space inside Linux

```bash
sudo apt clean                        # wipe downloaded .deb cache (often GBs)
sudo journalctl --vacuum-size=100M    # trim systemd logs
rm -rf /tmp/*                         # safe
du -hx ~/.cache /var/cache | sort -hr | head   # find other offenders
```

### Step 2 — compact the vhdx (this is what shrinks the Windows file)

Freeing space inside Linux does **nothing** to the Windows file size until you compact.
From **PowerShell as Administrator**:

```powershell
wsl --shutdown
wsl -l -v                             # find the exact distro name
wsl --manage {Distro} --set-sparse true
```

`--set-sparse true` also makes the vhdx auto-shrink going forward, so this stops
recurring. To compact manually instead:

```powershell
wsl --shutdown
diskpart
# in diskpart:
select vdisk file="C:\Users\{you}\AppData\Local\Packages\{UbuntuPkg}\LocalState\ext4.vhdx"
attach vdisk readonly
compact vdisk
detach vdisk
exit
```

## `jview` interactive window fails to connect

When opening the Basilisk interactive viewer with `jview`, the page defaults to the
WSL hostname, which the browser on Windows can't resolve. Change the WebSocket URL from:

```
ws://LAPTOP-Rcaraccio.localdomain:7100
```

to:

```
ws://localhost:7100
```

(WSL2 forwards `localhost` to the distro, but not the `.localdomain` hostname.)

If this still fails, try switching browser. 

## MPI runs limited by WSL core count

If an MPI build (`CC='mpicc -D_MPI={n}' make {case}.tst`) can't use as many ranks as
expected, WSL may be capped below your physical core count. Check what WSL sees:

```bash
nproc                                 # cores currently available to WSL
```

The cap lives in `.wslconfig` in your Windows user profile — from WSL it's at
`/mnt/c/Users/Riccardo/.wslconfig` (i.e. `C:\Users\Riccardo\.wslconfig` on Windows):

```ini
[wsl2]
memory=8GB
processors=8    # <-- raise to the number of cores you want WSL to use
```

Edit `processors` (and `memory` if needed), then restart WSL from **PowerShell** for
it to take effect:

```powershell
wsl --shutdown
```

Don't exceed your machine's physical core count. Check total cores on Windows with
`echo %NUMBER_OF_PROCESSORS%` (cmd) or in Task Manager.



