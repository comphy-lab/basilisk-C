# Tips and Troubleshooting

This document contains tips and solutions for common issues encountered when working with the Basilisk-C repository.

## Installing Darcs

Darcs is not available in the Ubuntu package repositories for newer Ubuntu versions. Instead, you can install it using Cabal:

```shell
# First install Cabal
sudo apt-get update
sudo apt-get install -y ghc cabal-install

# Then install Darcs
cabal update
cabal install darcs
```

After installation, make sure to add the Cabal binary directory to your PATH:

```shell
echo "$HOME/.cabal/bin" >> ~/.bashrc
source ~/.bashrc
```

## Installing OpenConnect VPN for Discoverer HPC on Ubuntu 24.04

If you encounter issues installing OpenConnect VPN on Ubuntu 24.04 (Noble), follow these steps:

1. Edit `/etc/apt/sources.list.d/ubuntu.sources` and replace "noble" with "jammy" at two instances:
   - Change `Suites: noble noble-updates noble-backports` to `Suites: jammy noble-updates noble-backports`
   - Change `Suites: noble-security` to `Suites: jammy-security`

2. Update and upgrade your system:
   ```shell
   sudo apt-get update
   sudo apt-get upgrade
   ```

3. Remove mailcap:
   ```shell
   sudo apt autoremove
   ```

4. Force install the required dependency:
   ```shell
   sudo apt install libtss2-mu0=3.2.0-1ubuntu1.1
   ```

5. Install OpenConnect:
   ```shell
   sudo apt-get install globalprotect-openconnect
   ```

6. **Important**: After installing OpenConnect, you may need to switch the sources back to "noble" to avoid dependency issues with other packages (like python3).

## GitHub Actions Workflow for Darcs Synchronization

We've set up a GitHub Actions workflow to automatically sync this repository with the Darcs repository. The workflow:

1. Runs daily at midnight UTC
2. Can be triggered manually
3. Installs Darcs using Cabal
4. Pulls all changes from the Darcs repository
5. Commits and pushes any changes to GitHub

The workflow file is located at `.github/workflows/sync-darcs.yml`.