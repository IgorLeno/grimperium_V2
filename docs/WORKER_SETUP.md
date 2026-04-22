# Worker Setup Guide

This guide describes how to set up a Grimperium compute worker on a fresh Linux machine (native or WSL2).

A worker claims molecule batches from the central server, runs CREST + MOPAC locally, and reports results back. It does **not** run the CLI or train ML models — only the chemistry pipeline.

---

## Prerequisites (all platforms)

- 8 GB RAM recommended (4 GB minimum for small molecules)
- Internet access for downloading tools and packages
- Coordinator must provide: `SERVER_URL`, `WORKER_ID`, `API_TOKEN`

---

## Linux (native)

Run the one-shot setup script from a terminal:

```bash
curl -sSL https://raw.githubusercontent.com/IgorLeno/grimperium_V2/main/scripts/setup_worker.sh | bash
```

The script will:
1. Install system build dependencies via `apt`
2. Install Poetry (Python package manager)
3. Clone the Grimperium repository to `~/grimperium_V2`
4. Install Python dependencies (`poetry install --extras "worker"`)
5. Download and install CREST + MOPAC binaries
6. Install MOPAC's required Qt/X11 system libraries
7. Verify the installation and print a summary

**Post-install verification:**

```bash
cd ~/grimperium_V2
poetry run python -c "from grimperium.worker import WorkerRunner; print('Worker OK')"
crest --version
mopac --version
```

**Configure and launch:**

```bash
cp /path/to/.env ~/grimperium_V2/.env   # get this file from the coordinator
cd ~/grimperium_V2
poetry run python -m grimperium.worker.runner
```

---

## Windows via WSL2

### Prerequisites

- Windows 10 build 19041+ or Windows 11
- PowerShell running as Administrator
- TeamViewer installed (for remote support if needed)

### Step 1 — Verify TeamViewer is set to auto-start

Open PowerShell as Administrator and run:

```powershell
Get-Service -Name TeamViewer | Select-Object Name, Status, StartType
```

`StartType` must be `Automatic`. If it is not:

```powershell
Set-Service -Name TeamViewer -StartupType Automatic
Start-Service -Name TeamViewer
```

### Step 2 — Install WSL2 with Ubuntu

```powershell
wsl --install -d Ubuntu
```

Restart the machine when prompted.

### Step 3 — Run the setup script inside Ubuntu

After reboot, open **Ubuntu** from the Start menu and run:

```bash
curl -sSL https://raw.githubusercontent.com/IgorLeno/grimperium_V2/main/scripts/setup_worker.sh | bash
```

### Step 4 — (Optional) Tune WSL2 memory limits

If the machine has less than 16 GB RAM, create `C:\Users\<your-username>\.wslconfig` with:

```ini
[wsl2]
memory=6GB
processors=4
```

Then restart WSL2:

```powershell
wsl --shutdown
```

### Step 5 — Configure and launch

Back inside Ubuntu:

```bash
# Copy the .env file provided by the coordinator
nano ~/grimperium_V2/.env

# Launch the worker
cd ~/grimperium_V2
poetry run python -m grimperium.worker.runner
```

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `ModuleNotFoundError: rich` | Old poetry lock without the `rich` extra | `poetry install --extras "worker"` |
| MOPAC exits immediately with no output | Missing Qt/X11 system libs | `sudo apt install -y libxcb-icccm4 libxcb-shape0 libxkbcommon-x11-0` |
| `crest: command not found` after reboot | PATH not reloaded | `source ~/.bashrc` |
| `poetry: command not found` | PATH not reloaded | `source ~/.bashrc` |
| Worker can't reach server | Firewall or wrong SERVER_URL | Confirm `SERVER_URL` in `.env`, check VPN/firewall |

---

## Updating an existing worker

```bash
cd ~/grimperium_V2
git pull origin main
poetry install --extras "worker"
```

No need to re-run `install_tools.sh` unless the CREST or MOPAC version has changed.
