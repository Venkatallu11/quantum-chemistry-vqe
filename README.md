# Quantum Hardware MCP Server

An MCP server that gives AI assistants (Claude, etc.) live data about IBM Quantum computers — queue depths, error rates, coherence times, best qubit picks, and historical snapshots.

Built with the [MCP Python SDK](https://github.com/modelcontextprotocol/python-sdk) and [qiskit-ibm-runtime](https://github.com/Qiskit/qiskit-ibm-runtime).

---

## Why researchers should care

Running quantum experiments is expensive in two ways: **time** (queue wait) and **quality** (gate errors corrupt results). This server exposes both dimensions as AI-queryable tools, so you can ask natural questions and get grounded answers:

- **Live qubit picks** — `best_qubits` ranks every qubit on a device right now by calibration data, so you can hand the compiler a pre-filtered register instead of hoping the transpiler gets lucky.
- **Historical lookup** — `device_on_date` retrieves the hardware snapshot for any past date from the local database, so your methods section can cite what the machine actually looked like when you ran the experiment.
- **Reproducibility** — snapshots are recorded every 6 hours automatically. If a reviewer asks "what was the CX error on the day you ran Figure 3?", you can answer exactly.
- **Instant comparison** — `compare_devices(sort_by="combined")` blends gate quality and queue depth into a single score, surfacing the best machine to use *right now*.

---

## Tools exposed

| Tool | What it does |
|---|---|
| `list_devices` | All IBM quantum computers you can access + status |
| `get_device_details` | Deep info on one machine: error rates, T1/T2, queue |
| `compare_devices` | Rank machines by CX error, queue depth, qubit count, or combined score |
| `queue_status` | Current queue snapshot — useful for picking the shortest wait |
| `device_history` | Snapshots for one machine over the last N days |
| `best_qubits` | Best n qubits on a machine right now, scored by calibration data |
| `device_on_date` | Historical stats for a machine on any past date (reproducibility) |

### `compare_devices` sort modes

| `sort_by` | What it optimises |
|---|---|
| `cx_error` | Lowest 2-qubit gate error — highest fidelity results |
| `queue` | Fewest pending jobs — fastest turnaround |
| `qubits` | Most qubits — largest circuits |
| `combined` | 70% quality + 30% availability, min-max normalised across current devices |

---

## Prerequisites

- Python 3.10 or newer
- An IBM Quantum account (free) — sign up at [quantum.ibm.com](https://quantum.ibm.com)
- Claude Desktop (to use the MCP integration)

---

## Setup

### 1. Clone the repo

```bash
git clone https://github.com/Lokesh-2025/quantum-hardware-mcp.git
cd quantum-hardware-mcp
```

### 2. Create a virtual environment and install

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 3. Add your IBM Quantum token

```bash
cp .env.example .env
```

Open `.env` and replace `your_token_here` with your token from [quantum.ibm.com/account](https://quantum.ibm.com/account).

> `.env` is in `.gitignore` — it will never be committed.

### 4. Connect to Claude Desktop

Edit `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "quantum-hardware": {
      "command": "/absolute/path/to/.venv/bin/python",
      "args": ["/absolute/path/to/quantum-hardware-mcp/server.py"]
    }
  }
}
```

Quit and reopen Claude Desktop. The hammer icon will show the quantum-hardware tools.

---

## Automatic snapshots

A background agent (`snapshot.py`) records device stats every 6 hours:

- **Locally** — a macOS LaunchAgent writes to `devices.db`, which feeds `device_history` and `device_on_date`.
- **GitHub Actions** — a scheduled workflow appends rows to `data/snapshots.csv` every 6 hours, building a public historical record.

---

## Project structure

```
quantum-hardware-mcp/
├── server.py          # MCP server — all 7 tools live here
├── snapshot.py        # Background agent — records device stats every 6h
├── report.py          # Daily Quantum Weatherman report (runs at 8am)
├── requirements.txt
├── .env.example
├── .github/
│   └── workflows/
│       └── snapshot.yml   # GitHub Actions: snapshot → CSV every 6h
├── data/
│   └── snapshots.csv      # Growing historical record (committed by CI)
├── reports/               # Daily reports + charts (git-ignored)
├── REPORTS.md             # Running summary log
└── devices.db             # Local SQLite snapshot store (git-ignored)
```

---

## Planned roadmap

- **Trend alerts** — notify when a device's error rate spikes above its 7-day average
- **Queue forecasting** — predict wait time from historical queue patterns
- **Multi-vendor** — add IonQ and Quantinuum device data alongside IBM
- **Circuit-aware ranking** — given a circuit's gate profile, recommend the best machine for that specific workload
- **Automated report push** — opt-in posting of the daily Weatherman to Slack / Discord

---

## License

MIT — see [LICENSE](LICENSE).
