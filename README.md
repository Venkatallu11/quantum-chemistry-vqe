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
| ---- | ------------ |
| `list_devices` | All IBM quantum computers you can access + status |
| `get_device_details` | Deep info on one machine: error rates, T1/T2, queue |
| `compare_devices` | Rank machines by CX error, queue depth, qubit count, or combined score |
| `queue_status` | Current queue snapshot — useful for picking the shortest wait |
| `device_history` | Snapshots for one machine over the last N days |
| `best_qubits` | Best n qubits on a machine right now, scored by calibration data |
| `device_on_date` | Historical stats for a machine on any past date (reproducibility) |
| `submit_job` | Compile and submit an OpenQASM 2.0 or 3.0 circuit to IBM hardware — returns a `job_id` |
| `job_status` | Check status of a submitted job (QUEUED / RUNNING / DONE / ERROR) |
| `job_results` | Retrieve bit-string measurement counts from a completed job |
| `cancel_job` | Cancel a queued or running job by `job_id` |
| `list_jobs` | List your most recent jobs with status and backend (newest first) |
| `run_grover` | Built-in Grover's search demo — builds the circuit, picks the least-busy backend, submits, returns `job_id` |
| `estimate_expectation` | Run the Estimator primitive — computes expectation values ⟨ψ\|O\|ψ⟩ for Pauli observables (needed for VQE, QAOA, quantum chemistry) |
| `circuit_report` | Dry-run analysis — transpiles your circuit and returns gate counts, qubit mapping, per-pair CX errors, and an estimated fidelity. No queue time. |

### `compare_devices` sort modes

| `sort_by` | What it optimises |
| --------- | ----------------- |
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

## Docker

### MCP server only

Runs just `server.py` in HTTP mode on port 3020. Use this if you want to connect the server to an existing agent or to Claude Desktop over HTTP.

**1. Build:**
```bash
docker build -t quantum-mcp .
```

**2. Run:**
```bash
docker run -p 3020:3020 -e IBM_QUANTUM_TOKEN=your_token quantum-mcp
```

The server will be available at `http://localhost:3020/sse`.

---

### Fullstack (MCP server + agent)

Runs both `server.py` and Jack's Node.js agent together. The agent connects to the server automatically — no manual config needed.

**1. Create a `.env` file in the project root:**
```bash
IBM_QUANTUM_TOKEN=your_ibm_token
LLM_PROVIDER=anthropic
ANTHROPIC_API_KEY=your_anthropic_key
```

**2. Build and start both containers:**
```bash
docker compose up --build
```

- MCP server: `http://localhost:3020`
- Agent API: `http://localhost:3021`

**3. Chat with the agent:**

In a separate terminal, run `node agent/chat.js` (or `npm run chat` from the `agent/` directory) — it will connect to the agent running in Docker.

**To stop:**
```bash
docker compose down
```

The SQLite snapshot database is persisted in a Docker volume (`mcp-data`) so historical data survives container restarts.

---

## HTTP Server Mode (Optional)

The MCP server can also run as an HTTP/SSE server for remote access or web-based AI assistants.

### Quick Start

```bash
# Start HTTP server on localhost (development)
python server.py --transport http

# Start on all interfaces for remote access
python server.py --transport http --host 0.0.0.0 --port 8080

# With custom CORS origins
python server.py --transport http --cors-origins "https://myapp.com,https://api.myapp.com"
```

### Configuration

HTTP server settings can be configured via command-line arguments or environment variables:

| Setting | CLI Argument | Environment Variable | Default |
| ------- | ------------ | -------------------- | ------- |
| Host | `--host` | `MCP_HTTP_HOST` | `127.0.0.1` |
| Port | `--port` | `MCP_HTTP_PORT` | `8000` |
| CORS Origins | `--cors-origins` | `MCP_CORS_ORIGINS` | `*` |

### Security: API Key Authentication

For production deployments, enable API key authentication:

**1. Generate a secure API key:**

```bash
python -c "import secrets; print(secrets.token_urlsafe(32))"
```

**2. Add to your `.env` file:**

```bash
MCP_API_KEY=your_generated_key_here
```

**3. Start the server:**

```bash
python server.py --transport http --host 0.0.0.0 --port 8080
```

**4. Include the key in client requests:**

Python:

```python
import requests

headers = {
    "X-API-Key": "your_generated_key_here",
    "Content-Type": "application/json"
}

response = requests.post(
    "http://your-server:8080/sse",
    headers=headers
)
```

JavaScript:

```javascript
fetch('http://your-server:8080/sse', {
  method: 'POST',
  headers: {
    'X-API-Key': 'your_generated_key_here',
    'Content-Type': 'application/json'
  }
})
```

MCP Client Configuration:

```json
{
  "mcpServers": {
    "quantum-hardware": {
      "url": "http://your-server:8080",
      "headers": {
        "X-API-Key": "your_generated_key_here"
      }
    }
  }
}
```

### Development vs Production

**Development Mode (No Authentication):**

- If `MCP_API_KEY` is not set, the server runs without authentication
- Convenient for local testing
- ⚠️ Only use on localhost, never expose to the internet

**Production Mode (With Authentication):**

- Set `MCP_API_KEY` environment variable
- Prevents unauthorized access
- Always use HTTPS in production (set up reverse proxy with nginx/caddy)
- Rotate API keys periodically

### Deployment Examples

#### Docker Deployment

See the **Docker** section above for `Dockerfile`, `agent/Dockerfile`, and `docker-compose.yml`.

#### Reverse Proxy (nginx)

```nginx
server {
    listen 443 ssl;
    server_name quantum-mcp.example.com;

    ssl_certificate /path/to/cert.pem;
    ssl_certificate_key /path/to/key.pem;

    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }
}
```

---

## Automatic snapshots

A background agent (`snapshot.py`) records device stats every 6 hours:

- **Locally** — a macOS LaunchAgent writes to `devices.db`, which feeds `device_history` and `device_on_date`.
- **GitHub Actions** — a scheduled workflow appends rows to `data/snapshots.csv` every 6 hours, building a public historical record.

---

## Project structure

```text
quantum-hardware-mcp/
├── server.py          # MCP server — all 13 tools live here
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
├── reports/fleet-report/  # Polished fleet analysis report + charts
├── test_bell_state.py     # End-to-end Bell state test on real hardware
├── REPORTS.md             # Running summary log
└── devices.db             # Local SQLite snapshot store (git-ignored)
```

---

## Running quantum circuits

Use the three job tools together to submit a circuit and get results:

```python
# 1. Pick the least-busy machine
queue_status()                            # → ibm_kingston (0 jobs)

# 2. Submit an OpenQASM 2.0 circuit
submit_job("ibm_kingston", """
OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg c[2];
h q[0];
cx q[0],q[1];
measure q[0] -> c[0];
measure q[1] -> c[1];
""", shots=1024)
# → { "job_id": "d8tvs9lbh0os73eq93ag", "status": "QUEUED", ... }

# 3. Check progress
job_status("d8tvs9lbh0os73eq93ag")        # → RUNNING → DONE

# 4. Get results
job_results("d8tvs9lbh0os73eq93ag")
# → { "counts": {"00": 509, "11": 487, "01": 26, "10": 2}, "total_shots": 1024 }
```

The circuit above is a Bell state — the ~97% correlation between `00` and `11` outcomes is quantum entanglement measured on real hardware. Run `test_bell_state.py` for an automated end-to-end test.

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
