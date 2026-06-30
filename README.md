# Quantum Hardware MCP Server

An MCP server that gives AI assistants (Claude, etc.) live data about IBM Quantum computers — queue depths, error rates, coherence times, and more.

Built with the [MCP Python SDK](https://github.com/modelcontextprotocol/python-sdk) and [qiskit-ibm-runtime](https://github.com/Qiskit/qiskit-ibm-runtime).

---

## Tools exposed

| Tool | What it does |
|---|---|
| `list_devices` | All IBM quantum computers you can access + status |
| `get_device_details` | Deep info on one machine: error rates, T1/T2, queue |
| `compare_devices` | Rank machines by CX error, queue depth, or qubit count |
| `queue_status` | Current queue snapshot — useful for picking the shortest wait |

---

## Prerequisites

- Python 3.10 or newer
- An IBM Quantum account (free) — sign up at [quantum.ibm.com](https://quantum.ibm.com)
- Claude Desktop (to use the MCP integration)

---

## Setup

### 1. Clone the repo

```bash
git clone https://github.com/YOUR_USERNAME/quantum-hardware-mcp.git
cd quantum-hardware-mcp
```

### 2. Create a virtual environment

```bash
python3 -m venv venv
source venv/bin/activate       # macOS / Linux
# venv\Scripts\activate        # Windows
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Get your IBM Quantum API token

1. Go to [quantum.ibm.com/account](https://quantum.ibm.com/account)
2. Log in (or create a free account)
3. Copy your API token from the "API token" section

### 5. Create your `.env` file

```bash
cp .env.example .env
```

Open `.env` and replace `your_token_here` with your real token:

```
IBM_QUANTUM_TOKEN=abc123...your_actual_token...xyz
```

> **Important:** `.env` is in `.gitignore`. It will never be committed. Never paste your token directly into `server.py`.

### 6. Test the server runs

```bash
python server.py
```

You should see no errors. Press `Ctrl+C` to stop it.
The server speaks over stdin/stdout (MCP stdio transport), so it won't print anything — a clean exit means it's working.

---

## Connect to Claude Desktop

### 1. Find your Claude Desktop config file

| OS | Path |
|---|---|
| macOS | `~/Library/Application Support/Claude/claude_desktop_config.json` |
| Windows | `%APPDATA%\Claude\claude_desktop_config.json` |

### 2. Add the server entry

Open the file (create it if it doesn't exist) and add:

```json
{
  "mcpServers": {
    "quantum-hardware": {
      "command": "/absolute/path/to/venv/bin/python",
      "args": ["/absolute/path/to/quantum-hardware-mcp/server.py"]
    }
  }
}
```

Replace the paths with your actual paths. To find them:

```bash
# From inside the project folder with venv activated:
which python          # → path to use for "command"
pwd                   # → prefix for "args" path
```

**Example on macOS:**

```json
{
  "mcpServers": {
    "quantum-hardware": {
      "command": "/Users/yourname/quantum-hardware-mcp/venv/bin/python",
      "args": ["/Users/yourname/quantum-hardware-mcp/server.py"]
    }
  }
}
```

### 3. Restart Claude Desktop

Quit and reopen Claude Desktop. You should see "quantum-hardware" in the MCP tools list (hammer icon).

---

## Usage examples

Once connected, ask Claude:

- *"List all IBM quantum computers I can access"*
- *"Get details for ibm_brisbane"*
- *"Which quantum computer has the lowest error rate right now?"*
- *"Which machine has the shortest queue right now?"*
- *"Compare all devices by queue depth"*

---

## Project structure

```
quantum-hardware-mcp/
├── server.py          # MCP server — all tool logic lives here
├── requirements.txt   # Python dependencies
├── .env.example       # Token template (safe to commit)
├── .env               # Your real token (git-ignored, never commit)
├── .gitignore
├── LICENSE            # MIT
└── README.md
```

---

## How it works (quick mental model)

```
Claude Desktop
     │  asks a question
     ▼
MCP Protocol (stdin/stdout)
     │  calls a tool
     ▼
server.py  ──── QiskitRuntimeService ────► IBM Quantum API
     │                                          │
     │  ◄─── live device data ─────────────────┘
     ▼
Claude Desktop (formats and answers)
```

The MCP protocol is like a plugin system: Claude Desktop spawns `server.py` as a child process and sends JSON messages over stdin/stdout. The server responds with tool results. No HTTP server needed.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `IBM_QUANTUM_TOKEN is not set` | Check your `.env` file exists and has the right key |
| `IBMNotAuthorizedError` | Token is wrong or expired — get a fresh one from quantum.ibm.com |
| Server not showing in Claude Desktop | Check the paths in `claude_desktop_config.json` are absolute |
| `compare_devices` is slow | Normal — it makes one API call per device to get calibration data |

---

## License

MIT — see [LICENSE](LICENSE).
