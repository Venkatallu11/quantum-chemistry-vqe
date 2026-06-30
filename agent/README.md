# Quantum Hardware MCP Agent

This Node.js application serves as an intelligent agent that leverages Large Language Models (LLMs) to interact with a `quantum-hardware-mcp-server` instance. The agent is designed to understand user questions about quantum hardware, select the appropriate tool from the `quantum-hardware-mcp-server`, execute it with the correct arguments, and then format the results into a natural, human-readable answer.

## Features

- **Official MCP SDK Integration:** Uses `@modelcontextprotocol/sdk` for proper MCP protocol communication
- **SSE Transport:** Communicates with quantum-hardware-mcp server via Server-Sent Events (SSE)
- **Multiple LLM Provider Support:** Choose from Google Gemini, Ollama (local), OpenAI, Anthropic Claude, or vLLM
- **Intelligent Tool Selection:** Utilizes LLMs to analyze user queries and select the most appropriate tool from the available tools provided by the `quantum-hardware-mcp-server`
- **Interactive Query Capability:** 🆕 Agent can ask clarifying questions during execution without exiting the ReAct loop (see Query Capability section)
- **Fallback to Direct Answer:** If the model cannot determine an appropriate tool to use, it will provide a direct answer based on its general knowledge
- **Dynamic Tool Execution:** Calls the selected tool on the `quantum-hardware-mcp-server` with the necessary arguments
- **Natural Language Response:** Formats the JSON or structured data returned by the tools into a clear and understandable natural language response
- **Chat History Context:** Considers conversation history to provide more relevant and accurate answers
- **MCP API Key Support:** Supports authentication with quantum-hardware-mcp server via X-API-Key header

## Prerequisites

- Node.js (v18 or higher) and npm installed
- Access to a running `quantum-hardware-mcp-server` instance
- An LLM provider (choose one):
  - **Google Gemini** - Requires API key (paid service)
  - **Ollama** - Free, runs locally, requires Ollama installation
  - **OpenAI** - Requires API key (paid service)
  - **Anthropic Claude** - Requires API key (paid service)
  - **vLLM** - Free, runs locally, high-performance inference engine

## Choosing an LLM Provider

### Google Gemini

- **Best for:** Production deployments, high accuracy
- **Pros:** Latest models, reliable API, good performance
- **Cons:** Requires paid API key, data sent to Google
- **Setup:** Get API key from [Google AI Studio](https://makersuite.google.com/app/apikey)

### Ollama (Recommended for Privacy)

- **Best for:** Local development, privacy-sensitive data, cost savings
- **Pros:** Free, runs locally, no API costs, keeps data private
- **Cons:** Requires local resources, model must be downloaded
- **Setup:** See [Ollama Setup](#ollama-setup) below

### OpenAI

- **Best for:** Advanced reasoning, GPT-4 capabilities
- **Pros:** Powerful models, reliable API
- **Cons:** Higher costs, data sent to OpenAI
- **Setup:** Get API key from [OpenAI Platform](https://platform.openai.com/api-keys)

### Anthropic Claude

- **Best for:** Long context, detailed analysis
- **Pros:** Excellent reasoning, large context window
- **Cons:** Requires paid API key, data sent to Anthropic
- **Setup:** Get API key from [Anthropic Console](https://console.anthropic.com/)

### vLLM

- **Best for:** High-performance local inference, production deployments
- **Pros:** Fast inference, efficient memory usage, runs locally, supports many models
- **Cons:** Requires GPU, more complex setup than Ollama
- **Setup:** See [vLLM Setup](#vllm-setup) below

## Installation

1. Navigate to the agent directory:

   ```bash
   cd /home/jwoehr/work/AI/MCP/quantum-hardware-mcp/agent
   ```

2. Install all dependencies (includes all LLM provider packages):

   ```bash
   npm install
   ```

   This installs:
   - Core dependencies (@modelcontextprotocol/sdk, express, cors, dotenv)
   - All LLM provider packages (@google/generative-ai, ollama, openai, @anthropic-ai/sdk)
   - Shared infrastructure (providers, config, concurrency)

3. Copy the example environment file and configure it:

   ```bash
   cp .env.example .env
   ```

4. Edit `.env` and configure your chosen provider (see [Configuration](#configuration) below)

## Configuration

### Basic Configuration

Create a `.env` file based on `.env.example`:

```bash
# Choose your LLM provider
LLM_PROVIDER=gemini  # Options: gemini, ollama, openai, anthropic, vllm

# Quantum Hardware MCP Server (Required)
# FastMCP default endpoint is /sse
QUANTUM_MCP_SERVER_URI=http://127.0.0.1:3020/sse

# MCP API Key (Required if server has authentication enabled)
MCP_API_KEY=your_mcp_api_key_here
```

### Provider-Specific Configuration

#### Gemini Configuration

```bash
LLM_PROVIDER=gemini
GEMINI_API_KEY=your_api_key_here
GEMINI_MODEL=gemini-2.5-pro
```

Available models:

- `gemini-2.5-pro` - Latest and most capable
- `gemini-1.5-pro` - Previous generation, still excellent
- `gemini-1.5-flash` - Faster, lower cost

#### Ollama (Local)

```bash
LLM_PROVIDER=ollama
OLLAMA_BASE_URL=http://localhost:11434
OLLAMA_MODEL=llama3.1:8b
```

Recommended models for quantum hardware queries:

- `llama3.1:8b` - Best balance of speed and quality (4.7GB)
- `llama3.1:70b` - Highest quality, slower (40GB)
- `mistral:7b` - Fast responses (4.1GB)
- `codellama:13b` - Better for technical queries (7.4GB)
- `qwen2.5-coder:7b` - Optimized for code (4.7GB)

#### OpenAI Configuration

```bash
LLM_PROVIDER=openai
OPENAI_API_KEY=your_api_key_here
OPENAI_MODEL=gpt-4o
```

Available models:

- `gpt-4o` - Latest multimodal model
- `gpt-4o-mini` - Faster, lower cost
- `gpt-4-turbo` - Previous generation
- `gpt-3.5-turbo` - Fastest, lowest cost

For OpenAI-compatible APIs (LocalAI, LM Studio):

```bash
OPENAI_BASE_URL=http://localhost:1234/v1
```

#### Anthropic Claude Configuration

```bash
LLM_PROVIDER=anthropic
ANTHROPIC_API_KEY=your_api_key_here
ANTHROPIC_MODEL=claude-3-5-sonnet-20241022
```

Available models:

- `claude-3-5-sonnet-20241022` - Latest, best balance
- `claude-3-opus-20240229` - Highest capability
- `claude-3-sonnet-20240229` - Good balance

#### vLLM Configuration

```bash
LLM_PROVIDER=vllm
VLLM_BASE_URL=http://localhost:8000/v1
VLLM_MODEL=meta-llama/Llama-3.1-8B-Instruct
```

## Ollama Setup

If you choose Ollama, follow these additional steps:

### 1. Install Ollama

**macOS:**

```bash
curl -fsSL https://ollama.com/install.sh | sh
```

**Linux:**

```bash
curl -fsSL https://ollama.com/install.sh | sh
```

**Windows:**
Download from [ollama.com](https://ollama.com/download)

### 2. Start Ollama Server

```bash
ollama serve
```

The server will start on `http://localhost:11434`

### 3. Pull a Model

Choose and download a model:

```bash
# Recommended: Good balance
ollama pull llama3.1:8b

# Alternative: Faster but less capable
ollama pull mistral:7b

# Alternative: Better for code/technical
ollama pull codellama:13b
```

### 4. Verify Model

List installed models:

```bash
ollama list
```

### 5. Configure Agent

Update your `.env`:

```bash
LLM_PROVIDER=ollama
OLLAMA_MODEL=llama3.1:8b  # Use the model you pulled
```

## vLLM Setup

vLLM is a high-throughput and memory-efficient inference engine for LLMs. It's ideal for production deployments requiring fast inference.

### 1. Install vLLM

**With pip (requires Python 3.8+):**

```bash
pip install vllm
```

**With Docker:**

```bash
docker pull vllm/vllm-openai:latest
```

### 2. Start vLLM Server

**Using pip installation:**

```bash
# Example: Run Llama 3.1 8B Instruct
vllm serve meta-llama/Llama-3.1-8B-Instruct \
  --host 0.0.0.0 \
  --port 8000
```

**Using Docker:**

```bash
docker run --runtime nvidia --gpus all \
  -v ~/.cache/huggingface:/root/.cache/huggingface \
  -p 8000:8000 \
  --ipc=host \
  vllm/vllm-openai:latest \
  --model meta-llama/Llama-3.1-8B-Instruct
```

### 3. Verify Server

Test the server is running:

```bash
curl http://localhost:8000/v1/models
```

### 4. Configure Agent

Update your `.env`:

```bash
LLM_PROVIDER=vllm
VLLM_BASE_URL=http://localhost:8000/v1
VLLM_MODEL=meta-llama/Llama-3.1-8B-Instruct
```

## Running the Agent

Start the agent server:

```bash
npm start
# or
node agent-server.js
```

You should see output similar to:

```text
╔══════════════════════════════════════════════════════════╗
║       Quantum Hardware MCP Agent Server                  ║
╚══════════════════════════════════════════════════════════╝

Server starting at http://localhost:3021
LLM Provider: ollama
✓ ollama provider initialized successfully
Model: llama3.1:8b
✓ Successfully fetched 15 tools from quantum-hardware-mcp-server.

✓ Quantum Hardware MCP Agent Server ready at http://localhost:3021
```

## Using the REPL Chat Interface

For an interactive command-line chat experience:

```bash
npm run chat
# or
node chat.js
```

This will start an interactive chat session where you can:

- Ask questions about quantum hardware and devices
- Execute queries and commands through natural language
- View results in a formatted, easy-to-read manner
- Maintain conversation context across multiple queries
- **Answer agent queries** when the agent needs clarification (🆕 see [Query Capability](#query-capability))

**Chat Commands:**

- Type your quantum hardware-related questions naturally
- `/help` - Display available commands and examples
- `/clear` - Clear the chat history
- `/exit` or `/quit` - End the chat session

**Example Questions:**

- "List available quantum backends"
- "Get status of quantum device ibm_brisbane"
- "Show queue information for backend"
- "What quantum computers are available?"
- "Check calibration data for a device"
- "Get properties of quantum backend"

**Example Interactive Query:**

```text
You: Show me quantum device information
Agent: 🤔 Which quantum backend would you like information about?
      💡 Suggestions: 1. ibm_brisbane  2. ibm_kyoto  3. ibm_osaka
You: ibm_brisbane
Agent: 🤖 Here is the information for ibm_brisbane: ...
```

## API Endpoint

### `POST /chat`

This is the main endpoint for interacting with the agent.

**Request Body:**

```json
{
    "question": "Your question here",
    "history": []
}
```

- `question` (string, required): The user's question or prompt
- `history` (array, optional): Array of previous conversation turns in format:

  ```json
  [
    { "role": "user", "content": "Previous question" },
    { "role": "assistant", "content": "Previous answer" }
  ]
  ```

**Success Response (200 OK):**

```json
{
    "answer": "The formatted, natural language answer to your question."
}
```

**Error Response (500 Internal Server Error):**

```json
{
    "answer": "Sorry, there was an error processing your request."
}
```

## Query Capability

🆕 **New Feature**: The agent can now ask clarifying questions during execution!

The ReAct loop now supports **interactive queries**, allowing the agent to ask for clarification, request additional information, or confirm actions without exiting the loop.

### How It Works

When the agent needs more information, it can pause execution and ask you a question:

```text
You: "Get quantum device status"
Agent: 🤔 "Which quantum backend would you like to check?"
       💡 Suggestions: ibm_brisbane, ibm_kyoto, ibm_osaka
You: "ibm_brisbane"
Agent: 🤖 "Backend ibm_brisbane is online with 127 qubits..."
```

### Use Cases

- **Backend Selection**: When multiple backends match your criteria
- **Confirmation**: Before performing operations
- **Clarification**: When your request is ambiguous
- **Additional Info**: When more details are needed

### Example Scenarios

#### Scenario 1: Ambiguous Request

```text
You: "Show me device information"
Agent: "Which device: ibm_brisbane, ibm_kyoto, or ibm_osaka?"
You: "ibm_brisbane"
Agent: [shows device information]
```

#### Scenario 2: Confirmation

```text
You: "Submit a quantum job"
Agent: "This will use quantum resources. Confirm backend: ibm_brisbane?"
You: "Yes"
Agent: [submits job]
```

#### Scenario 3: Progressive Refinement

```text
You: "Check quantum backend"
Agent: "Which backend?"
You: "ibm_brisbane"
Agent: "What information: status, calibration, or properties?"
You: "status"
Agent: [shows status for ibm_brisbane]
```

### Query Capability Configuration

```bash
# Maximum iterations (includes queries)
MAX_ITERATIONS=10

# Session timeout (default: 5 minutes)
SESSION_TIMEOUT_MS=300000
```

## Troubleshooting

### Provider Package Not Found

If you see an error about missing provider package after a fresh install:

```text
Provider "ollama" is not available. Please install the required dependency:
npm install ollama
```

Solution: This shouldn't happen with a fresh `npm install` as all providers are installed as optional dependencies. If it does occur:

1. Verify node_modules exists: `ls -la node_modules/@google`
2. Reinstall dependencies: `rm -rf node_modules package-lock.json && npm install`
3. Check that you're running from the agent directory

### Ollama Connection Refused

```text
Cannot connect to Ollama at http://localhost:11434
```

Solution: Make sure Ollama is running:

```bash
ollama serve
```

### Ollama Model Not Found

```text
Model "llama3.1:8b" is not available
```

Solution: Pull the model:

```bash
ollama pull llama3.1:8b
```

### vLLM Connection Issues

```text
Cannot connect to vLLM server at http://localhost:8000/v1
```

Solution: Ensure vLLM is running and accessible:

```bash
# Check if vLLM is running
curl http://localhost:8000/v1/models

# Start vLLM if not running
vllm serve meta-llama/Llama-3.1-8B-Instruct --port 8000
```

### API Key Issues

For cloud providers (Gemini, OpenAI, Anthropic), ensure:

1. API key is correctly set in `.env`
2. API key is valid and has not expired
3. You have sufficient API credits/quota

### MCP Server Connection

```text
Failed to fetch tools from quantum-hardware-mcp-server
```

Solution: Ensure `quantum-hardware-mcp-server` is running and accessible at the configured URI. If authentication is enabled, ensure `MCP_API_KEY` is set correctly.

### MCP Authentication Error

```text
Unauthorized or authentication failed
```

Solution: Check that `MCP_API_KEY` in your `.env` matches the API key configured in the quantum-hardware-mcp server.

## Concurrency and Timeouts

The agent is designed to handle multiple concurrent `/chat` requests safely. To protect downstream services, requests are automatically queued and timed out based on these settings:

- `LLM_CONCURRENCY` (default: 4) - Maximum number of simultaneous LLM generation requests. If using Ollama locally, you may want to lower this to 1 to match Ollama's default single-model queueing.
- `MCP_CONCURRENCY` (default: 8) - Maximum number of simultaneous tool execution requests to the MCP server.
- `LLM_TIMEOUT_MS` (default: 60000) - How long to wait for the LLM to reply before returning a 504 error.
- `MCP_TIMEOUT_MS` (default: 30000) - How long to wait for an MCP tool execution to complete.

## Environment Variables Reference

| Variable | Required | Default | Description |
| -------- | -------- | ------- | ----------- |
| `LLM_PROVIDER` | Yes | `gemini` | LLM provider: gemini, ollama, openai, anthropic, or vllm |
| `QUANTUM_MCP_SERVER_URI` | Yes | - | URI of the Quantum Hardware MCP server |
| `MCP_API_KEY` | No | - | API key for quantum-hardware-mcp server authentication |
| `PORT` | No | `3021` | Port for the agent server |
| `QUANTUM_AGENT_URL` | No | `http://localhost:3021` | URL for chat client |

### Gemini Variables

| Variable | Required | Description |
| -------- | -------- | ----------- |
| `GEMINI_API_KEY` | Yes | Google Gemini API key |
| `GEMINI_MODEL` | Yes | Model name (e.g., gemini-2.5-pro) |

### Ollama Variables

| Variable | Required | Default | Description |
| -------- | -------- | ------- | ----------- |
| `OLLAMA_MODEL` | Yes | - | Model name (e.g., llama3.1:8b) |
| `OLLAMA_BASE_URL` | No | `http://localhost:11434` | Ollama server URL |
| `OLLAMA_TEMPERATURE` | No | `0.7` | Generation temperature (0.0-1.0) |
| `OLLAMA_KEEP_ALIVE` | No | `5m` | Keep model loaded duration |

### OpenAI Variables

| Variable | Required | Default | Description |
| -------- | -------- | ------- | ----------- |
| `OPENAI_API_KEY` | Yes | - | OpenAI API key |
| `OPENAI_MODEL` | Yes | - | Model name (e.g., gpt-4o) |
| `OPENAI_BASE_URL` | No | - | Custom API endpoint (for compatible APIs) |
| `OPENAI_TEMPERATURE` | No | `0.7` | Generation temperature (0.0-1.0) |

### Anthropic Variables

| Variable | Required | Default | Description |
| -------- | -------- | ------- | ----------- |
| `ANTHROPIC_API_KEY` | Yes | - | Anthropic API key |
| `ANTHROPIC_MODEL` | Yes | - | Model name (e.g., claude-3-5-sonnet-20241022) |
| `ANTHROPIC_TEMPERATURE` | No | `0.7` | Generation temperature (0.0-1.0) |
| `ANTHROPIC_MAX_TOKENS` | No | `4096` | Maximum tokens in response |

### vLLM Variables

| Variable | Required | Default | Description |
| -------- | -------- | ------- | ----------- |
| `VLLM_BASE_URL` | Yes | - | vLLM server URL (e.g., <http://localhost:8000/v1>) |
| `VLLM_MODEL` | Yes | - | Model name as loaded in vLLM |
| `VLLM_API_KEY` | No | `EMPTY` | API key (optional for local vLLM) |
| `VLLM_TEMPERATURE` | No | `0.7` | Generation temperature (0.0-1.0) |
| `VLLM_MAX_TOKENS` | No | `4096` | Maximum tokens in response |
| `VLLM_TOP_P` | No | `0.95` | Top-p sampling parameter |

## Switching Providers

To switch between providers:

1. Update `LLM_PROVIDER` in your `.env` file (all provider packages are already installed)
2. Configure the provider-specific variables
3. Restart the agent server

Example - switching from Gemini to Ollama:

```bash
# Update .env
LLM_PROVIDER=ollama
OLLAMA_MODEL=llama3.1:8b

# Restart
npm start
```

All provider packages are installed by default as optional dependencies, so you don't need to install them separately.

## Architecture

The agent uses a provider abstraction layer that allows it to work with any LLM:

```text
┌─────────────────┐
│   User Query    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Agent Server   │
└────────┬────────┘
         │
    ┌────┴─────────────────┐
    │                      │
    ▼                      ▼
┌─────────┐          ┌──────────┐
│   LLM   │          │Quantum   │
│Provider │          │Hardware  │
└─────────┘          │MCP Server│
    │                └──────────┘
    ├─ Gemini
    ├─ Ollama
    ├─ OpenAI
    ├─ Anthropic
    └─ vLLM
```

## License

ISC

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues.
