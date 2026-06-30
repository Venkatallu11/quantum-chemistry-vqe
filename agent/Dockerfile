FROM node:22-slim

WORKDIR /app

# Install dependencies first — separate layer for cache efficiency
COPY package.json package-lock.json ./
RUN npm ci --omit=dev

# Copy agent code
COPY . .

# All secrets (ANTHROPIC_API_KEY etc.) must be passed at runtime via -e or --env-file
# QUANTUM_MCP_SERVER_URI is set by docker-compose to point at the mcp-server container
EXPOSE 3021

CMD ["node", "agent-server.js"]
