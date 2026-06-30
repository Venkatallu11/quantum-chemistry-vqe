# Shared LLM Provider Infrastructure

This directory contains the shared LLM provider abstraction layer used by the quantum-hardware-mcp-agent.

## Overview

The shared provider infrastructure enables the agent to support multiple LLM providers (Gemini, Ollama, OpenAI, Anthropic, vLLM) through a unified interface.

## Directory Structure

```
shared/
├── providers/                  # Provider implementations
│   ├── base-provider.js        # Abstract base class
│   ├── gemini-provider.js      # Google Gemini adapter
│   ├── ollama-provider.js      # Ollama (local) adapter
│   ├── openai-provider.js      # OpenAI adapter
│   ├── anthropic-provider.js   # Anthropic Claude adapter
│   └── provider-factory.js     # Factory for creating providers
└── config/                     # Configuration management
    └── provider-config.js      # Configuration validation
```

## Supported Providers

| Provider | Type | Cost | Privacy | Setup Difficulty |
|----------|------|------|---------|------------------|
| **Gemini** | Cloud API | Paid | Cloud | Easy |
| **Ollama** | Local | Free | Local | Medium |
| **OpenAI** | Cloud API | Paid | Cloud | Easy |
| **Anthropic** | Cloud API | Paid | Cloud | Easy |

## Usage in Agent

The quantum-hardware-mcp-agent uses the providers like this:

```javascript
const ProviderFactory = require('./shared/providers/provider-factory');
const ProviderConfig = require('./shared/config/provider-config');

// Validate and get configuration
const { provider, config } = ProviderConfig.validate();

// Create provider instance
const llmProvider = await ProviderFactory.createProvider(provider, config);

// Use provider
const chat = await llmProvider.createChat(history);
const result = await llmProvider.sendMessage(chat, prompt);
const text = await llmProvider.extractTextResponse(result);
```

## Provider Interface

All providers implement these methods:

```javascript
class BaseLLMProvider {
    async initialize(config)           // Setup provider with config
    async createChat(history)          // Create chat session
    async sendMessage(chat, prompt)    // Send message in session
    async extractTextResponse(result)  // Extract text from response
    transformHistory(history)          // Convert history format
    getMetadata()                      // Get provider info
}
```

## Configuration

Each agent configures the provider through environment variables:

```bash
# Choose provider
LLM_PROVIDER=ollama

# Configure chosen provider
OLLAMA_MODEL=llama3.1:8b
OLLAMA_BASE_URL=http://localhost:11434
```

See each agent's `.env.example` for complete configuration options.

## Adding a New Provider

To add a new LLM provider:

1. **Create provider class** in `providers/`:

   ```javascript
   // providers/new-provider.js
   const BaseLLMProvider = require('./base-provider');
   
   class NewProvider extends BaseLLMProvider {
       async initialize(config) { /* ... */ }
       async createChat(history) { /* ... */ }
       async sendMessage(chat, prompt) { /* ... */ }
       async extractTextResponse(result) { /* ... */ }
       transformHistory(history) { /* ... */ }
   }
   
   module.exports = NewProvider;
   ```

2. **Register in factory** (`provider-factory.js`):

   ```javascript
   static PROVIDERS = {
       gemini: './gemini-provider',
       ollama: './ollama-provider',
       openai: './openai-provider',
       anthropic: './anthropic-provider',
       newprovider: './new-provider'  // Add this
   };
   ```

3. **Add configuration** (`config/provider-config.js`):

   ```javascript
   static PROVIDER_REQUIREMENTS = {
       // ... existing
       newprovider: ['NEWPROVIDER_API_KEY', 'NEWPROVIDER_MODEL']
   };
   
   static getProviderConfig(provider) {
       const configs = {
           // ... existing
           newprovider: {
               apiKey: process.env.NEWPROVIDER_API_KEY,
               model: process.env.NEWPROVIDER_MODEL
           }
       };
       return configs[provider];
   }
   ```

4. **Update package.json** (root and each agent):

   ```json
   "optionalDependencies": {
       "newprovider-sdk": "^1.0.0"
   }
   ```

5. **Document** in each agent's `.env.example` and README

## Benefits of Shared Infrastructure

✅ **Single Source of Truth** - One implementation, multiple consumers
✅ **Consistency** - All agents behave identically
✅ **Easy Maintenance** - Fix bugs once, all agents benefit
✅ **Easy Extension** - Add providers once, all agents can use them
✅ **Reduced Code** - ~1,200 lines shared vs ~4,800 duplicated
✅ **Easier Testing** - Test providers once with confidence

## Deployment

The `shared/` directory is included within the agent directory:

```bash
# Directory structure
agent/
├── shared/              # Shared infrastructure
│   ├── providers/
│   ├── config/
│   └── concurrency/
├── agent-server.js      # References ./shared/
├── chat.js
└── ...
```

## Version Compatibility

All providers in this directory follow semantic versioning:

- **Major version** - Breaking changes to provider interface
- **Minor version** - New providers or features (backward compatible)
- **Patch version** - Bug fixes

Current version: **1.0.0**

## Testing

Test the shared providers independently:

```bash
# Run provider tests (if available)
cd shared
npm test

# Test with a specific agent
cd ../ibmi-mcp-agent
npm start
```

## Troubleshooting

### "Provider not available" Error

```
Provider "ollama" is not available. Please install the required dependency:
npm install ollama
```

**Solution:** Install from root:

```bash
npm install ollama
```

### "MODULE_NOT_FOUND" for Provider

Ensure the shared directory is accessible:

```bash
ls -la shared/providers/
```

### Import Path Issues

The agent uses relative paths from the agent directory:

```javascript
require('./shared/providers/provider-factory')  // Correct
require('../shared/providers/provider-factory')  // Wrong (old multi-agent structure)
```

## License

ISC

## Contributing

When modifying shared providers:

1. Test with ALL agents
2. Maintain backward compatibility
3. Update version appropriately
4. Document changes in all agent READMEs
