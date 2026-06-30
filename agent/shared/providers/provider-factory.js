/**
 * Provider Factory
 * Creates and initializes LLM provider instances
 */
const ProviderConfig = require('../config/provider-config');

class ProviderFactory {
    static PROVIDERS = {
        gemini: './gemini-provider',
        ollama: './ollama-provider',
        openai: './openai-provider',
        anthropic: './anthropic-provider',
        vllm: './vllm-provider'
    };

    /**
     * Create and initialize a provider
     * @param {string} providerName - Name of the provider to create
     * @param {Object} config - Provider configuration
     * @returns {Promise<BaseLLMProvider>} Initialized provider instance
     * @throws {Error} If provider not found or initialization fails
     */
    static async createProvider(providerName, config) {
        const providerPath = this.PROVIDERS[providerName.toLowerCase()];
        
        if (!providerPath) {
            const available = Object.keys(this.PROVIDERS).join(', ');
            throw new Error(
                `Unknown provider: "${providerName}". Available providers: ${available}`
            );
        }

        // Check if provider dependencies are installed
        if (!ProviderConfig.isProviderAvailable(providerName)) {
            throw new Error(
                `Provider "${providerName}" is not available. ` +
                `Please install the required dependency:\n` +
                this.getInstallCommand(providerName)
            );
        }

        try {
            const ProviderClass = require(providerPath);
            const provider = new ProviderClass();
            await provider.initialize(config);
            
            console.log(`✓ ${providerName} provider initialized successfully`);
            return provider;
        } catch (error) {
            if (error.code === 'MODULE_NOT_FOUND') {
                throw new Error(
                    `Failed to load ${providerName} provider. ` +
                    `Please install the required dependency:\n` +
                    this.getInstallCommand(providerName)
                );
            }
            throw new Error(`Failed to initialize ${providerName} provider: ${error.message}`);
        }
    }

    /**
     * Get npm install command for a provider
     * @param {string} providerName - Provider name
     * @returns {string} npm install command
     */
    static getInstallCommand(providerName) {
        const commands = {
            gemini: 'npm install @google/generative-ai',
            ollama: 'npm install ollama',
            openai: 'npm install openai',
            anthropic: 'npm install @anthropic-ai/sdk',
            vllm: 'npm install openai'
        };
        return commands[providerName] || 'npm install';
    }

    /**
     * Get list of available providers
     * @returns {Array<string>} List of provider names
     */
    static getAvailableProviders() {
        return Object.keys(this.PROVIDERS);
    }

    /**
     * Get list of installed providers
     * @returns {Array<string>} List of installed provider names
     */
    static getInstalledProviders() {
        return this.getAvailableProviders().filter(provider => 
            ProviderConfig.isProviderAvailable(provider)
        );
    }
}

module.exports = ProviderFactory;
