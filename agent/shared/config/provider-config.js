/**
 * Provider Configuration Management
 * Handles validation and configuration for all LLM providers
 */
class ProviderConfig {
    static PROVIDER_REQUIREMENTS = {
        gemini: ['GEMINI_API_KEY', 'GEMINI_MODEL'],
        ollama: ['OLLAMA_MODEL'],  // BASE_URL optional, defaults to localhost
        openai: ['OPENAI_API_KEY', 'OPENAI_MODEL'],
        anthropic: ['ANTHROPIC_API_KEY', 'ANTHROPIC_MODEL'],
        vllm: ['VLLM_BASE_URL', 'VLLM_MODEL']  // API_KEY optional, defaults to "EMPTY"
    };

    /**
     * Validate environment configuration and return provider info
     * @returns {Object} Provider configuration
     * @throws {Error} If validation fails
     */
    static validate() {
        let provider = process.env.LLM_PROVIDER;
        
        // Backward compatibility: detect Gemini from old env vars
        if (!provider && process.env.GEMINI_API_KEY) {
            console.warn('⚠️  DEPRECATION: Please set LLM_PROVIDER=gemini explicitly in your .env file');
            provider = 'gemini';
        }
        
        // Default to Gemini for backward compatibility
        provider = provider || 'gemini';
        provider = provider.toLowerCase();
        
        const required = this.PROVIDER_REQUIREMENTS[provider];
        
        if (!required) {
            const available = Object.keys(this.PROVIDER_REQUIREMENTS).join(', ');
            throw new Error(
                `Unknown LLM provider: "${provider}". Available providers: ${available}`
            );
        }

        // Validate provider-specific requirements
        const missing = required.filter(v => !process.env[v]);
        
        if (missing.length > 0) {
            throw new Error(
                `Missing required environment variables for ${provider} provider: ${missing.join(', ')}\n` +
                `Please check your .env file and ensure all required variables are set.`
            );
        }

        // Note: Each agent validates its own server URI requirements
        // (IBM_MCP_SERVER_URI, MONGODB_MCP_SERVER_URI, PHP_MCP_SERVER_URI, etc.)

        return {
            provider,
            config: this.getProviderConfig(provider)
        };
    }

    /**
     * Get provider-specific configuration from environment variables
     * @param {string} provider - Provider name
     * @returns {Object} Provider configuration
     */
    static getProviderConfig(provider) {
        const configs = {
            gemini: {
                apiKey: process.env.GEMINI_API_KEY,
                model: process.env.GEMINI_MODEL
            },
            ollama: {
                baseUrl: process.env.OLLAMA_BASE_URL || 'http://localhost:11434',
                model: process.env.OLLAMA_MODEL,
                temperature: parseFloat(process.env.OLLAMA_TEMPERATURE || '0.7'),
                keepAlive: process.env.OLLAMA_KEEP_ALIVE || '5m'
            },
            openai: {
                apiKey: process.env.OPENAI_API_KEY,
                model: process.env.OPENAI_MODEL,
                baseUrl: process.env.OPENAI_BASE_URL,  // Optional for compatible APIs
                temperature: parseFloat(process.env.OPENAI_TEMPERATURE || '0.7')
            },
            anthropic: {
                apiKey: process.env.ANTHROPIC_API_KEY,
                model: process.env.ANTHROPIC_MODEL,
                temperature: parseFloat(process.env.ANTHROPIC_TEMPERATURE || '0.7'),
                maxTokens: parseInt(process.env.ANTHROPIC_MAX_TOKENS || '4096')
            },
            vllm: {
                baseUrl: process.env.VLLM_BASE_URL,
                model: process.env.VLLM_MODEL,
                apiKey: process.env.VLLM_API_KEY || 'EMPTY',  // vLLM doesn't require API key by default
                temperature: parseFloat(process.env.VLLM_TEMPERATURE || '0.7'),
                maxTokens: parseInt(process.env.VLLM_MAX_TOKENS || '4096'),
                topP: parseFloat(process.env.VLLM_TOP_P || '0.95')
            }
        };

        return configs[provider];
    }

    /**
     * Get list of available providers
     * @returns {Array<string>} List of provider names
     */
    static getAvailableProviders() {
        return Object.keys(this.PROVIDER_REQUIREMENTS);
    }

    /**
     * Check if a specific provider is available (dependencies installed)
     * @param {string} provider - Provider name
     * @returns {boolean} True if provider dependencies are available
     */
    static isProviderAvailable(provider) {
        try {
            switch (provider) {
                case 'gemini':
                    require.resolve('@google/generative-ai');
                    return true;
                case 'ollama':
                    require.resolve('ollama');
                    return true;
                case 'openai':
                    require.resolve('openai');
                    return true;
                case 'anthropic':
                    require.resolve('@anthropic-ai/sdk');
                    return true;
                case 'vllm':
                    require.resolve('openai');  // vLLM uses openai package
                    return true;
                default:
                    return false;
            }
        } catch (error) {
            return false;
        }
    }
}

module.exports = ProviderConfig;
