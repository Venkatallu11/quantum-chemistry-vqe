/**
 * Ollama Provider
 * Implements the BaseLLMProvider interface for Ollama
 */
const BaseLLMProvider = require('./base-provider');

class OllamaProvider extends BaseLLMProvider {
    /**
     * Initialize the Ollama provider
     * @param {Object} config - Configuration object
     * @param {string} config.baseUrl - Ollama server URL (default: http://localhost:11434)
     * @param {string} config.model - Model name (e.g., llama3.1:8b)
     * @param {number} config.temperature - Temperature for generation (default: 0.7)
     * @param {string} config.keepAlive - Keep model loaded duration (default: 5m)
     */
    async initialize(config) {
        if (!config.model) {
            throw new Error('Ollama model name is required');
        }

        try {
            const { Ollama } = require('ollama');
            this.client = new Ollama({ 
                host: config.baseUrl || 'http://localhost:11434' 
            });
            this.model = config.model;
            this.temperature = config.temperature || 0.7;
            this.keepAlive = config.keepAlive || '5m';
            this.baseUrl = config.baseUrl || 'http://localhost:11434';

            // Verify the model is available
            await this.checkModelAvailability();
        } catch (error) {
            if (error.code === 'MODULE_NOT_FOUND') {
                throw new Error(
                    'Failed to initialize Ollama provider. ' +
                    'Please ensure ollama package is installed: ' +
                    'npm install ollama'
                );
            }
            throw error;
        }
    }

    /**
     * Check if the specified model is available
     * @throws {Error} If model is not available
     */
    async checkModelAvailability() {
        try {
            const models = await this.client.list();
            const modelExists = models.models.some(m => m.name === this.model);
            
            if (!modelExists) {
                console.warn(
                    `⚠️  Model "${this.model}" not found locally. ` +
                    `Available models: ${models.models.map(m => m.name).join(', ')}`
                );
                console.warn(
                    `To pull the model, run: ollama pull ${this.model}`
                );
                throw new Error(
                    `Model "${this.model}" is not available. ` +
                    `Please pull it first: ollama pull ${this.model}`
                );
            }
        } catch (error) {
            if (error.code === 'ECONNREFUSED') {
                throw new Error(
                    `Cannot connect to Ollama at ${this.baseUrl}. ` +
                    `Please ensure Ollama is running: ollama serve`
                );
            }
            throw error;
        }
    }

    /**
     * Create a new chat session with history
     * @param {Array} history - Conversation history in standardized format
     * @returns {Promise<Object>} Chat session object
     */
    async createChat(history) {
        const messages = this.transformHistory(history);
        return {
            messages,
            model: this.model,
            options: {
                temperature: this.temperature
            }
        };
    }

    /**
     * Send a message in the chat session
     * @param {Object} chat - Chat session object
     * @param {string} prompt - Message prompt
     * @returns {Promise<Object>} Response object
     */
    async sendMessage(chat, prompt) {
        const llmTimeoutMs = parseInt(process.env.LLM_TIMEOUT_MS) || 60000;
        const signal = AbortSignal.timeout(llmTimeoutMs);

        // Add the new user message to the conversation
        const messages = [...chat.messages, {
            role: 'user',
            content: prompt
        }];

        // Call Ollama chat API
        const response = await this.client.chat({
            model: this.model,
            messages: messages,
            stream: false,
            options: {
                temperature: this.temperature
            },
            keep_alive: this.keepAlive,
            signal
        });

        // Create new chat object with updated messages
        const newChat = {
            ...chat,
            messages: [...messages, {
                role: 'assistant',
                content: response.message.content
            }]
        };

        return { response, chat: newChat };
    }

    /**
     * Extract text from Ollama response
     * @param {Object} result - Ollama response object
     * @returns {Promise<string>} Extracted text
     */
    async extractTextResponse(result) {
        // `result` is already the raw Ollama response object (sendMessage returns
        // { response, chat } and the caller unwraps it before calling here).
        // Do NOT add an extra `.response` — that would be undefined.
        return result.message.content;
    }

    /**
     * Transform standardized history to Ollama format
     * @param {Array} history - Standardized history format
     * @returns {Array} Ollama-specific history format
     */
    transformHistory(history) {
        if (!history || history.length === 0) {
            return [];
        }

        return history.map(entry => {
            // If already in standard format with role and content
            if (entry.role && entry.content) {
                return {
                    role: entry.role === 'model' ? 'assistant' : entry.role,
                    content: entry.content
                };
            }

            // If in Gemini format with parts
            if (entry.role && entry.parts) {
                return {
                    role: entry.role === 'model' ? 'assistant' : entry.role,
                    content: entry.parts[0]?.text || ''
                };
            }

            // Fallback
            return entry;
        });
    }

    /**
     * Get provider metadata
     * @returns {Object} Provider information
     */
    getMetadata() {
        return {
            name: 'Ollama',
            provider: 'ollama',
            model: this.model,
            baseUrl: this.baseUrl,
            version: '1.0.0'
        };
    }

    /**
     * List available models from Ollama
     * @returns {Promise<Array>} List of available models
     */
    async listModels() {
        try {
            const models = await this.client.list();
            return models.models;
        } catch (error) {
            console.error('Failed to list Ollama models:', error);
            return [];
        }
    }
}

module.exports = OllamaProvider;
