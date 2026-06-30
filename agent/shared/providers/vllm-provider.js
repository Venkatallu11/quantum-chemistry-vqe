/**
 * vLLM Provider
 * Implements the BaseLLMProvider interface for vLLM (OpenAI-compatible API)
 * vLLM is a high-throughput and memory-efficient inference engine for LLMs
 */
const BaseLLMProvider = require('./base-provider');

class VLLMProvider extends BaseLLMProvider {
    /**
     * Initialize the vLLM provider
     * @param {Object} config - Configuration object
     * @param {string} config.baseUrl - vLLM server URL (e.g., http://localhost:8000/v1)
     * @param {string} config.model - Model name (must match the model loaded in vLLM)
     * @param {string} config.apiKey - Optional API key (default: "EMPTY" for local vLLM)
     * @param {number} config.temperature - Temperature for generation (default: 0.7)
     * @param {number} config.maxTokens - Maximum tokens in response (default: 4096)
     * @param {number} config.topP - Top-p sampling parameter (default: 0.95)
     */
    async initialize(config) {
        if (!config.baseUrl) {
            throw new Error('vLLM base URL is required');
        }
        if (!config.model) {
            throw new Error('vLLM model name is required');
        }

        try {
            const { OpenAI } = require('openai');
            
            // vLLM uses OpenAI-compatible API
            this.client = new OpenAI({
                apiKey: config.apiKey || 'EMPTY', // vLLM doesn't require API key by default
                baseURL: config.baseUrl
            });
            
            this.model = config.model;
            this.baseUrl = config.baseUrl;
            this.temperature = config.temperature || 0.7;
            this.maxTokens = config.maxTokens || 4096;
            this.topP = config.topP || 0.95;
            
            // Verify connection to vLLM server
            await this.verifyConnection();
        } catch (error) {
            if (error.code === 'MODULE_NOT_FOUND') {
                throw new Error(
                    'Failed to initialize vLLM provider. ' +
                    'Please ensure openai package is installed: ' +
                    'npm install openai'
                );
            }
            throw error;
        }
    }

    /**
     * Verify connection to vLLM server
     * @returns {Promise<void>}
     * @throws {Error} If connection fails
     */
    async verifyConnection() {
        try {
            // Try to list models to verify connection
            await this.client.models.list();
        } catch (error) {
            throw new Error(
                `Cannot connect to vLLM server at ${this.baseUrl}. ` +
                `Please ensure vLLM is running and accessible. ` +
                `Error: ${error.message}`
            );
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
            temperature: this.temperature,
            maxTokens: this.maxTokens,
            topP: this.topP
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

        // Call vLLM chat completion API (OpenAI-compatible)
        const response = await this.client.chat.completions.create({
            model: this.model,
            messages: messages,
            temperature: this.temperature,
            max_tokens: this.maxTokens,
            top_p: this.topP
        }, { signal });

        // Create new chat object with updated messages
        const assistantMessage = response.choices[0].message;
        const newChat = {
            ...chat,
            messages: [...messages, assistantMessage]
        };

        return { response, chat: newChat };
    }

    /**
     * Extract text from vLLM response
     * @param {Object} result - vLLM response object
     * @returns {Promise<string>} Extracted text
     */
    async extractTextResponse(result) {
        // Handle both wrapped and unwrapped response formats
        const response = result.response || result;
        if (!response || !response.choices || !response.choices[0]) {
            throw new Error('Invalid vLLM response format: missing choices array');
        }
        return response.choices[0].message.content;
    }

    /**
     * Transform standardized history to vLLM format
     * @param {Array} history - Standardized history format
     * @returns {Array} vLLM-specific history format (OpenAI-compatible)
     */
    transformHistory(history) {
        if (!history || history.length === 0) {
            return [];
        }

        return history.map(entry => {
            // If already in OpenAI format (role + content)
            if (entry.role && entry.content && !entry.parts) {
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
            name: 'vLLM',
            provider: 'vllm',
            model: this.model,
            baseUrl: this.baseUrl,
            version: '1.0.0'
        };
    }
}

module.exports = VLLMProvider;

// Made with Bob
