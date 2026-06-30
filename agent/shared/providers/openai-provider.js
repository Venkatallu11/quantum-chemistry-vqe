/**
 * OpenAI Provider
 * Implements the BaseLLMProvider interface for OpenAI (and compatible APIs)
 */
const BaseLLMProvider = require('./base-provider');

class OpenAIProvider extends BaseLLMProvider {
    /**
     * Initialize the OpenAI provider
     * @param {Object} config - Configuration object
     * @param {string} config.apiKey - OpenAI API key
     * @param {string} config.model - Model name (e.g., gpt-4o, gpt-3.5-turbo)
     * @param {string} config.baseUrl - Optional base URL for compatible APIs
     * @param {number} config.temperature - Temperature for generation (default: 0.7)
     */
    async initialize(config) {
        if (!config.apiKey) {
            throw new Error('OpenAI API key is required');
        }
        if (!config.model) {
            throw new Error('OpenAI model name is required');
        }

        try {
            const { OpenAI } = require('openai');
            
            const clientConfig = {
                apiKey: config.apiKey
            };
            
            // Support for OpenAI-compatible APIs (LocalAI, LM Studio, etc.)
            if (config.baseUrl) {
                clientConfig.baseURL = config.baseUrl;
            }
            
            this.client = new OpenAI(clientConfig);
            this.model = config.model;
            this.temperature = config.temperature || 0.7;
            this.baseUrl = config.baseUrl;
        } catch (error) {
            if (error.code === 'MODULE_NOT_FOUND') {
                throw new Error(
                    'Failed to initialize OpenAI provider. ' +
                    'Please ensure openai package is installed: ' +
                    'npm install openai'
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
            temperature: this.temperature
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

        // Call OpenAI chat completion API
        const response = await this.client.chat.completions.create({
            model: this.model,
            messages: messages,
            temperature: this.temperature
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
     * Extract text from OpenAI response
     * @param {Object} result - OpenAI response object
     * @returns {Promise<string>} Extracted text
     */
    async extractTextResponse(result) {
        return result.response.choices[0].message.content;
    }

    /**
     * Transform standardized history to OpenAI format
     * @param {Array} history - Standardized history format
     * @returns {Array} OpenAI-specific history format
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
            name: 'OpenAI',
            provider: 'openai',
            model: this.model,
            baseUrl: this.baseUrl || 'https://api.openai.com/v1',
            version: '1.0.0'
        };
    }
}

module.exports = OpenAIProvider;
