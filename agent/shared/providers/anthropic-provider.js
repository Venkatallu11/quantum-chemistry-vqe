/**
 * Anthropic Provider
 * Implements the BaseLLMProvider interface for Anthropic Claude
 */
const BaseLLMProvider = require('./base-provider');

class AnthropicProvider extends BaseLLMProvider {
    /**
     * Initialize the Anthropic provider
     * @param {Object} config - Configuration object
     * @param {string} config.apiKey - Anthropic API key
     * @param {string} config.model - Model name (e.g., claude-3-5-sonnet-20241022)
     * @param {number} config.temperature - Temperature for generation (default: 0.7)
     * @param {number} config.maxTokens - Maximum tokens in response (default: 4096)
     */
    async initialize(config) {
        if (!config.apiKey) {
            throw new Error('Anthropic API key is required');
        }
        if (!config.model) {
            throw new Error('Anthropic model name is required');
        }

        try {
            const Anthropic = require('@anthropic-ai/sdk');
            this.client = new Anthropic({
                apiKey: config.apiKey
            });
            this.model = config.model;
            this.temperature = config.temperature || 0.7;
            this.maxTokens = config.maxTokens || 4096;
        } catch (error) {
            if (error.code === 'MODULE_NOT_FOUND') {
                throw new Error(
                    'Failed to initialize Anthropic provider. ' +
                    'Please ensure @anthropic-ai/sdk is installed: ' +
                    'npm install @anthropic-ai/sdk'
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
            temperature: this.temperature,
            maxTokens: this.maxTokens
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

        // Call Anthropic messages API
        const response = await this.client.messages.create({
            model: this.model,
            max_tokens: this.maxTokens,
            temperature: this.temperature,
            messages: messages
        }, { signal });

        // Create new chat object with updated messages
        const assistantContent = response.content[0].text;
        const newChat = {
            ...chat,
            messages: [...messages, {
                role: 'assistant',
                content: assistantContent
            }]
        };

        return { response, chat: newChat };
    }

    /**
     * Extract text from Anthropic response
     * @param {Object} result - Anthropic response object
     * @returns {Promise<string>} Extracted text
     */
    async extractTextResponse(result) {
        return result.response.content[0].text;
    }

    /**
     * Transform standardized history to Anthropic format
     * @param {Array} history - Standardized history format
     * @returns {Array} Anthropic-specific history format
     */
    transformHistory(history) {
        if (!history || history.length === 0) {
            return [];
        }

        return history.map(entry => {
            // If already in Anthropic format (role + content string)
            if (entry.role && entry.content && typeof entry.content === 'string') {
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
            name: 'Anthropic Claude',
            provider: 'anthropic',
            model: this.model,
            version: '1.0.0'
        };
    }
}

module.exports = AnthropicProvider;
