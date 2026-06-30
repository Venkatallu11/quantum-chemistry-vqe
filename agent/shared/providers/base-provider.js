/**
 * Base Provider Class
 * Abstract base class that all LLM providers must extend
 */
class BaseLLMProvider {
    /**
     * Initialize the provider with configuration
     * @param {Object} config - Provider-specific configuration
     * @returns {Promise<void>}
     */
    async initialize(config) {
        throw new Error('initialize() must be implemented by provider subclass');
    }

    /**
     * Create a new chat session with history
     * @param {Array} history - Conversation history in standardized format
     * @returns {Promise<Object>} Chat session object
     */
    async createChat(history) {
        throw new Error('createChat() must be implemented by provider subclass');
    }

    /**
     * Send a message in the chat session
     * @param {Object} chat - Chat session object
     * @param {string} prompt - Message prompt
     * @returns {Promise<{response: Object, chat: Object}>} Response and new chat object
     */
    async sendMessage(chat, prompt) {
        throw new Error('sendMessage() must be implemented by provider subclass');
    }

    /**
     * Extract text from provider response
     * @param {Object} result - Provider response object
     * @returns {Promise<string>} Extracted text
     */
    async extractTextResponse(result) {
        throw new Error('extractTextResponse() must be implemented by provider subclass');
    }

    /**
     * Transform standardized history to provider format
     * @param {Array} history - Standardized history format
     * @returns {Array} Provider-specific history format
     */
    transformHistory(history) {
        throw new Error('transformHistory() must be implemented by provider subclass');
    }

    /**
     * Get provider metadata
     * @returns {Object} Provider information
     */
    getMetadata() {
        return {
            name: this.constructor.name,
            version: '1.0.0'
        };
    }

    /**
     * Standardize history format from various input formats
     * @param {Array} history - History in any format
     * @returns {Array} Standardized history format
     */
    standardizeHistory(history) {
        if (!history || history.length === 0) {
            return [];
        }

        return history.map(entry => {
            // If already in standardized format
            if (entry.role && entry.content) {
                return entry;
            }

            // If in Gemini format with parts
            if (entry.role && entry.parts) {
                return {
                    role: entry.role === 'model' ? 'assistant' : entry.role,
                    content: entry.parts[0]?.text || ''
                };
            }

            // Default: return as-is
            return entry;
        });
    }
}

module.exports = BaseLLMProvider;
