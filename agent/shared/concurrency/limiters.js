// Use dynamic import since p-limit is an ESM module
let pLimit;

/**
 * Initializes the pLimit instance dynamically.
 */
async function initPLimit() {
    if (!pLimit) {
        const module = await import('p-limit');
        pLimit = module.default;
    }
}

// LLM Limiters mapped by provider name
const llmLimiters = {};

// Single MCP Limiter
let mcpLimiter = null;

/**
 * Gets or creates a concurrency limiter for the specified LLM provider.
 * @param {string} providerName - The LLM provider (e.g., 'gemini', 'ollama')
 * @returns {Promise<Function>} The limiter function
 */
async function getLLMLimiter(providerName) {
    await initPLimit();
    if (!llmLimiters[providerName]) {
        const concurrency = parseInt(process.env.LLM_CONCURRENCY) || 4;
        llmLimiters[providerName] = pLimit(concurrency);
    }
    return llmLimiters[providerName];
}

/**
 * Gets or creates the concurrency limiter for MCP server calls.
 * @returns {Promise<Function>} The limiter function
 */
async function getMCPLimiter() {
    await initPLimit();
    if (!mcpLimiter) {
        const concurrency = parseInt(process.env.MCP_CONCURRENCY) || 8;
        mcpLimiter = pLimit(concurrency);
    }
    return mcpLimiter;
}

module.exports = {
    getLLMLimiter,
    getMCPLimiter
};
