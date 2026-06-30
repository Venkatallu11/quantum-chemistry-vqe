const crypto = require('crypto');

/**
 * Generates a unique request ID.
 * Uses crypto.randomUUID if available, falls back to a timestamp + random string.
 * @returns {string} Unique request ID
 */
function generateRequestId() {
    if (crypto.randomUUID) {
        return crypto.randomUUID();
    }
    return `${Date.now()}-${Math.random().toString(36).substring(2, 9)}`;
}

/**
 * Creates a logger bound to a specific request ID.
 * @param {string} reqId - The request ID
 * @returns {Object} Logger object with log, error, warn methods
 */
function createLogger(reqId) {
    const prefix = `[req=${reqId}]`;
    return {
        log: (...args) => console.log(prefix, ...args),
        error: (...args) => console.error(prefix, ...args),
        warn: (...args) => console.warn(prefix, ...args)
    };
}

/**
 * Express middleware to attach a request ID and logger to the request.
 * @param {Object} req - Express request
 * @param {Object} res - Express response
 * @param {Function} next - Next middleware function
 */
function requestLoggerMiddleware(req, res, next) {
    req.id = generateRequestId();
    req.logger = createLogger(req.id);
    next();
}

module.exports = {
    generateRequestId,
    createLogger,
    requestLoggerMiddleware
};

// Made with Bob
