const express = require('express');
const cors = require('cors');
const crypto = require('crypto');
const { Client } = require('@modelcontextprotocol/sdk/client/index.js');
const { SSEClientTransport } = require('@modelcontextprotocol/sdk/client/sse.js');
const ProviderFactory = require('./shared/providers/provider-factory');
const ProviderConfig = require('./shared/config/provider-config');
const { requestLoggerMiddleware, createLogger } = require('./lib/request-logger');
const { getLLMLimiter, getMCPLimiter } = require('./shared/concurrency/limiters');
require('dotenv').config();

// Qiskit specialist model config — configurable via .env
// Default: IBM's own Mistral model tuned for Qiskit (free on HuggingFace, runs via Ollama)
const QISKIT_CODE_MODEL = process.env.QISKIT_CODE_MODEL || 'hf.co/Qiskit/mistral-small-3.2-24b-qiskit-GGUF:latest';
const QISKIT_CODE_MODEL_URL = process.env.QISKIT_CODE_MODEL_URL || 'http://localhost:11434';
const QISKIT_CODE_MODEL_ENABLED = !!process.env.QISKIT_CODE_MODEL || !!process.env.QISKIT_CODE_MODEL_URL;

// --- Session Management for Query Support ---
/**
 * ReAct Session class to maintain state between query interactions
 */
class ReActSession {
    constructor(sessionId, question, history) {
        this.sessionId = sessionId;
        this.originalQuestion = question;
        this.history = history;
        this.iterations = 0;
        this.toolsUsed = [];
        this.queriesAsked = [];
        this.chat = null;
        this.currentPrompt = null;
        this.createdAt = Date.now();
    }
    
    isExpired(maxAgeMs = 300000) { // 5 minutes default
        return Date.now() - this.createdAt > maxAgeMs;
    }
    
    toJSON() {
        return {
            sessionId: this.sessionId,
            originalQuestion: this.originalQuestion,
            iterations: this.iterations,
            toolsUsed: this.toolsUsed,
            queriesAsked: this.queriesAsked,
            createdAt: this.createdAt
        };
    }
}

// In-memory session store (can be replaced with Redis for production)
const activeSessions = new Map();

/**
 * Generates a cryptographically secure session ID
 * @returns {string} Session ID
 */
const generateSessionId = () => {
    return crypto.randomBytes(16).toString('hex');
};

/**
 * Cleans up expired sessions
 */
const cleanupExpiredSessions = () => {
    const now = Date.now();
    let cleaned = 0;
    for (const [sessionId, session] of activeSessions.entries()) {
        if (session.isExpired()) {
            activeSessions.delete(sessionId);
            cleaned++;
        }
    }
    if (cleaned > 0) {
        console.log(`[Session Cleanup] Removed ${cleaned} expired session(s)`);
    }
};

// Run cleanup every minute
setInterval(cleanupExpiredSessions, 60000);

// Validate required environment variables and get provider config
let providerName, providerConfig;
try {
    const config = ProviderConfig.validate();
    providerName = config.provider;
    providerConfig = config.config;
} catch (error) {
    console.error('Configuration Error:', error.message);
    console.error('\nPlease check your .env file and ensure all required variables are set.');
    console.error('See .env.example for reference.');
    process.exit(1);
}

const app = express();
const port = process.env.PORT || 3021;

app.use(express.json());
app.use(cors());
app.use(requestLoggerMiddleware);

// --- LLM Provider Setup (abstracted) ---
let llmProvider;

// --- Quantum Hardware MCP Server setup ---
const QUANTUM_MCP_SERVER_URI = process.env.QUANTUM_MCP_SERVER_URI;
const MCP_API_KEY = process.env.MCP_API_KEY;

// Validate Quantum MCP server URI
if (!QUANTUM_MCP_SERVER_URI) {
    console.error('Missing required environment variable: QUANTUM_MCP_SERVER_URI');
    process.exit(1);
}

// MCP Client instance
let mcpClient = null;
let availableTools = [];
let formattedToolList = '';

// --- Skill Loading ---
let loadedSkills = '';

const loadSkills = async (logger = console) => {
    const skillsPath = process.env.QUANTUM_MCP_AGENT_SKILLS_PATH;
    if (!skillsPath) {
        return;
    }

    try {
        const fs = require('fs/promises');
        const path = require('path');
        
        // Check if directory exists
        try {
            const stats = await fs.stat(skillsPath);
            if (!stats.isDirectory()) {
                logger.warn(`[Skills] Path is not a directory: ${skillsPath}`);
                return;
            }
        } catch (err) {
            logger.warn(`[Skills] Skills path does not exist: ${skillsPath}`);
            return;
        }

        const entries = await fs.readdir(skillsPath, { withFileTypes: true });
        let skillsContent = [];

        for (const entry of entries) {
            if (entry.isDirectory()) {
                const skillMdPath = path.join(skillsPath, entry.name, 'SKILL.md');
                try {
                    const content = await fs.readFile(skillMdPath, 'utf8');
                    skillsContent.push(`--- SKILL: ${entry.name} ---\n${content}\n`);
                    logger.log(`[Skills] Loaded skill: ${entry.name}`);
                } catch (err) {
                    // SKILL.md doesn't exist or isn't readable, skip
                }
            }
        }

        if (skillsContent.length > 0) {
            loadedSkills = `\n\nAVAILABLE SKILLS:\n${skillsContent.join('\n')}\n`;
            logger.log(`[Skills] Successfully loaded ${skillsContent.length} skills.`);
        }
    } catch (error) {
        logger.error('[Skills] Error loading skills:', error.message);
    }
};

/**
 * Initialize MCP client connection
 * @returns {Promise<void>}
 */
const initializeMCPClient = async (logger = console) => {
    try {
        logger.log(`[MCP] Connecting to: ${QUANTUM_MCP_SERVER_URI}`);
        
        // Create SSE transport with authentication
        // Headers must be in requestInit, which the SDK merges with _commonHeaders()
        const transportOptions = {};
        
        if (MCP_API_KEY) {
            transportOptions.requestInit = {
                headers: {
                    'X-API-Key': MCP_API_KEY
                }
            };
            logger.log(`[MCP] Using X-API-Key authentication (key length: ${MCP_API_KEY.length})`);
        } else {
            logger.warn('[MCP] No MCP_API_KEY configured - server may reject connection if auth is required');
        }

        const transport = new SSEClientTransport(
            new URL(QUANTUM_MCP_SERVER_URI),
            transportOptions
        );

        // Create MCP client
        mcpClient = new Client(
            {
                name: 'quantum-hardware-mcp-agent',
                version: '1.0.0',
            },
            {
                capabilities: {},
            }
        );

        // Connect to server
        logger.log('[MCP] Attempting connection...');
        await mcpClient.connect(transport);
        logger.log('✓ Connected to quantum-hardware-mcp server via SSE');

        // Wait a moment for initialization to complete
        // The MCP protocol requires an initialization handshake
        await new Promise(resolve => setTimeout(resolve, 100));
        logger.log('[MCP] Initialization handshake complete');

        // Fetch available tools
        const toolsResponse = await mcpClient.listTools();
        availableTools = toolsResponse.tools.map(tool => ({
            name: tool.name,
            description: tool.description,
            inputSchema: tool.inputSchema,
        }));
        formattedToolList = JSON.stringify(availableTools, null, 2);
        logger.log(`✓ Successfully fetched ${availableTools.length} tools from quantum-hardware-mcp-server.`);

    } catch (error) {
        logger.error('FATAL: Failed to initialize MCP client:', error.message);
        logger.error('Error details:', error);
        if (error.message.includes('401')) {
            logger.error('\n⚠️  Authentication failed (401 Unauthorized)');
            logger.error('   Please check:');
            logger.error('   1. MCP_API_KEY is set in your .env file');
            logger.error('   2. The API key matches the quantum-hardware-mcp server configuration');
            logger.error('   3. Headers must be in requestInit.headers for SSEClientTransport');
        }
        throw error;
    }
};

/**
 * Calls an MCP tool and returns the result
 * @param {string} toolName - Name of the tool to call
 * @param {Object} toolArguments - Arguments for the tool
 * @param {Object} logger - Logger instance
 * @returns {Promise<Object>} Tool result data
 * @throws {Error} If tool call fails
 */
const callMCPTool = async (toolName, toolArguments, logger) => {
    try {
        const result = await mcpClient.callTool({
            name: toolName,
            arguments: toolArguments,
        });

        // Extract content from result
        if (result.content && Array.isArray(result.content)) {
            // If content is an array, extract text from each item
            const textContent = result.content
                .filter(item => item.type === 'text')
                .map(item => item.text)
                .join('\n');
            
            // Try to parse as JSON if possible
            try {
                return JSON.parse(textContent);
            } catch (e) {
                return { text: textContent };
            }
        }

        return result;
    } catch (error) {
        logger.error(`[MCP Tool] Error calling ${toolName}:`, error.message);
        throw error;
    }
};

// --- Qiskit Specialist Model ---
/**
 * Call the Qiskit-tuned local model for quantum code questions.
 * Runs via Ollama — no API key, no cost, fully local.
 * Model is configurable via QISKIT_CODE_MODEL in .env.
 */
const callQiskitCodeModel = async (prompt, logger) => {
    try {
        const { Ollama } = require('ollama');
        const client = new Ollama({ host: QISKIT_CODE_MODEL_URL });
        logger.log(`[Qiskit Model] Calling ${QISKIT_CODE_MODEL} at ${QISKIT_CODE_MODEL_URL}`);
        const response = await client.chat({
            model: QISKIT_CODE_MODEL,
            messages: [{ role: 'user', content: prompt }],
            stream: false,
        });
        return response.message.content;
    } catch (error) {
        if (error.code === 'ECONNREFUSED') {
            throw new Error(
                `Cannot connect to Ollama at ${QISKIT_CODE_MODEL_URL}. ` +
                `Please ensure Ollama is running: ollama serve`
            );
        }
        throw error;
    }
};

// --- ReAct Loop Implementation ---
/**
 * Runs the Reason + Act (ReAct) loop to autonomously handle multi-step tool execution
 * @param {Object} chat - LLM chat instance
 * @param {string} question - User's question
 * @param {string} formattedTools - Formatted tool list JSON string
 * @param {Object} req - Express request object
 * @param {Object} res - Express response object
 * @param {ReActSession} session - Optional session for resuming from query
 * @returns {Promise<void>}
 */
const runReActLoop = async (chat, question, formattedTools, req, res, session = null) => {
    let isComplete = false;
    let iterations = session ? session.iterations : 0;
    const maxIterations = parseInt(process.env.MAX_ITERATIONS) || 10;
    let finalAnswer = "I was unable to complete the task within the maximum number of steps.";
    let toolsUsed = session ? session.toolsUsed : [];
    
    let currentPrompt = session ? session.currentPrompt : `You are an expert quantum hardware autonomous agent.
Your goal is to answer the user's question by taking actions and analyzing their results.

Here is a list of available tools:
${formattedTools}

${loadedSkills ? `You have the following skills available to guide your actions:\n${loadedSkills}\n` : ''}
${QISKIT_CODE_MODEL_ENABLED ? `You have a Qiskit specialist model available (${QISKIT_CODE_MODEL}). This model is specifically trained on Qiskit and OpenQASM code and will produce better quantum code than you can. You MUST use it for any code-related request.` : ''}

User's request: "${question}"

Decide what to do next based on the request.
- If you need to use a tool to gather hardware data or run a circuit, return a JSON object with "action": "tool", "toolName": "<name>", and "toolArguments": {<args>}.
- If the user is asking about writing, debugging, or explaining ANY quantum code (OpenQASM, Qiskit Python, circuit design), you MUST use the Qiskit specialist model — return a JSON object with "action": "model_call", and "prompt": "<the exact question or code to send to the specialist>". Do NOT answer code questions yourself.
- If you need to ask the user a question or request clarification, return a JSON object with "action": "query", "query": "<your question to the user>", and optionally "suggestions": ["option1", "option2", ...].
- If you have gathered enough information and can answer the user's request (or if no tools are needed), return a JSON object with "action": "answer", and "finalAnswer": "<your natural language answer>".

Your response MUST be exactly ONE valid JSON object and nothing else. Do not use markdown blocks like \`\`\`json.`;

    while (!isComplete && iterations < maxIterations) {
        iterations++;
        req.logger.log(`[ReAct] Iteration ${iterations}...`);
        
        const llmLimiter = await getLLMLimiter(providerName);
        const sendResult = await llmLimiter(() => llmProvider.sendMessage(chat, currentPrompt));
        chat = sendResult.chat;
        const result = sendResult.response || sendResult;
        let textResponse = await llmProvider.extractTextResponse(result);
        
        let action;
        try {
            // Robust JSON extraction to handle nested code blocks in the output
            let jsonString = textResponse;
            try {
                action = JSON.parse(jsonString);
            } catch (e) {
                // If direct parse fails, try extracting between first { and last }
                const start = textResponse.indexOf('{');
                const end = textResponse.lastIndexOf('}');
                if (start !== -1 && end !== -1 && end > start) {
                    action = JSON.parse(textResponse.substring(start, end + 1));
                } else {
                    throw new Error("No JSON object found");
                }
            }
        } catch (error) {
            req.logger.warn('[ReAct] Failed to parse action JSON, falling back to treating response as final answer.');
            finalAnswer = textResponse;
            break;
        }

        if (action.action === 'answer' || (!action.toolName && !action.query && action.finalAnswer)) {
            finalAnswer = action.finalAnswer || action.answer || textResponse;
            isComplete = true;
            req.logger.log('[ReAct] Final answer reached.');
        } else if (action.action === 'query' || action.query) {
            // Handle query action - ask user for input
            req.logger.log(`[ReAct] Query to user: ${action.query}`);
            
            // Create or update session
            if (!session) {
                session = new ReActSession(generateSessionId(), question, req.body.history || []);
            }
            
            session.chat = chat;
            session.currentPrompt = currentPrompt;
            session.iterations = iterations;
            session.toolsUsed = toolsUsed;
            session.queriesAsked.push({
                query: action.query,
                iteration: iterations,
                timestamp: Date.now()
            });
            
            activeSessions.set(session.sessionId, session);
            req.logger.log(`[ReAct] Session ${session.sessionId} saved, awaiting user response`);
            
            // Return query response to client
            return res.json({
                status: 'query',
                query: action.query,
                context: {
                    sessionId: session.sessionId,
                    iteration: iterations,
                    toolsUsed: toolsUsed
                },
                suggestions: action.suggestions || []
            });
        } else if (action.action === 'model_call' || action.prompt && !action.toolName) {
            // Call the Qiskit specialist model for quantum code questions
            req.logger.log(`[ReAct] Qiskit model call: ${(action.prompt || '').substring(0, 80)}...`);
            if (!QISKIT_CODE_MODEL_ENABLED) {
                currentPrompt = `The Qiskit specialist model is not configured. Set QISKIT_CODE_MODEL in .env to enable it.\n\nFall back to answering from general knowledge or use a tool. Respond with ONLY valid JSON.`;
            } else {
                try {
                    const modelResponse = await callQiskitCodeModel(action.prompt, req.logger);
                    toolsUsed.push(`qiskit_model:${QISKIT_CODE_MODEL}`);
                    currentPrompt = `The Qiskit specialist model responded:\n${modelResponse}\n\nDecide what to do next:
- To use a tool, return JSON: { "action": "tool", "toolName": "<name>", "toolArguments": {<args>} }
- To call the Qiskit model again, return JSON: { "action": "model_call", "prompt": "<question>" }
- To ask the user a question, return JSON: { "action": "query", "query": "<question>", "suggestions": ["opt1", "opt2"] }
- To provide the final answer, return JSON: { "action": "answer", "finalAnswer": "<answer text>" }
Respond with ONLY valid JSON.`;
                } catch (error) {
                    req.logger.error(`[ReAct] Qiskit model call failed:`, error.message);
                    currentPrompt = `The Qiskit specialist model failed: ${error.message}\n\nFall back to answering from general knowledge or use a tool. Respond with ONLY valid JSON.`;
                }
            }
        } else if (action.action === 'tool' || action.toolName) {
            req.logger.log(`[ReAct] Tool selected: ${action.toolName}`);
            const toolArgs = action.toolArguments || action.arguments || {};
            try {
                const mcpLimiter = await getMCPLimiter();
                const resultData = await mcpLimiter(() => callMCPTool(action.toolName, toolArgs, req.logger));
                toolsUsed.push(action.toolName);
                currentPrompt = `Tool "${action.toolName}" returned:\n${JSON.stringify(resultData, null, 2)}\n\nDecide what to do next:
- To use another tool, return JSON: { "action": "tool", "toolName": "<name>", "toolArguments": {<args>} }
- To call the Qiskit specialist model for code questions, return JSON: { "action": "model_call", "prompt": "<question>" }
- To ask the user a question, return JSON: { "action": "query", "query": "<question>", "suggestions": ["opt1", "opt2"] }
- To provide the final answer, return JSON: { "action": "answer", "finalAnswer": "<answer text>" }
Respond with ONLY valid JSON.`;
            } catch (error) {
                req.logger.error(`[ReAct] Tool execution failed:`, error.message);
                currentPrompt = `Tool "${action.toolName}" failed with error: ${error.message}\n\nDecide what to do next (try another tool, ask the user, or provide an answer). Respond with ONLY valid JSON.`;
            }
        } else {
            req.logger.warn('[ReAct] Unrecognized action format, ending loop.');
            finalAnswer = textResponse;
            break;
        }
    }

    if (!isComplete) {
        req.logger.log('[ReAct] Max iterations reached.');
    }

    // Return final answer with metadata
    return res.json({
        status: 'complete',
        answer: finalAnswer,
        metadata: {
            iterations,
            toolsUsed,
            queriesAsked: session ? session.queriesAsked.length : 0
        }
    });
};

// --- Chat Endpoint ---
/**
 * Handles chat requests by selecting and calling appropriate tools from the Quantum Hardware MCP server
 * @param {Object} req - Express request object
 * @param {string} req.body.question - The user's question
 * @param {Array} req.body.history - Chat history for context
 * @param {Object} req.body.resumeContext - Optional context for resuming from a query
 * @param {Object} res - Express response object
 * @returns {Promise<void>} Sends JSON response with answer or error
 */
app.post('/chat', async (req, res) => {
    const { question, history, resumeContext } = req.body;

    if (availableTools.length === 0) {
        return res.status(500).json({
            status: 'error',
            answer: 'Tool list is not available. Please check the connection to the quantum-hardware-mcp-server.'
        });
    }

    try {
        // Check if resuming from a query
        if (resumeContext?.sessionId) {
            const session = activeSessions.get(resumeContext.sessionId);
            
            if (!session) {
                return res.status(400).json({
                    status: 'error',
                    error: 'Session expired or not found. Please start a new conversation.'
                });
            }
            
            if (session.isExpired()) {
                activeSessions.delete(resumeContext.sessionId);
                return res.status(400).json({
                    status: 'error',
                    error: 'Session expired. Please start a new conversation.'
                });
            }
            
            req.logger.log(`[Chat] Resuming session ${session.sessionId} with user answer: "${question}"`);
            
            // Update session with user's answer
            const lastQuery = session.queriesAsked[session.queriesAsked.length - 1];
            if (lastQuery) {
                lastQuery.answer = question;
                lastQuery.answeredAt = Date.now();
            }
            
            // Create prompt with user's answer
            session.currentPrompt = `User answered your query "${resumeContext.query || 'your question'}" with: "${question}"\n\nDecide what to do next:
- To use another tool, return JSON: { "action": "tool", "toolName": "<name>", "toolArguments": {<args>} }
- To ask another question, return JSON: { "action": "query", "query": "<question>", "suggestions": ["opt1", "opt2"] }
- To provide the final answer, return JSON: { "action": "answer", "finalAnswer": "<answer text>" }
Respond with ONLY valid JSON.`;
            
            // Resume ReAct loop with restored session
            await runReActLoop(session.chat, question, formattedToolList, req, res, session);
            
            // Clean up session only if we sent a final answer (not a query)
            if (res.headersSent) {
                const responseWasQuery = session.queriesAsked.length > resumeContext.iteration;
                
                if (!responseWasQuery) {
                    activeSessions.delete(resumeContext.sessionId);
                    req.logger.log(`[Chat] Session ${session.sessionId} completed and cleaned up`);
                } else {
                    req.logger.log(`[Chat] Session ${session.sessionId} waiting for next user response`);
                }
            }
            
            return;
        }
        
        // New conversation - start fresh ReAct loop
        req.logger.log(`[Chat] Starting new conversation: "${question}"`);
        
        // Standardize history format if needed
        const standardizedHistory = llmProvider.standardizeHistory(history || []);
        
        // Create chat with provider
        const chat = await llmProvider.createChat(standardizedHistory);

        // Run the ReAct loop to autonomously handle multi-step actions
        await runReActLoop(chat, question, formattedToolList, req, res);

    } catch (error) {
        req.logger.error('Error in LLM agent:', error);
        res.status(500).json({
            status: 'error',
            answer: 'Sorry, there was an error processing your request.'
        });
    }
});

// --- Server Startup ---
app.listen(port, async () => {
    const startupLogger = createLogger('startup');
    startupLogger.log('╔══════════════════════════════════════════════════════════╗');
    startupLogger.log('║       Quantum Hardware MCP Agent Server                  ║');
    startupLogger.log('╚══════════════════════════════════════════════════════════╝');
    startupLogger.log(`\nServer starting at http://localhost:${port}`);
    startupLogger.log(`LLM Provider: ${providerName}`);
    
    try {
        // Load Skills
        await loadSkills(startupLogger);

        // Initialize LLM provider
        llmProvider = await ProviderFactory.createProvider(providerName, providerConfig);
        const metadata = llmProvider.getMetadata();
        startupLogger.log(`Model: ${metadata.model || 'unknown'}`);
        
        // Initialize MCP client and fetch tools
        await initializeMCPClient(startupLogger);
        
        startupLogger.log(`\n✓ Quantum Hardware MCP Agent Server ready at http://localhost:${port}\n`);
    } catch (error) {
        startupLogger.error('\nFATAL: Failed to start server:', error.message);
        process.exit(1);
    }
});

// Graceful shutdown
process.on('SIGINT', async () => {
    console.log('\n\nShutting down gracefully...');
    if (mcpClient) {
        await mcpClient.close();
    }
    process.exit(0);
});

process.on('SIGTERM', async () => {
    console.log('\n\nShutting down gracefully...');
    if (mcpClient) {
        await mcpClient.close();
    }
    process.exit(0);
});

// Made with Bob
