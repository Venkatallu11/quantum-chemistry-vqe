#!/usr/bin/env node

const readline = require('readline');
const fs = require('fs').promises;
const path = require('path');
require('dotenv').config();

const AGENT_URL = process.env.QUANTUM_AGENT_URL || 'http://localhost:3021';
const chatHistory = [];

// Create readline interface for user input
const rl = readline.createInterface({
    input: process.stdin,
    output: process.stdout,
    prompt: '\n⚛️  Quantum Query> '
});

/**
 * Processes file references in the input text and replaces them with file contents
 * @param {string} input - User input that may contain @/path/to/file references
 * @returns {Promise<string>} Input with file contents inserted
 */
async function processFileReferences(input) {
    // Match @/path/to/file or @./relative/path patterns
    const fileRefPattern = /@(\.?\/[^\s]+)/g;
    const matches = [...input.matchAll(fileRefPattern)];
    
    if (matches.length === 0) {
        return input;
    }
    
    let processedInput = input;
    const fileContents = [];
    
    for (const match of matches) {
        const filePath = match[1];
        const fullMatch = match[0];
        
        try {
            // Resolve relative paths from current working directory
            const resolvedPath = path.resolve(process.cwd(), filePath);
            const content = await fs.readFile(resolvedPath, 'utf-8');
            
            // Format the file content with clear delimiters
            const formattedContent = `\n\n[File: ${filePath}]\n\`\`\`\n${content}\n\`\`\`\n`;
            fileContents.push({ match: fullMatch, content: formattedContent, path: filePath });
            
            console.log(`📎 Attached file: ${filePath} (${content.length} bytes)`);
        } catch (error) {
            const errorMsg = `\n[Error reading file ${filePath}: ${error.message}]\n`;
            fileContents.push({ match: fullMatch, content: errorMsg, path: filePath });
            console.log(`⚠️  Warning: Could not read file ${filePath}: ${error.message}`);
        }
    }
    
    // Replace all file references with their contents
    for (const { match, content } of fileContents) {
        processedInput = processedInput.replace(match, content);
    }
    
    return processedInput;
}

/**
 * Sends a chat message to the Quantum Hardware agent server
 * @param {string} question - The user's question
 * @param {Object} resumeContext - Optional context for resuming from a query
 * @returns {Promise<string>} The agent's answer
 */
async function chat(question, resumeContext = null) {
    try {
        const requestBody = {
            question,
            history: chatHistory
        };
        
        // Add resume context if continuing from a query
        if (resumeContext) {
            requestBody.resumeContext = resumeContext;
        }
        
        const response = await fetch(`${AGENT_URL}/chat`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify(requestBody),
        });

        if (!response.ok) {
            throw new Error(`Server responded with status: ${response.status}`);
        }

        const data = await response.json();
        
        // Handle query response
        if (data.status === 'query') {
            console.log('\n🤔 Agent Query:', data.query);
            
            // Show suggestions if available
            if (data.suggestions && data.suggestions.length > 0) {
                console.log('\n💡 Suggestions:');
                data.suggestions.forEach((s, i) => console.log(`  ${i + 1}. ${s}`));
            }
            
            // Get user's answer using a promise-based approach
            const userAnswer = await new Promise((resolve) => {
                rl.question('\n📝 Your answer> ', (answer) => {
                    resolve(answer.trim());
                });
            });
            
            // Resume with user's answer
            return await chat(userAnswer, {
                sessionId: data.context.sessionId,
                query: data.query
            });
        }
        
        // Handle final answer (both new format and backward compatible)
        const answer = data.answer || (data.status === 'complete' ? data.answer : 'No response');
        
        // Add to chat history in standardized format
        chatHistory.push({
            role: 'user',
            content: question
        });
        chatHistory.push({
            role: 'assistant',
            content: answer
        });

        return answer;
    } catch (error) {
        if (error.code === 'ECONNREFUSED') {
            return `Error: Cannot connect to agent server at ${AGENT_URL}. Make sure the agent server is running.`;
        }
        return `Error: ${error.message}`;
    }
}

/**
 * Standardize history entry to common format
 * @param {Object} entry - History entry
 * @returns {Object} Standardized entry
 */
function standardizeHistoryEntry(entry) {
    // Already standardized
    if (entry.role && entry.content && !entry.parts) {
        return entry;
    }
    
    // Gemini format with parts
    if (entry.role && entry.parts) {
        return {
            role: entry.role === 'model' ? 'assistant' : entry.role,
            content: entry.parts[0]?.text || ''
        };
    }
    
    return entry;
}

/**
 * Displays welcome message and instructions
 */
function displayWelcome() {
    console.log('\n╔══════════════════════════════════════════════════════════╗');
    console.log('║       Quantum Hardware Agent Chat Interface              ║');
    console.log('╚══════════════════════════════════════════════════════════╝');
    console.log('\nConnected to:', AGENT_URL);
    console.log('\nCommands:');
    console.log('  - Type your quantum hardware-related questions');
    console.log('  - Use @/path/to/file to insert file contents inline');
    console.log('  - Type "/exit" or "/quit" to end the session');
    console.log('  - Type "/clear" to clear chat history');
    console.log('  - Type "/help" to see this message again');
    console.log('\nExamples:');
    console.log('  - "List available quantum backends"');
    console.log('  - "Get status of quantum device ibm_brisbane"');
    console.log('  - "Show queue information for backend"');
    console.log('  - "What quantum computers are available?"');
    console.log('  - "Check calibration data for a device"');
    console.log('  - "Get properties of quantum backend"');
}

/**
 * Processes user commands
 * @param {string} input - User input
 * @returns {boolean} True if chat should continue, false to exit
 */
async function processInput(input) {
    const trimmed = input.trim().toLowerCase();

    switch (trimmed) {
        case '/exit':
        case '/quit':
            console.log('\n👋 Goodbye!\n');
            return false;

        case '/clear':
            chatHistory.length = 0;
            console.log('\n✓ Chat history cleared.\n');
            return true;

        case '/help':
            displayWelcome();
            return true;

        case '':
            return true;

        default:
            console.log('\n⏳ Processing...\n');
            // Process file references before sending to chat
            const processedInput = await processFileReferences(input);
            const answer = await chat(processedInput);
            console.log('🤖 Answer:', answer);
            return true;
    }
}

/**
 * Main chat loop
 */
async function main() {
    displayWelcome();
    
    rl.prompt();

    rl.on('line', async (line) => {
        const shouldContinue = await processInput(line);
        
        if (!shouldContinue) {
            rl.close();
            process.exit(0);
        }
        
        rl.prompt();
    });

    rl.on('close', () => {
        console.log('\n👋 Chat session ended.\n');
        process.exit(0);
    });

    // Handle Ctrl+C gracefully
    process.on('SIGINT', () => {
        console.log('\n\n👋 Chat session interrupted.\n');
        process.exit(0);
    });
}

// Start the chat
main().catch(error => {
    console.error('Fatal error:', error);
    process.exit(1);
});

// Made with Bob
