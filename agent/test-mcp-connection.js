#!/usr/bin/env node

/**
 * Test script to verify MCP server connection and authentication
 */

const { Client } = require('@modelcontextprotocol/sdk/client/index.js');
const { SSEClientTransport } = require('@modelcontextprotocol/sdk/client/sse.js');
require('dotenv').config();

const QUANTUM_MCP_SERVER_URI = process.env.QUANTUM_MCP_SERVER_URI;
const MCP_API_KEY = process.env.MCP_API_KEY;

async function testConnection() {
    console.log('='.repeat(60));
    console.log('MCP Connection Test');
    console.log('='.repeat(60));
    console.log(`Server URI: ${QUANTUM_MCP_SERVER_URI}`);
    console.log(`API Key configured: ${MCP_API_KEY ? 'Yes' : 'No'}`);
    if (MCP_API_KEY) {
        console.log(`API Key length: ${MCP_API_KEY.length}`);
        console.log(`API Key preview: ${MCP_API_KEY.substring(0, 10)}...`);
    }
    console.log('='.repeat(60));

    try {
        // Test 1: Create transport with authentication
        console.log('\n[Test 1] Creating SSE transport...');
        
        // The SDK merges requestInit.headers into the fetch call
        const transportOptions = {};
        
        if (MCP_API_KEY) {
            // Headers must go in requestInit, which the SDK merges with _commonHeaders()
            transportOptions.requestInit = {
                headers: {
                    'X-API-Key': MCP_API_KEY
                }
            };
            console.log('  ✓ Added X-API-Key header to requestInit');
            console.log('  ✓ Headers:', JSON.stringify(transportOptions.requestInit.headers));
        }

        const transport = new SSEClientTransport(
            new URL(QUANTUM_MCP_SERVER_URI),
            transportOptions
        );
        console.log('  ✓ Transport created');

        // Test 2: Create client
        console.log('\n[Test 2] Creating MCP client...');
        const client = new Client(
            {
                name: 'quantum-hardware-mcp-test',
                version: '1.0.0',
            },
            {
                capabilities: {},
            }
        );
        console.log('  ✓ Client created');

        // Test 3: Connect
        console.log('\n[Test 3] Connecting to server...');
        await client.connect(transport);
        console.log('  ✓ Connected successfully!');

        // Test 4: List tools
        console.log('\n[Test 4] Fetching tools...');
        const toolsResponse = await client.listTools();
        console.log(`  ✓ Found ${toolsResponse.tools.length} tools:`);
        toolsResponse.tools.forEach((tool, i) => {
            console.log(`    ${i + 1}. ${tool.name}: ${tool.description}`);
        });

        // Test 5: Close connection
        console.log('\n[Test 5] Closing connection...');
        await client.close();
        console.log('  ✓ Connection closed');

        console.log('\n' + '='.repeat(60));
        console.log('✅ ALL TESTS PASSED');
        console.log('='.repeat(60));
        process.exit(0);

    } catch (error) {
        console.error('\n' + '='.repeat(60));
        console.error('❌ TEST FAILED');
        console.error('='.repeat(60));
        console.error('Error:', error.message);
        console.error('\nFull error:', error);
        
        if (error.message.includes('401')) {
            console.error('\n⚠️  Authentication Error (401 Unauthorized)');
            console.error('Possible causes:');
            console.error('  1. MCP_API_KEY is incorrect');
            console.error('  2. Server expects different authentication method');
            console.error('  3. API key format is wrong');
            console.error('\nTroubleshooting:');
            console.error('  - Verify MCP_API_KEY matches server configuration');
            console.error('  - Check server logs for authentication details');
            console.error('  - Try without MCP_API_KEY if server doesn\'t require auth');
        } else if (error.message.includes('404')) {
            console.error('\n⚠️  Endpoint Not Found (404)');
            console.error('Possible causes:');
            console.error('  1. Wrong endpoint path (should be /sse for FastMCP)');
            console.error('  2. Server not running on specified port');
            console.error('\nTroubleshooting:');
            console.error('  - Verify QUANTUM_MCP_SERVER_URI is correct');
            console.error('  - Check if server is running: curl http://127.0.0.1:8111/sse');
        } else if (error.code === 'ECONNREFUSED') {
            console.error('\n⚠️  Connection Refused');
            console.error('Possible causes:');
            console.error('  1. Server is not running');
            console.error('  2. Wrong port number');
            console.error('\nTroubleshooting:');
            console.error('  - Start the quantum-hardware-mcp server');
            console.error('  - Verify port 8111 is correct');
        }
        
        process.exit(1);
    }
}

// Run test
testConnection();

// Made with Bob
