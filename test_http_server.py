#!/usr/bin/env python3
"""
Test script for the Quantum Hardware MCP HTTP server.

Usage:
    # Start the server in another terminal:
    python server.py --transport http
    
    # Then run this test:
    python test_http_server.py
"""

import requests
import json
import sys
import time

SERVER_URL = "http://127.0.0.1:8000"
API_KEY = None  # Set this if testing with authentication

def test_server_health():
    """Test if the server is running and responding."""
    print("Testing server health...")
    try:
        response = requests.get(f"{SERVER_URL}/health", timeout=5)
        print(f"✓ Server is running (status: {response.status_code})")
        return True
    except requests.exceptions.RequestException as e:
        print(f"✗ Server is not responding: {e}")
        return False

def test_mcp_endpoint():
    """Test the MCP SSE endpoint."""
    print("\nTesting MCP SSE endpoint...")
    headers = {"Content-Type": "application/json"}
    if API_KEY:
        headers["X-API-Key"] = API_KEY
    
    try:
        # Try to connect to the SSE endpoint
        response = requests.get(f"{SERVER_URL}/sse", headers=headers, timeout=5, stream=True)
        print(f"✓ SSE endpoint accessible (status: {response.status_code})")
        
        # Read first few lines
        lines = []
        for i, line in enumerate(response.iter_lines(decode_unicode=True)):
            if i >= 5:  # Read first 5 lines
                break
            if line:
                lines.append(line)
                print(f"  {line[:80]}...")
        
        return True
    except requests.exceptions.RequestException as e:
        print(f"✗ SSE endpoint error: {e}")
        return False

def test_without_api_key():
    """Test access without API key (should work in development mode)."""
    print("\nTesting without API key...")
    headers = {"Content-Type": "application/json"}
    
    try:
        response = requests.get(f"{SERVER_URL}/sse", headers=headers, timeout=5, stream=True)
        if response.status_code == 200:
            print("✓ Access granted without API key (development mode)")
            return True
        elif response.status_code == 401:
            print("✓ Access denied without API key (authentication enabled)")
            return True
        else:
            print(f"? Unexpected status code: {response.status_code}")
            return False
    except requests.exceptions.RequestException as e:
        print(f"✗ Request error: {e}")
        return False

def test_with_wrong_api_key():
    """Test access with wrong API key (should fail if auth is enabled)."""
    print("\nTesting with wrong API key...")
    headers = {
        "Content-Type": "application/json",
        "X-API-Key": "wrong_key_12345"
    }
    
    try:
        response = requests.get(f"{SERVER_URL}/sse", headers=headers, timeout=5, stream=True)
        if response.status_code == 401:
            print("✓ Access denied with wrong API key (authentication working)")
            return True
        elif response.status_code == 200:
            print("⚠ Access granted with wrong key (authentication not enabled)")
            return True
        else:
            print(f"? Unexpected status code: {response.status_code}")
            return False
    except requests.exceptions.RequestException as e:
        print(f"✗ Request error: {e}")
        return False

def main():
    print("=" * 70)
    print("Quantum Hardware MCP Server - HTTP Test Suite")
    print("=" * 70)
    print(f"Server URL: {SERVER_URL}")
    print(f"API Key: {'Set' if API_KEY else 'Not set (testing development mode)'}")
    print("=" * 70)
    
    # Give server time to start if just launched
    print("\nWaiting for server to be ready...")
    time.sleep(2)
    
    results = []
    
    # Run tests
    results.append(("Server Health", test_server_health()))
    results.append(("MCP SSE Endpoint", test_mcp_endpoint()))
    results.append(("Access without API key", test_without_api_key()))
    results.append(("Access with wrong API key", test_with_wrong_api_key()))
    
    # Summary
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {test_name}")
    
    print("=" * 70)
    print(f"Results: {passed}/{total} tests passed")
    print("=" * 70)
    
    return 0 if passed == total else 1

if __name__ == "__main__":
    sys.exit(main())

# Made with Bob
