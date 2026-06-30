# Project: Quantum Hardware MCP Server

## What this is
An open-source MCP server (Python) that gives AI assistants live data
about IBM Quantum computers. Built by a new grad learning in public.

## Goals (in order)
1. Working MCP server: list_devices, get_device_details, compare_devices, queue_status
2. Connects to Claude Desktop
3. Published on GitHub (MIT license) + MCP registries
4. v2 later: autonomous daily-report agent on top

## How to work (efficiency rules)
- Before starting ANY task, re-read this file and confirm in one line: "Checked CLAUDE.md, following the rules."
- ONE feature per session. Finish, test, stop. Never start extra work I didn't ask for.
- Keep answers short. No long explanations unless I ask. Code > talk.
- Before any task: state a 2-3 step plan in one line each. Wait for my OK if it touches more than 2 files.
- Don't re-read files you haven't changed. Don't re-run tests that already passed.
- Smallest possible change that works. No refactors, no "improvements" I didn't request.
- Minimal dependencies: only the MCP SDK and qiskit-ibm-runtime. Ask before adding anything else.
- If stuck after 3 attempts at the same error, STOP and explain the problem in plain English instead of looping.

## Learning rules (I'm new)
- Comment code heavily, prefer simple over clever.
- When I ask "explain", use plain English and short analogies, no jargon.
- After finishing a feature, give me: 1-line summary + one interview question about it.

## Safety rules
- NEVER put the IBM API token in code. Use .env, keep .env in .gitignore. Check .gitignore before every commit.
- Never run destructive commands (rm -rf, force push) without asking.
- git commit after each working feature with a clear message. Never commit broken code.

## Definition of done (per feature)
Code works + tested + committed + README updated + 1-line summary to me.
