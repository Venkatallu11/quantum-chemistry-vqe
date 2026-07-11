#!/usr/bin/env python3
"""
azure_backend.py — connect to the real Azure Quantum workspace and compare
available targets before submitting anything.
============================================================================
Mirrors the IBM side's list_devices/compare_devices pattern: never hardcode
a target name and submit blind. Always connect, list what's ACTUALLY
available right now (provider, target id, availability, queue time), and
let the caller pick from that real list.

Auth: DeviceCodeCredential — prints a URL + one-time code to the console,
you approve it in any browser, and azure-identity caches the token locally
(~/.IdentityService) so you don't have to re-auth every run.

Run standalone to verify connectivity and print the live target list:
    python vqe/azure_backend.py
"""
import os
from dotenv import load_dotenv
from azure.identity import DeviceCodeCredential, TokenCachePersistenceOptions
from azure.quantum import Workspace

load_dotenv()


def connect_workspace() -> Workspace:
    """Connect to the Azure Quantum workspace named in .env. Raises a clear
    error if the required env vars are missing instead of silently
    defaulting to anything."""
    resource_id = os.environ.get("AZURE_QUANTUM_RESOURCE_ID")
    location = os.environ.get("AZURE_QUANTUM_LOCATION")
    tenant_id = os.environ.get("AZURE_QUANTUM_TENANT_ID")
    if not resource_id or not location:
        raise RuntimeError(
            "AZURE_QUANTUM_RESOURCE_ID and AZURE_QUANTUM_LOCATION must be set "
            "in .env (see .env.example)."
        )

    credential = DeviceCodeCredential(
        tenant_id=tenant_id,
        cache_persistence_options=TokenCachePersistenceOptions(name="azure-quantum-vqe"),
    )
    return Workspace(resource_id=resource_id, location=location, credential=credential)


def compare_targets(workspace: Workspace, provider: str | None = None) -> list[dict]:
    """List every target actually available in this workspace right now,
    with live availability and queue time. Optionally filter by provider
    id (e.g. "quantinuum"). Returns the raw comparison data; does NOT pick
    one for you.
    """
    rows = []
    for t in workspace.get_targets():
        if provider and t.provider_id.lower() != provider.lower():
            continue
        try:
            queue_time = t.average_queue_time()
        except Exception:
            queue_time = None
        rows.append({
            "provider_id": t.provider_id,
            "target_id": t.name,
            "current_availability": t.current_availability,
            "average_queue_time_s": queue_time,
        })
    return rows


def print_comparison(rows: list[dict]) -> None:
    if not rows:
        print("  No targets found (check provider access on this workspace).")
        return
    print(f"\n  {'provider':<15} {'target':<30} {'availability':<15} {'avg queue (s)':<15}")
    print("  " + "-" * 78)
    for r in rows:
        q = "n/a" if r["average_queue_time_s"] is None else str(r["average_queue_time_s"])
        print(f"  {r['provider_id']:<15} {r['target_id']:<30} "
              f"{str(r['current_availability']):<15} {q:<15}")


def main():
    print("\n" + "=" * 70)
    print("  Azure Quantum — workspace connection + live target comparison")
    print("=" * 70)
    print("\n  Connecting (approve the device-code login in your browser if prompted)...")
    ws = connect_workspace()
    print(f"  Connected to workspace: {ws.name} ({ws.location})")

    rows = compare_targets(ws)
    print_comparison(rows)
    print(f"\n  {len(rows)} target(s) total.\n")
    return rows


if __name__ == "__main__":
    main()
