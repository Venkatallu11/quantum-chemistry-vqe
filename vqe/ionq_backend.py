#!/usr/bin/env python3
"""
ionq_backend.py — connect to the real IonQ cloud API and list what's
actually available before submitting anything.
============================================================================
Mirrors azure_backend.py's compare_targets pattern (which itself mirrors
the IBM side's list_devices/compare_devices pattern): never hardcode a
target and submit blind. Always connect, list what's ACTUALLY live right
now, and let the caller pick from that real list.

Auth: reads IONQ_API_KEY from .env (never hardcode a token in code).

Run standalone to verify connectivity and print the live backend list:
    python vqe/ionq_backend.py
"""
import os
from dotenv import load_dotenv
from qiskit_ionq import IonQProvider

load_dotenv()


def connect_provider() -> IonQProvider:
    """Connect to the IonQ cloud API using the token in .env. Raises a
    clear error if the required env var is missing instead of silently
    defaulting to anything."""
    token = os.environ.get("IONQ_API_KEY")
    if not token:
        raise RuntimeError(
            "IONQ_API_KEY must be set in .env (see .env.example)."
        )
    return IonQProvider(token)


def list_backends(provider: IonQProvider) -> list[dict]:
    """List every backend actually available to this account right now.
    Returns the raw comparison data; does NOT pick one for you."""
    rows = []
    for backend in provider.backends():
        rows.append({
            "name": backend.name,
            "num_qubits": backend.num_qubits,
            "simulator": backend.name == "ionq_simulator",
            "available": backend.status(),
        })
    return rows


def print_backends(rows: list[dict]) -> None:
    if not rows:
        print("  No backends found (check API key / account access).")
        return
    print(f"\n  {'name':<20} {'qubits':<10} {'simulator':<12} {'available':<10}")
    print("  " + "-" * 54)
    for r in rows:
        print(f"  {r['name']:<20} {r['num_qubits']:<10} "
              f"{str(r['simulator']):<12} {str(r['available']):<10}")


def get_simulator(provider: IonQProvider):
    """Return the free cloud simulator backend. This project only ever
    targets ionq_simulator (free, no QPU time charged) -- ionq_qpu is
    real trapped-ion hardware and deliberately out of scope here."""
    return provider.get_backend("ionq_simulator")


def main():
    print("\n" + "=" * 70)
    print("  IonQ Cloud — API connection + live backend listing")
    print("=" * 70)
    provider = connect_provider()
    print("\n  Connected.")

    rows = list_backends(provider)
    print_backends(rows)
    print(f"\n  {len(rows)} backend(s) total.\n")
    return rows


if __name__ == "__main__":
    main()
