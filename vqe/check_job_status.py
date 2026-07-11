#!/usr/bin/env python3
"""Quick status check of recent jobs in the Azure Quantum workspace."""
import os
import sys
from dotenv import load_dotenv
load_dotenv()

sys.path.insert(0, os.path.dirname(__file__))
from azure_backend import connect_workspace


def main():
    ws = connect_workspace()
    jobs = list(ws.list_jobs())
    jobs.sort(key=lambda j: j.details.creation_time, reverse=True)
    print(f"\n  {len(jobs)} job(s) in workspace, most recent first:\n")
    for j in jobs[:10]:
        print(f"  id={j.id}  target={j.details.target}  status={j.details.status}  "
              f"created={j.details.creation_time}")


if __name__ == "__main__":
    main()
