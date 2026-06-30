"""
report.py
---------
Generates the daily Quantum Weatherman report.

Reads the last 7 days of snapshots from devices.db, writes a short
plain-English summary, generates a trend chart as PNG, and saves both
to reports/YYYY-MM-DD/. Also appends the summary to REPORTS.md.

Run manually:
    .venv/bin/python report.py

The LaunchAgent fires this automatically at 8am every day.
"""

import os
import sqlite3
from collections import defaultdict
from datetime import datetime, timezone

import matplotlib
matplotlib.use("Agg")   # no display needed — we're saving to a file
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

BASE_DIR  = os.path.dirname(__file__)
DB_PATH   = os.path.join(BASE_DIR, "devices.db")
REPORT_DIR = os.path.join(BASE_DIR, "reports")
REPORTS_MD = os.path.join(BASE_DIR, "REPORTS.md")


# --------------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------------

def load_snapshots(days: int = 7) -> list[dict]:
    """Pull all rows from the last `days` days, oldest first."""
    with sqlite3.connect(DB_PATH) as con:
        con.row_factory = sqlite3.Row
        rows = con.execute(
            """
            SELECT ts, name, num_qubits, operational,
                   pending_jobs, avg_cx_error, avg_readout_error
            FROM   device_snapshots
            WHERE  ts >= datetime('now', ? || ' days')
            ORDER  BY ts ASC
            """,
            (f"-{days}",),
        ).fetchall()
    return [dict(r) for r in rows]


def group_by_device(snapshots: list[dict]) -> dict[str, list[dict]]:
    groups = defaultdict(list)
    for s in snapshots:
        groups[s["name"]].append(s)
    return dict(groups)


# --------------------------------------------------------------------------
# Analysis helpers
# --------------------------------------------------------------------------

def daily_averages(snaps: list[dict]) -> dict[str, dict]:
    """
    For one device's snapshots, compute daily averages of queue depth
    and CX error. Returns {date_str: {pending_jobs, avg_cx_error}}.
    """
    by_day = defaultdict(list)
    for s in snaps:
        by_day[s["ts"][:10]].append(s)   # "YYYY-MM-DD"

    result = {}
    for day, rows in sorted(by_day.items()):
        pending = [r["pending_jobs"] for r in rows if r["pending_jobs"] is not None]
        cx      = [r["avg_cx_error"] for r in rows if r["avg_cx_error"] is not None]
        result[day] = {
            "pending_jobs":  sum(pending) / len(pending) if pending else None,
            "avg_cx_error":  sum(cx) / len(cx)           if cx      else None,
        }
    return result


def best_machine(by_device: dict[str, list[dict]]) -> tuple[str, dict]:
    """
    Rank devices by their most recent snapshot.
    Primary sort: lowest CX error. Tiebreak: shortest queue.
    """
    latest = {name: snaps[-1] for name, snaps in by_device.items()}
    ranked = sorted(
        latest.items(),
        key=lambda kv: (
            kv[1]["avg_cx_error"]  if kv[1]["avg_cx_error"]  is not None else float("inf"),
            kv[1]["pending_jobs"]  if kv[1]["pending_jobs"]   is not None else float("inf"),
        ),
    )
    return ranked[0] if ranked else ("unknown", {})


def compute_trends(by_device: dict[str, list[dict]]) -> dict[str, dict]:
    """
    For each device compare the first and last daily average over the window.
    Returns per-device dicts with cx_change and q_change.
    """
    trends = {}
    for name, snaps in by_device.items():
        avgs = daily_averages(snaps)
        days = sorted(avgs.keys())
        if len(days) < 2:
            continue
        first, last = avgs[days[0]], avgs[days[-1]]
        trends[name] = {
            "first_cx": first["avg_cx_error"],
            "last_cx":  last["avg_cx_error"],
            "cx_change": (
                last["avg_cx_error"] - first["avg_cx_error"]
                if first["avg_cx_error"] and last["avg_cx_error"] else None
            ),
            "first_q": first["pending_jobs"],
            "last_q":  last["pending_jobs"],
            "q_change": (
                last["pending_jobs"] - first["pending_jobs"]
                if first["pending_jobs"] is not None and last["pending_jobs"] is not None
                else None
            ),
        }
    return trends


# --------------------------------------------------------------------------
# Report text
# --------------------------------------------------------------------------

def build_report_text(today: str, top_name: str, top_snap: dict,
                      trends: dict) -> str:
    lines = [f"## Quantum Weatherman — {today}\n"]

    # Best machine
    cx = top_snap.get("avg_cx_error")
    q  = top_snap.get("pending_jobs", "?")
    cx_str = f", avg CX error {cx:.4f}" if cx else ""
    lines.append(f"**Best machine today:** {top_name} ({q} jobs in queue{cx_str})\n")

    # Classify each device as improved / worsened / steady
    improved, worsened, steady = [], [], []
    for name, t in sorted(trends.items()):
        cx_c = t["cx_change"]
        q_c  = t["q_change"]

        got_better = cx_c is not None and t["first_cx"] and cx_c < -0.05 * t["first_cx"]
        got_worse  = (
            (cx_c is not None and t["first_cx"] and cx_c > 0.10 * t["first_cx"])
            or (q_c is not None and q_c > 20)
        )

        if got_better:
            improved.append(
                f"{name} (CX error {t['first_cx']:.4f} → {t['last_cx']:.4f})"
            )
        elif got_worse:
            detail = []
            if cx_c is not None and t["first_cx"]:
                detail.append(f"CX error {t['first_cx']:.4f} → {t['last_cx']:.4f}")
            if q_c is not None and q_c > 0:
                detail.append(f"queue +{int(q_c)}")
            worsened.append(f"{name} ({', '.join(detail)})" if detail else name)
        else:
            steady.append(name)

    if improved:
        lines.append(f"**Improved:** {', '.join(improved)}\n")
    if worsened:
        lines.append(f"**Getting worse:** {', '.join(worsened)}\n")
    if steady:
        lines.append(f"**Steady:** {', '.join(steady)}\n")

    lines.append("![Trend chart](chart.png)\n")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# Chart
# --------------------------------------------------------------------------

def generate_chart(by_device: dict[str, list[dict]], output_path: str) -> None:
    """
    Two-panel chart:
      Top    — queue depth (pending_jobs) per device over 7 days
      Bottom — avg CX error per device over 7 days
    """
    fig, (ax_q, ax_cx) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    fig.suptitle("IBM Quantum — 7-Day Trends", fontsize=14, fontweight="bold")

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    has_cx_data = False

    for idx, (name, snaps) in enumerate(sorted(by_device.items())):
        avgs  = daily_averages(snaps)
        days  = sorted(avgs.keys())
        dates = [datetime.strptime(d, "%Y-%m-%d") for d in days]
        color = colors[idx % len(colors)]

        # Queue depth panel
        pending = [avgs[d]["pending_jobs"] for d in days]
        if any(p is not None for p in pending):
            ax_q.plot(dates, pending, marker="o", label=name,
                      color=color, linewidth=2, markersize=5)

        # CX error panel
        cx = [avgs[d]["avg_cx_error"] for d in days]
        if any(c is not None for c in cx):
            has_cx_data = True
            ax_cx.plot(dates, cx, marker="o", label=name,
                       color=color, linewidth=2, markersize=5)

    ax_q.set_ylabel("Pending Jobs")
    ax_q.set_title("Queue Depth", fontsize=11)
    ax_q.legend(fontsize=9)
    ax_q.grid(True, alpha=0.3)
    ax_q.set_ylim(bottom=0)

    ax_cx.set_ylabel("Avg CX Error")
    ax_cx.set_title("Gate Error Rate  (lower = better)", fontsize=11)
    ax_cx.grid(True, alpha=0.3)
    # AutoDateLocator picks the right tick interval whether we have 1 day or 7.
    ax_cx.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax_cx.xaxis.set_major_formatter(mdates.ConciseDateFormatter(mdates.AutoDateLocator()))
    plt.xticks(rotation=30, ha="right")

    if has_cx_data:
        ax_cx.legend(fontsize=9)
    else:
        ax_cx.text(0.5, 0.5, "No CX error data yet\n(run snapshot.py a few times)",
                   transform=ax_cx.transAxes, ha="center", va="center",
                   color="grey", fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

def main() -> None:
    today = datetime.now(timezone.utc).strftime("%Y-%m-%d")

    snapshots = load_snapshots(days=7)
    if not snapshots:
        print("No snapshots in the last 7 days — run snapshot.py first.")
        return

    by_device   = group_by_device(snapshots)
    trends      = compute_trends(by_device)
    top_name, top_snap = best_machine(by_device)

    # Create today's output folder
    out_dir = os.path.join(REPORT_DIR, today)
    os.makedirs(out_dir, exist_ok=True)

    # Write report.md
    report_text = build_report_text(today, top_name, top_snap, trends)
    with open(os.path.join(out_dir, "report.md"), "w") as f:
        f.write(report_text)

    # Generate chart.png
    generate_chart(by_device, os.path.join(out_dir, "chart.png"))

    # Append to REPORTS.md running log
    with open(REPORTS_MD, "a") as f:
        f.write("\n---\n\n")
        f.write(report_text)

    print(
        f"[{datetime.now(timezone.utc).isoformat()}] "
        f"Report saved → {out_dir}"
    )


if __name__ == "__main__":
    main()
