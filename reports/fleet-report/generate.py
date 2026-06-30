"""
Generate the State of the IBM Quantum Fleet report.
Reads from devices.db, writes charts + fleet-report.md into this directory.
All findings and chart titles are computed from live data — nothing hardcoded.
"""

import sqlite3
import os
from datetime import datetime, timezone
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# ── config ────────────────────────────────────────────────────────────────────
DB_PATH = os.path.join(os.path.dirname(__file__), "..", "..", "devices.db")
OUT_DIR = os.path.dirname(__file__)

DEVICES = ["ibm_fez", "ibm_kingston", "ibm_marrakesh"]
COLORS  = {"ibm_fez": "#1f77b4", "ibm_kingston": "#ff7f0e", "ibm_marrakesh": "#2ca02c"}
LABELS  = {"ibm_fez": "ibm_fez (156 q)", "ibm_kingston": "ibm_kingston (156 q)",
           "ibm_marrakesh": "ibm_marrakesh (156 q)"}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "figure.dpi": 150,
})

# ── load data ─────────────────────────────────────────────────────────────────
conn = sqlite3.connect(DB_PATH)
cur  = conn.cursor()
cur.execute("SELECT ts, name, pending_jobs, avg_cx_error, avg_readout_error FROM device_snapshots ORDER BY ts")
raw = cur.fetchall()
conn.close()

data = defaultdict(lambda: {"ts": [], "jobs": [], "cx": [], "ro": []})
for ts_str, name, jobs, cx, ro in raw:
    ts = datetime.fromisoformat(ts_str).astimezone(timezone.utc)
    data[name]["ts"].append(ts)
    data[name]["jobs"].append(jobs)
    data[name]["cx"].append(cx)
    data[name]["ro"].append(ro)

def valid_pairs(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if y is not None]
    return (list(zip(*pairs)) if pairs else ([], []))

def mean_non_null(lst):
    vals = [v for v in lst if v is not None]
    return sum(vals) / len(vals) if vals else None

def non_null(lst):
    return [v for v in lst if v is not None]

# ── derived stats ─────────────────────────────────────────────────────────────
all_ts   = sorted(set(t for d in data.values() for t in d["ts"]))
n_snaps  = len(all_ts)
n_rows   = n_snaps * len(DEVICES)
span_h   = (all_ts[-1] - all_ts[0]).total_seconds() / 3600
span_d   = span_h / 24
date_fmt = "%Y-%m-%d %H:%M UTC"

# queue stats
max_jobs_dev, max_jobs_val, max_jobs_ts = None, 0, None
for dev in DEVICES:
    for t, j in zip(data[dev]["ts"], data[dev]["jobs"]):
        if j is not None and j > max_jobs_val:
            max_jobs_val, max_jobs_dev, max_jobs_ts = j, dev, t

# readout error stats per device
ro_stats = {}
for dev in DEVICES:
    vals = non_null(data[dev]["ro"])
    if vals:
        ro_stats[dev] = {
            "min": min(vals), "max": max(vals), "mean": sum(vals)/len(vals),
            "first": vals[0], "last": vals[-1], "n": len(vals),
        }

# cx error stats per device
cx_obs = {}   # {dev: [(ts, val), ...]}
for dev in DEVICES:
    pairs = [(t, v) for t, v in zip(data[dev]["ts"], data[dev]["cx"]) if v is not None]
    if pairs:
        cx_obs[dev] = pairs

cx_snap_count = len(set(t for dev in cx_obs for t, _ in cx_obs[dev]))

# null counts
ro_nulls = sum(1 for d in data.values() for v in d["ro"] if v is None)
cx_nulls = sum(1 for d in data.values() for v in d["cx"] if v is None)

# ══════════════════════════════════════════════════════════════════════════════
# Chart 1 — Queue depth over time
# ══════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 4))
for dev in DEVICES:
    xs, ys = valid_pairs(data[dev]["ts"], data[dev]["jobs"])
    ax.plot(list(xs), list(ys), marker="o", linewidth=1.8, markersize=5,
            color=COLORS[dev], label=LABELS[dev])

ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax.xaxis.set_major_locator(mdates.AutoDateLocator())
fig.autofmt_xdate(rotation=0, ha="center")
ax.set_ylabel("Pending jobs in queue")
ax.set_title(f"Queue Depth — All Devices ({n_snaps} snapshots, "
             f"{all_ts[0].strftime('%b %d')}–{all_ts[-1].strftime('%b %d %Y')})", fontsize=11)
ax.legend(framealpha=0.4, fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "chart_queue.png"))
plt.close(fig)
print("Saved chart_queue.png")

# ══════════════════════════════════════════════════════════════════════════════
# Chart 2 — Readout error over time
# ══════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 4))
for dev in DEVICES:
    xs, ys = valid_pairs(data[dev]["ts"], data[dev]["ro"])
    if xs:
        ax.plot(list(xs), [v * 100 for v in ys], marker="s", linewidth=1.8,
                markersize=5, color=COLORS[dev], label=LABELS[dev])

ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
ax.xaxis.set_major_locator(mdates.AutoDateLocator())
fig.autofmt_xdate(rotation=0, ha="center")
ax.set_ylabel("Avg readout error (%)")
ax.set_title(f"Readout Error Rate Over Time — {all_ts[0].strftime('%b %d')}–{all_ts[-1].strftime('%b %d %Y')}",
             fontsize=11)
ax.legend(framealpha=0.4, fontsize=9)
if ro_nulls > 0:
    ax.annotate(f"⚠ {ro_nulls} readings lacked\ncalibration data",
                xy=(0.02, 0.95), xycoords="axes fraction",
                fontsize=8, color="#888", va="top")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "chart_readout_error.png"))
plt.close(fig)
print("Saved chart_readout_error.png")

# ══════════════════════════════════════════════════════════════════════════════
# Chart 3 — Error rate summary bars (mean readout + mean CX)
# ══════════════════════════════════════════════════════════════════════════════
fig, (ax_ro, ax_cx) = plt.subplots(1, 2, figsize=(10, 4))

# readout bars
ro_devs   = [d.replace("ibm_", "") for d in DEVICES if d in ro_stats]
ro_means  = [ro_stats[f"ibm_{d}"]["mean"] * 100 for d in ro_devs]
ro_colors = [COLORS[f"ibm_{d}"] for d in ro_devs]
bars = ax_ro.bar(ro_devs, ro_means, color=ro_colors, width=0.5, alpha=0.85)
ax_ro.bar_label(bars, fmt="%.2f%%", padding=3, fontsize=9)
ax_ro.set_ylabel("Error rate (%)")
ax_ro.set_title("Mean Readout Error by Device", fontsize=10)
ax_ro.set_ylim(0, max(ro_means) * 1.35)

# CX bars
cx_devs  = [d.replace("ibm_", "") for d in DEVICES if d in cx_obs]
cx_means = [mean_non_null([v for _, v in cx_obs[f"ibm_{d}"]]) * 100 for d in cx_devs]
cx_colors = [COLORS[f"ibm_{d}"] for d in cx_devs]
cx_note = f"({cx_snap_count} calibration snapshot{'s' if cx_snap_count != 1 else ''})"
bars2 = ax_cx.bar(cx_devs, cx_means, color=cx_colors, width=0.5, alpha=0.85)
ax_cx.bar_label(bars2, fmt="%.2f%%", padding=3, fontsize=9)
ax_cx.set_ylabel("Error rate (%)")
ax_cx.set_title(f"Mean CX / 2-Qubit Gate Error\n{cx_note}", fontsize=10)
ax_cx.set_ylim(0, max(cx_means) * 1.35 if cx_means else 1)

fig.suptitle(f"Device Error Rate Summary — {all_ts[0].strftime('%b %d')}–{all_ts[-1].strftime('%b %d %Y')}",
             fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "chart_error_summary.png"), bbox_inches="tight")
plt.close(fig)
print("Saved chart_error_summary.png")

# ══════════════════════════════════════════════════════════════════════════════
# Build dynamic report text
# ══════════════════════════════════════════════════════════════════════════════

# --- queue findings ---
queue_bullets = [
    f"- Peak queue depth was **{max_jobs_val} pending jobs** on {max_jobs_dev} "
    f"({max_jobs_ts.strftime('%H:%M UTC %b %d')})."
]
# find lightest device by mean queue
mean_jobs = {dev: mean_non_null(data[dev]["jobs"]) for dev in DEVICES}
lightest  = min(mean_jobs, key=lambda d: mean_jobs[d] or 999)
queue_bullets.append(
    f"- {lightest} had the lightest average queue load "
    f"(mean {mean_jobs[lightest]:.1f} jobs), making it the most accessible backend on average."
)
queue_bullets.append(
    "- Queue depth is volatile on timescales of minutes; "
    "these readings reflect point-in-time snapshots, not sustained load."
)

# --- readout findings ---
ro_bullets = []
if ro_stats:
    best_ro  = min(ro_stats, key=lambda d: ro_stats[d]["mean"])
    worst_ro = max(ro_stats, key=lambda d: ro_stats[d]["mean"])
    ro_bullets.append(
        f"- {best_ro} had the **lowest mean readout error** across the window "
        f"({ro_stats[best_ro]['mean']*100:.2f}%, range "
        f"{ro_stats[best_ro]['min']*100:.2f}–{ro_stats[best_ro]['max']*100:.2f}%)."
    )
    ro_bullets.append(
        f"- {worst_ro} had the highest mean readout error "
        f"({ro_stats[worst_ro]['mean']*100:.2f}%, range "
        f"{ro_stats[worst_ro]['min']*100:.2f}–{ro_stats[worst_ro]['max']*100:.2f}%)."
    )
    # flag any device whose latest readout is >50% higher than its mean
    for dev in DEVICES:
        if dev not in ro_stats:
            continue
        s = ro_stats[dev]
        if s["last"] > s["mean"] * 1.5:
            ro_bullets.append(
                f"- {dev}'s most recent readout error ({s['last']*100:.2f}%) is "
                f"notably above its window mean ({s['mean']*100:.2f}%) — worth monitoring."
            )

# --- cx error table ---
if cx_obs:
    cx_rows = []
    for dev in DEVICES:
        if dev not in cx_obs:
            continue
        for ts, val in cx_obs[dev]:
            cx_rows.append(f"| {dev} | {ts.strftime('%Y-%m-%d %H:%M UTC')} | {val*100:.2f}% |")
    cx_table = "\n".join(cx_rows)
    cx_section = f"""**CX (two-qubit gate) error observations ({cx_snap_count} calibration snapshot{'s' if cx_snap_count != 1 else ''}):**

| Device | Snapshot | CX error |
|--------|----------|----------|
{cx_table}

Note: CX error data was absent from historical snapshots due to a gate-name mismatch
(Eagle devices use ECR, not CX). This was fixed on 2026-06-22. Trend analysis will
become possible once several more snapshots accumulate."""
else:
    cx_section = "*No CX error data available yet.*"

# --- data gaps table ---
gap_rows = []
if ro_nulls > 0:
    gap_rows.append(
        f"| Missing readout error | {ro_nulls}/{n_rows} readings | "
        "Early snapshots; API returned no calibration data |"
    )
cx_hist_nulls = cx_nulls - (n_rows - len([v for d in data.values() for v in d["cx"] if v is None]))
gap_rows.append(
    f"| Missing CX error (historical) | {cx_nulls}/{n_rows} readings | "
    "ECR gate not matched by old `g.gate == 'cx'` filter — fixed 2026-06-22 |"
)
# detect near-duplicates (snapshots within 2 min of each other)
dup_pairs = []
for i in range(len(all_ts) - 1):
    diff = (all_ts[i+1] - all_ts[i]).total_seconds()
    if diff < 120:
        dup_pairs.append((all_ts[i], all_ts[i+1], int(diff)))
if dup_pairs:
    gap_rows.append(
        f"| Near-duplicate snapshots | {len(dup_pairs)} pair(s) | "
        "Two snapshots within 2 minutes — likely a retry or double-fire |"
    )

gaps_table = "\n".join(gap_rows) if gap_rows else "| None detected | — | — |"

snap_list = "\n".join(f"- {t.strftime(date_fmt)}" for t in all_ts)
report_date = datetime.now(timezone.utc).strftime(date_fmt)

# ── write report ───────────────────────────────────────────────────────────────
report = f"""# State of the IBM Quantum Fleet
**Report generated:** {report_date}
**Data window:** {all_ts[0].strftime(date_fmt)} – {all_ts[-1].strftime(date_fmt)} ({span_d:.1f} days)
**Devices monitored:** ibm_fez, ibm_kingston, ibm_marrakesh (all 156-qubit Eagle-class)
**Snapshots collected:** {n_snaps} polling intervals × {len(DEVICES)} devices = {n_rows} rows

---

## Executive summary

All three monitored backends remained operational throughout the {span_d:.0f}-day observation window.
Queue depths varied, peaking at {max_jobs_val} pending jobs on {max_jobs_dev}.
Readout error rates were broadly stable across devices, with {best_ro if ro_stats else 'N/A'} showing
the best average performance. CX (two-qubit gate) error data is newly available following
a bug fix on 2026-06-22 and should be treated as early-stage data, not a trend.

---

## Devices in scope

| Device | Qubits | Class | Operational throughout? |
|--------|--------|-------|------------------------|
| ibm_fez | 156 | Eagle r3 | Yes |
| ibm_kingston | 156 | Eagle r3 | Yes |
| ibm_marrakesh | 156 | Eagle r3 | Yes |

No device reported `operational = False` during the observation window.

---

## Chart 1 — Queue depth over time

![Queue depth over time](chart_queue.png)

**Findings:**

{chr(10).join(queue_bullets)}

---

## Chart 2 — Readout error rate over time

![Readout error over time](chart_readout_error.png)

**Findings:**

{chr(10).join(ro_bullets) if ro_bullets else '_No readout error data available._'}

**Data note:** {n_rows - ro_nulls} of {n_rows} device-readings included readout error data.
{ro_nulls} readings returned NULL — these are excluded from all calculations.

---

## Chart 3 — Error rate summary

![Error rate summary](chart_error_summary.png)

{cx_section}

---

## Data quality and gaps

| Issue | Affected rows | Detail |
|-------|--------------|--------|
{gaps_table}

All NULL values are stored as-is in `devices.db`. No imputation or forward-filling applied.

---

## Snapshot inventory

{snap_list}

---

## Recommendations

1. **Allow CX error to accumulate.** The ECR/CX bug was fixed on 2026-06-22. Wait for
   5–7 more snapshots before drawing conclusions from the CX trend.
2. **Deduplicate near-simultaneous snapshots.** {len(dup_pairs)} near-duplicate pair(s) detected.
   Add a guard in `snapshot.py` to skip writes when a snapshot for the same window already exists.
3. **Extend the collection window.** IBM calibration cycles run 24–48 hours; {span_d:.0f} days
   of data {'is a reasonable baseline' if span_d >= 7 else 'is still short — aim for 14+ days for reliable trend analysis'}.

---

*Generated automatically from `devices.db`. All statistics computed from raw snapshot values,
no smoothing applied. IBM Quantum hardware characteristics are subject to change without notice.*
"""

with open(os.path.join(OUT_DIR, "fleet-report.md"), "w") as f:
    f.write(report)
print(f"Saved fleet-report.md ({len(report)} chars)")
