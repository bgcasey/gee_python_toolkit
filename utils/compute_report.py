# ---
# title:   GEE Compute Usage Report
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Collects Earth Engine compute usage information and
#   writes it to a plain-text report. Two signals are
#   captured:
#
#   1. Per-algorithm EECU profiles (ee.profilePrinting)
#      for interactive computations, i.e., anything that
#      forces evaluation such as reduceRegion().getInfo().
#      Rows with the highest EECU-s are the compute choke
#      points.
#   2. Total batch EECU-seconds for export tasks, read
#      from the task status after completion.
#
#   Earth Engine is lazy: building an image expression
#   uses no compute. Only evaluated results and batch
#   tasks consume EECUs, so profile sections must contain
#   a getInfo()-style call to produce output. To find
#   choke points cheaply, run the pipeline on a small
#   test AOI with profiling on, then extrapolate.
# ---

import datetime
import io
import os
import time
from contextlib import contextmanager

import ee


class ComputeReport:
    """Collect EE compute usage and write a txt report.

    Args:
        name (str): Report name, used in the output file
            name (e.g., the script name).
        out_dir (str): Directory for the report file.
            Defaults to 'gee_compute_reports' in the
            current working directory.
        enabled (bool): If False, all methods are no-ops
            so calling code needs no conditionals.
    """

    def __init__(self, name, out_dir=None, enabled=True):
        self.name = name
        self.enabled = enabled
        if out_dir is None:
            out_dir = os.path.join(
                os.getcwd(), "gee_compute_reports"
            )
        self.out_dir = out_dir
        self._blocks = []

    @contextmanager
    def section(self, name):
        """Profile a block of code and record EECU usage.

        Wraps the block in ee.profilePrinting so every
        computation evaluated inside it (getInfo, etc.)
        is profiled per algorithm.

        Args:
            name (str): Section label used in the report.
        """
        if not self.enabled:
            yield
            return

        buf = io.StringIO()
        start = time.time()
        try:
            with ee.profilePrinting(destination=buf):
                yield
        finally:
            elapsed = time.time() - start
            profile = buf.getvalue().strip()
            if not profile:
                profile = (
                    "No server-side computation was "
                    "evaluated in this section. Add a "
                    "getInfo()-style call to profile it."
                )
            self._blocks.append(
                f"--- Section: {name} ---\n"
                f"Wall time: {elapsed:.1f} s\n"
                f"EECU profile (highest compute first):\n"
                f"{profile}\n"
            )

    def log_task(self, task, poll_interval=30):
        """Wait for an export task and record its EECU use.

        Blocks until the task finishes. Batch EECU totals
        are only available on completed tasks.

        Args:
            task (ee.batch.Task): A started export task.
            poll_interval (int): Seconds between checks.
        """
        if not self.enabled:
            return

        description = task.config.get(
            "description", task.id
        )
        print(
            f"Waiting for task '{description}' to record "
            f"compute usage..."
        )
        while task.active():
            time.sleep(poll_interval)

        status = task.status()
        eecu = status.get("batch_eecu_usage_seconds")
        start_ms = status.get("start_timestamp_ms")
        update_ms = status.get("update_timestamp_ms")
        runtime = (
            f"{(update_ms - start_ms) / 1000:.0f} s"
            if start_ms and update_ms
            else "unknown"
        )
        eecu_text = (
            f"{eecu:.1f}" if eecu is not None else "unavailable"
        )

        block = (
            f"--- Export task: {description} ---\n"
            f"State: {status['state']}\n"
            f"Runtime: {runtime}\n"
            f"Batch EECU-seconds: {eecu_text}\n"
        )
        if status["state"] == "FAILED":
            block += (
                f"Error: {status.get('error_message')}\n"
            )
        self._blocks.append(block)

    def write(self):
        """Write collected blocks to a timestamped txt file.

        Returns:
            str: Path to the report file, or None if the
            report is disabled or empty.
        """
        if not self.enabled or not self._blocks:
            return None

        os.makedirs(self.out_dir, exist_ok=True)
        timestamp = datetime.datetime.now().strftime(
            "%Y%m%d_%H%M%S"
        )
        path = os.path.join(
            self.out_dir,
            f"{self.name}_compute_{timestamp}.txt",
        )

        header = (
            f"{'=' * 60}\n"
            f"GEE Compute Usage Report: {self.name}\n"
            f"Generated: "
            f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S}\n"
            f"{'=' * 60}\n\n"
            "How to read this report:\n"
            "- EECU-s = Earth Engine Compute Unit seconds.\n"
            "- In section profiles, rows are sorted by\n"
            "  compute; the top rows are the choke points.\n"
            "- 'Count' is how many times an algorithm ran;\n"
            "  high counts suggest repeated work that\n"
            "  could be cached or restructured.\n"
            "- Batch EECU-seconds is the total compute\n"
            "  consumed by an export task.\n\n"
        )

        with open(path, "w") as f:
            f.write(header)
            f.write("\n".join(self._blocks))

        print(f"Compute report written to: {path}")
        return path
