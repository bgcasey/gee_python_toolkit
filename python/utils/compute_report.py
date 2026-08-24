# ---
# title:   GEE Compute Usage Report
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Collects Earth Engine compute usage information and
#   writes it to a plain-text report. Three signals are
#   captured:
#
#   1. Per-algorithm EECU profiles (ee.profilePrinting)
#      for interactive computations, i.e., anything that
#      forces evaluation such as reduceRegion().getInfo().
#      Rows with the highest EECU-s are the compute choke
#      points.
#   2. Total batch EECU-seconds for export tasks, read
#      from the task status after completion.
#   3. Warnings and errors for troubleshooting: exceptions
#      raised inside a profiled section (with traceback),
#      Python warnings emitted there (e.g., EE deprecation
#      notices), failed export tasks, and manual notes added
#      with log_warning(). All are echoed to the console as
#      they happen and collected into an "Issues detected"
#      summary at the top of the report.
#
#   Only evaluated results and batch tasks consume
#   EECUs, so profile sections must contain a
#   getInfo()-style call to produce output. To find
#   choke points cheaply, run the pipeline on a
#   small test AOI with profiling on, then
#   extrapolate.
# ---

import datetime
import io
import os
import sys
import time
import traceback
import warnings
from contextlib import contextmanager

import ee


# Longest captured section output kept in the report, in
# characters. Sections that dump whole feature collections can
# run to megabytes; the tail is dropped with a marker.
MAX_SECTION_OUTPUT = 20000


class _Tee:
    """Write to several streams at once.

    Used to mirror a section's stdout into the report while it
    still reaches the console live.
    """

    def __init__(self, *streams):
        self._streams = streams

    def write(self, text):
        for stream in self._streams:
            stream.write(text)
        return len(text)

    def flush(self):
        for stream in self._streams:
            stream.flush()


class ComputeReport:
    """Collect EE compute usage and write a txt report.

    Besides EECU profiles, the report captures anything that
    helps troubleshoot a run: exceptions raised inside a
    profiled section, Python warnings emitted there (e.g., EE
    deprecation notices), failed export tasks, and any manual
    notes added with log_warning(). All of these are echoed to
    the console as they happen and collected into an "Issues
    detected" summary at the top of the report file.

    Args:
        name (str): Report name, used in the output file
            name (e.g., the script name).
        out_dir (str): Directory for the report file.
            Defaults to the repo-root 'gee_compute_reports'
            directory, resolved from this module's location
            so it is independent of the working directory
            VS Code runs a script from.
        enabled (bool): If False, all methods are no-ops
            so calling code needs no conditionals.
    """

    def __init__(self, name, out_dir=None, enabled=True):
        self.name = name
        self.enabled = enabled
        if out_dir is None:
            # Repo root is two levels up from this module
            # (python/utils/compute_report.py), so reports
            # always land in the project's top-level
            # gee_compute_reports/ regardless of cwd.
            repo_root = os.path.dirname(
                os.path.dirname(
                    os.path.dirname(os.path.abspath(__file__))
                )
            )
            out_dir = os.path.join(
                repo_root, "gee_compute_reports"
            )
        self.out_dir = out_dir
        self._blocks = []
        # Short, high-signal lines surfaced in the report
        # header so failures/warnings are easy to spot.
        self._issues = []

    def _record_issue(self, level, message):
        """Note a warning/error and echo it to the console.

        Args:
            level (str): "WARNING" or "ERROR".
            message (str): One-line summary of the issue.
        """
        line = f"[{level}] {message}"
        self._issues.append(line)
        print(line, file=sys.stderr)

    def log_warning(self, message):
        """Record a manual troubleshooting note.

        Use for conditions the script detects itself (empty
        collection, suspicious stats, skipped step) that
        would not otherwise raise. No-op when disabled.

        Args:
            message (str): The note to record.
        """
        if not self.enabled:
            return
        self._record_issue("WARNING", message)
        self._blocks.append(f"--- Warning ---\n{message}\n")

    @contextmanager
    def section(self, name, raise_on_error=True):
        """Profile a block of code and record EECU usage.

        Wraps the block in ee.profilePrinting so every
        computation evaluated inside it (getInfo, etc.)
        is profiled per algorithm.

        Python warnings emitted in the block are captured and
        recorded. If the block raises (e.g., an EEException
        from a getInfo call), the error and traceback are
        recorded against this section and echoed to the
        console.

        Args:
            name (str): Section label used in the report.
            raise_on_error (bool): If True (default), a failure
                *in the profiled block* is re-raised after being
                recorded, so the script still fails loudly. A
                failure to retrieve the profile itself is never
                re-raised: it is recorded as a warning and the
                run continues. Set False for
                optional diagnostic blocks (e.g., a min/max
                preview) whose failure should not abort the
                run; the error is recorded and execution
                continues past the block. When False, code
                after the block must tolerate a result the
                block never assigned (initialize it first).
        """
        if not self.enabled:
            yield
            return

        buf = io.StringIO()
        printed_buf = io.StringIO()
        start = time.time()
        status = "OK"
        error_text = None
        # Distinguishes a failure in the caller's block from one
        # raised by ee.profilePrinting on exit. The assignment is
        # reached only when the block completed normally: if the
        # block raises, the exception is thrown in at the yield.
        body_ok = False
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            # Mirror the block's stdout into the report so the
            # numbers a section prints (band min/max, band names,
            # ...) are kept alongside its EECU profile, instead
            # of only scrolling past in the terminal.
            real_stdout = sys.stdout
            sys.stdout = _Tee(real_stdout, printed_buf)
            try:
                with ee.profilePrinting(destination=buf):
                    yield
                    body_ok = True
            except Exception as exc:
                if body_ok:
                    # The profiled code worked; only the
                    # profile fetch failed. Earth Engine can hand
                    # back a profile ID that getProfiles() will
                    # not resolve ("Profile not found"). Observed
                    # intermittently in sections that evaluate
                    # many separate getInfo calls - it is not a
                    # timeout (seen in a 6 s / 14-call section,
                    # while slower single-call sections were
                    # fine), and which section trips it varies
                    # between runs. Batching a section into one
                    # call avoids it. Either way this is a
                    # diagnostics problem, never a reason to
                    # abort a run whose real work succeeded.
                    status = "OK (profile unavailable)"
                    self._record_issue(
                        "WARNING",
                        f"Section '{name}' ran successfully but "
                        f"its EECU profile could not be "
                        f"retrieved: {type(exc).__name__}: {exc}",
                    )
                else:
                    status = "FAILED"
                    error_text = (
                        f"{type(exc).__name__}: {exc}\n"
                        f"{traceback.format_exc()}"
                    )
                    self._record_issue(
                        "ERROR",
                        f"Section '{name}' failed: "
                        f"{type(exc).__name__}: {exc}",
                    )
                    if raise_on_error:
                        raise
            finally:
                sys.stdout = real_stdout
                elapsed = time.time() - start
                profile = buf.getvalue().strip()
                if not profile:
                    if status.startswith("OK (profile"):
                        profile = (
                            "Unavailable: the section ran but "
                            "Earth Engine would not return its "
                            "profile. Seen intermittently in "
                            "sections that evaluate many "
                            "separate getInfo calls; batch them "
                            "into one call to profile it."
                        )
                    else:
                        profile = (
                            "No server-side computation was "
                            "evaluated in this section. Add a "
                            "getInfo()-style call to profile it."
                        )
                block = (
                    f"--- Section: {name} ---\n"
                    f"Status: {status}\n"
                    f"Wall time: {elapsed:.1f} s\n"
                )
                printed = printed_buf.getvalue().strip()
                if printed:
                    if len(printed) > MAX_SECTION_OUTPUT:
                        dropped = len(printed) - MAX_SECTION_OUTPUT
                        printed = (
                            printed[:MAX_SECTION_OUTPUT]
                            + f"\n... [{dropped} more characters "
                            "truncated]"
                        )
                    block += f"Output:\n{printed}\n"
                for w in caught:
                    warn_msg = f"{w.category.__name__}: {w.message}"
                    self._record_issue(
                        "WARNING",
                        f"Section '{name}': {warn_msg}",
                    )
                    block += f"Warning: {warn_msg}\n"
                if error_text is not None:
                    block += f"Error:\n{error_text}\n"
                block += (
                    f"EECU profile (highest compute first):\n"
                    f"{profile}\n"
                )
                self._blocks.append(block)

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
            error_message = status.get("error_message")
            block += f"Error: {error_message}\n"
            self._record_issue(
                "ERROR",
                f"Export '{description}' failed: "
                f"{error_message}",
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
            "  consumed by an export task.\n"
            "- Issues detected (below) lists warnings and\n"
            "  errors captured during the run; each points\n"
            "  to the section or export where it occurred.\n\n"
        )

        if self._issues:
            header += (
                f"{'-' * 60}\n"
                f"Issues detected ({len(self._issues)}):\n"
            )
            header += "\n".join(self._issues) + "\n"
            header += f"{'-' * 60}\n\n"
        else:
            header += "No warnings or errors were captured.\n\n"

        with open(path, "w") as f:
            f.write(header)
            f.write("\n".join(self._blocks))

        print(f"Compute report written to: {path}")
        return path
