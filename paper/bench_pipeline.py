"""Compare base-depth mate-fix pipelines on identical input and output settings.

Build the supplied binaries first. This measures whole processes, not just
pileup generation. All outputs must match before a run is accepted.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bam", type=Path)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--pipeline", type=Path, required=True)
    parser.add_argument("--fork", type=Path,
                        help="Optional pipeline binary using the lean Seqair fork")
    parser.add_argument("--goose", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True,
                        help="New directory for timings, logs and results.json")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    programs = {"optimized rust + seqair": args.baseline.resolve(),
                "rust + goose pipeline": args.pipeline.resolve()}
    if args.fork:
        programs["rust + goose pipeline + seqair fork"] = args.fork.resolve()
    programs["goose v0.1.3"] = args.goose.resolve()
    report = {"timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "platform": platform.platform(), "affinity": sorted(os.sched_getaffinity(0)),
              "input": str(args.bam.resolve()), "input_sha256": digest(args.bam),
              "input_bytes": args.bam.stat().st_size, "runs": [],
              "note": "Two sequential runs each, reversed order on round two; no cache dropping. "
                      "Input hashing precedes timing. Output hashing follows timing. "
                      "Binary provenance and compiler flags must be recorded by the caller."}
    hashes = {name: digest(binary) for name, binary in programs.items()}
    expected = None
    for round_ in range(2):
        for name in programs if round_ == 0 else reversed(programs):
            stem = args.output / f"{len(report['runs']) + 1}"
            output = stem.with_suffix(".tsv")
            command = [str(programs[name]), "base-depth", str(args.bam.resolve()),
                       "-m", "-t", "32", "-c", "1000000", "-o", str(output)]
            if name != "goose v0.1.3":
                command.append("--seqair")
            with stem.with_suffix(".log").open("w") as log:
                start = time.perf_counter()
                subprocess.run(["/usr/bin/time", "-f", "%e %U %S %M", "-o",
                                str(stem.with_suffix(".time")), *command],
                               stderr=log, check=True)
                elapsed = time.perf_counter() - start
            wall, user, system, rss = map(float, stem.with_suffix(".time").read_text().split())
            checksum = digest(output)
            if expected is None:
                expected = checksum
            run = dict(stage=name, command=command, wall_s=elapsed, time_wall_s=wall,
                       user_s=user, system_s=system, peak_rss_kib=int(rss),
                       binary_sha256=hashes[name], output_sha256=checksum,
                       output_bytes=output.stat().st_size)
            report["runs"].append(run)
            (args.output / "results.json").write_text(json.dumps(report, indent=2) + "\n")
            if checksum != expected:
                raise RuntimeError(f"Output mismatch: {name}; retained {output}")
            output.unlink()
            print(f"{name}: {elapsed:.2f}s, {rss / 1024**2:.2f} GiB; output matched", flush=True)
    for name in programs:
        runs = [run for run in report["runs"] if run["stage"] == name]
        print(name, "mean seconds:", statistics.mean(run["wall_s"] for run in runs),
              "mean peak GiB:", statistics.mean(run["peak_rss_kib"] / 1024**2 for run in runs))


if __name__ == "__main__":
    main()
