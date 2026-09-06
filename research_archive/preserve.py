#!/usr/bin/env python3
"""Inventory and append-only preservation of noncanonical research files.

Bulk copying is deliberately non-destructive. Original working files stay put.
Run with the repository venv; --save creates a new dated checkpoint. A manifest
classifies every discovered untracked/ignored file, including exclusions.
This utility is a preservation aid, not a mathematical verifier or secret scanner
with a completeness guarantee. Review its summary before publishing.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
EXCLUDED_PARTS = {".venv", "pdfenv", "__pycache__", ".pytest_cache",
                  ".mypy_cache", "pip-cache", "site-packages", "python-wheels"}
SECRET_PATTERNS = {
    "github_token": rb"\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,})\b",
    "openai_style_key": rb"\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{35,}\b",
    "aws_access_key": rb"\b(?:AKIA|ASIA)[A-Z0-9]{16}\b",
    "private_key": rb"-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----",
}
SECRETS = {name: re.compile(pattern) for name, pattern in SECRET_PATTERNS.items()}
TMP_REF = re.compile(r"(?:/home/math/quadra/)?(?:computations/)?tmp/[A-Za-z0-9_.+/-]+")


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


def paths(*args):
    return {p.decode("utf-8", "surrogateescape") for p in git(*args).split(b"\0") if p}


def classify(name, data_head):
    p = Path(name)
    if p.parts[0] == "research_archive":
        return "archive_self", "append-only archives are not recursively copied"
    if any(part in EXCLUDED_PARTS for part in p.parts):
        return "environment_or_cache", "environment, downloaded software, or regenerable cache"
    if p.suffix.lower() in {".pyc", ".pyo", ".pyd", ".o", ".so", ".a", ".whl", ".aux", ".toc", ".synctex"}:
        return "build_product", "compiled/package product; source must be retained separately"
    if data_head.startswith(b"\x7fELF") or data_head[:4] in {
            b"\xfe\xed\xfa\xce", b"\xce\xfa\xed\xfe", b"\xfe\xed\xfa\xcf", b"\xcf\xfa\xed\xfe"}:
        return "build_product", "compiled executable (magic bytes), not an output log"
    if p.name in {".env", ".netrc", ".pypirc", "credentials", "credentials.json"}:
        return "credential_candidate", "credential-like filename; do not publish without review"
    return "research", "unreviewed research/sources/outputs retained without requiring polishing"


def inventory():
    untracked = paths("ls-files", "--others", "--exclude-standard", "-z")
    ignored = paths("ls-files", "--others", "--ignored", "--exclude-standard", "-z")
    modified = paths("diff", "HEAD", "--name-only", "--diff-filter=AM", "-z")
    records = []
    for name in sorted(untracked | ignored | modified):
        if name.startswith("research_archive/"):
            continue
        path = ROOT / name
        if not path.exists() and not path.is_symlink():
            continue
        record = {"path": name, "origin": "modified_tracked" if name in modified
                  else "ignored" if name in ignored else "untracked"}
        if path.is_symlink():
            record.update(category="symlink_reference", size=0,
                          reason="not followed; target recorded, not copied", target=str(path.readlink()))
            records.append(record)
            continue
        if not path.is_file():
            continue
        stat = path.stat()
        record.update(size=stat.st_size, mtime_ns=stat.st_mtime_ns)
        with path.open("rb") as stream:
            head = stream.read(128)
        category, reason = classify(name, head)
        record.update(category=category, reason=reason)
        if category == "research":
            data = path.read_bytes()
            flags = [label for label, pattern in SECRETS.items() if pattern.search(data)]
            if flags:
                record.update(category="credential_candidate", reason="strong secret-pattern match",
                              secret_pattern_names=flags)
            else:
                record["sha256"] = hashlib.sha256(data).hexdigest()
                record["evidence_status"] = "ARCHIVAL SNAPSHOT: status inherited from source; not newly verified"
                if len(data) > 25 * 1024 * 1024:
                    record.update(category="oversized_research", reason="requires reviewed durable storage; no silent omission")
        records.append(record)
    return records


def dependencies(records):
    by_path = {r["path"]: r for r in records}
    tracked = paths("ls-files", "-z")
    out = []
    for name in sorted(tracked):
        if name.startswith("research_archive/") or Path(name).suffix not in {".py", ".cpp", ".cc", ".h", ".md", ".sh"}:
            continue
        path = ROOT / name
        if not path.is_file():
            continue
        content = path.read_text(errors="replace")
        for number, line in enumerate(content.splitlines(), 1):
            for match in TMP_REF.finditer(line):
                ref = match.group().removeprefix("/home/math/quadra/").rstrip(".")
                entry = by_path.get(ref)
                candidate = ROOT / ref
                if ref in tracked:
                    resolution = "tracked"
                elif entry:
                    resolution = entry["category"]
                elif candidate.is_dir():
                    resolution = "directory_or_import_search_path"
                elif candidate.exists():
                    resolution = "exists_not_in_snapshot_inventory"
                else:
                    resolution = "missing_or_dynamic_reference_REVIEW"
                hint = "output_or_example" if any(s in line for s in
                    ("--output", "write_text", "write_bytes", "mkdir", "--log", " > ")) else "READ_OR_OUTPUT_REVIEW"
                out.append({"source": name, "line": number, "reference": ref,
                            "resolution": resolution, "role_hint": hint,
                            "line_excerpt": line[:400]})
    return out


def restore(snapshot, name):
    manifest = json.loads((snapshot / "manifest.json").read_text())
    item = next((r for r in manifest["files"] if r["path"] == name and "archived_path" in r), None)
    if item is None:
        raise SystemExit("Requested path is not a preserved payload in this snapshot")
    dest = ROOT / name
    if not dest.resolve().is_relative_to(ROOT):
        raise SystemExit("Restore target outside repository")
    data = (snapshot / item["archived_path"]).read_bytes()
    if hashlib.sha256(data).hexdigest() != item["sha256"]:
        raise SystemExit("Archive hash mismatch")
    if dest.exists():
        if dest.is_file() and hashlib.sha256(dest.read_bytes()).hexdigest() == item["sha256"]:
            print("Already present with matching content:", name)
            return
        raise SystemExit("Refusing to overwrite an existing different working file")
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(snapshot / item["archived_path"], dest)
    print("Restored missing research dependency:", name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", action="store_true")
    parser.add_argument("--snapshot", help="existing snapshot for --restore-path")
    parser.add_argument("--restore-path")
    parser.add_argument("--verify", action="store_true", help="verify every saved payload against its hash")
    args = parser.parse_args()
    if args.verify:
        if not args.snapshot:
            parser.error("--verify requires --snapshot")
        snapshot = (ROOT / args.snapshot).resolve()
        manifest = json.loads((snapshot / "manifest.json").read_text())
        count = 0
        for item in manifest["files"]:
            if "archived_path" not in item:
                continue
            data = (snapshot / item["archived_path"]).read_bytes()
            if len(data) != item["size"] or hashlib.sha256(data).hexdigest() != item["sha256"]:
                raise SystemExit("Missing/corrupt archive payload: " + item["path"])
            count += 1
        print(json.dumps({"verified_payloads": count, "snapshot": args.snapshot}))
        return
    if args.restore_path:
        if not args.snapshot:
            parser.error("--restore-path requires --snapshot")
        restore((ROOT / args.snapshot).resolve(), args.restore_path)
        return
    records = inventory()
    counts = Counter(r["category"] for r in records)
    sizes = Counter()
    for r in records:
        sizes[r["category"]] += r["size"]
    summary = {"counts": dict(counts), "bytes": dict(sizes),
               "credential_candidates": [{"path": r["path"], "reason": r["reason"]}
                   for r in records if r["category"] == "credential_candidate"],
               "oversized_research": [{"path": r["path"], "bytes": r["size"]}
                   for r in records if r["category"] == "oversized_research"]}
    if args.save:
        if summary["credential_candidates"] or summary["oversized_research"]:
            raise SystemExit("Review credential/oversized candidates before saving: " + json.dumps(summary))
        stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d/%H%M%SZ")
        target = ROOT / "research_archive" / stamp
        target.mkdir(parents=True, exist_ok=False)
        for r in records:
            if r["category"] != "research":
                continue
            origin = ROOT / r["path"]
            destination = target / "files" / r["path"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            # Re-read at the copying boundary: actively running agents may write.
            before = origin.stat()
            shutil.copy2(origin, destination)
            data = destination.read_bytes()
            r.update(sha256=hashlib.sha256(data).hexdigest(), size=len(data),
                     archived_path=destination.relative_to(target).as_posix())
            if before.st_mtime_ns != origin.stat().st_mtime_ns:
                r["concurrent_write_warning"] = "source changed during copy; resnapshot at next checkpoint"
            if any(pattern.search(data) for pattern in SECRETS.values()):
                raise SystemExit("Secret-pattern match at copy boundary: " + r["path"])
        deps = dependencies(records)
        manifest = {"created_utc": datetime.now(timezone.utc).isoformat(),
                    "git_head": git("rev-parse", "HEAD").decode().strip(),
                    "status": "PRESERVATION ONLY; no verification upgrade",
                    "summary": summary, "files": records}
        # Package names/versions only: never serialize environment variables,
        # credential-bearing index URLs, or pip configuration.
        packages = json.loads(subprocess.check_output(
            [sys.executable, "-m", "pip", "list", "--format=json"], cwd=ROOT))
        (target / "runtime.json").write_text(json.dumps({
            "python": sys.version, "packages": packages,
            "status": "observed environment metadata, not a tested lockfile"
        }, indent=2) + "\n")
        (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        (target / "ignored_references.json").write_text(json.dumps(deps, indent=2) + "\n")
        summary.update(snapshot=target.relative_to(ROOT).as_posix(),
                       ignored_reference_count=len(deps),
                       reference_resolutions=dict(Counter(d["resolution"] for d in deps)))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
