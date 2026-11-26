"""Create or verify the pre-fit source manifest and receipt, from repository root."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

study = Path("anchoredRI/tools/subsample")
workflow = Path(".github/workflows/anchoredri-task-h-pilot.yml")
manifest = study / "SOURCE.sha256"
launch = study / "launch-pilot.json"
approved = "e5a4577c1fe73869c6a24ce64c8759f190de3197"
package = ["anchoredRI/" + name for name in ("DESCRIPTION", "NAMESPACE", "R", "src", "inst", "man", "LICENSE")]
base_trees = study / "package-trees.txt"


def git(*args):
    return subprocess.check_output(["git", *args], text=True).strip()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def files():
    sources = [path for path in study.iterdir() if path.is_file() and path not in (manifest, launch)]
    sources.extend([workflow, study / "archived-shares" / "fits.csv"])
    sources.extend(Path(p) for p in git("ls-files", "--", *package).splitlines())
    return sorted(set(sources))


mode = sys.argv[1]
if mode == "create":
    lines = [git("rev-parse", f"{approved}:{p}") + " " + p for p in package]
    base_trees.write_text("\n".join(lines) + "\n")
elif mode != "verify":
    raise SystemExit("Expected create or verify")
for line in base_trees.read_text().splitlines():
    blob, path = line.split(" ", 1)
    if blob != git("rev-parse", "HEAD:" + path):
        raise SystemExit("Package tree mismatch " + path)

content = "".join(sha(p) + "  " + str(p) + "\n" for p in files())
if mode == "create":
    manifest.write_text(content)
else:
    receipt = json.loads(launch.read_text())
    parent = git("rev-parse", "HEAD^")
    if (git("diff", "--name-only", parent, "HEAD") != str(launch)
        or receipt["implementation_commit"] != parent or receipt["approved_protocol"] != approved
        or receipt["manifest_sha256"] != sha(manifest) or manifest.read_text() != content):
        raise SystemExit("Implementation freeze or launch receipt mismatch")
    print(f"Verified {len(files())} source hashes at {parent}")
with (study / "pilot-roster.csv").open() as stream:
    rows = list(csv.DictReader(stream))
if len({r["Index"] for r in rows}) != 20 or len({(r["Analyte"], r["pathology_stratum"]) for r in rows}) != 20:
    raise SystemExit("Expected one pilot per analyte and stratum")
