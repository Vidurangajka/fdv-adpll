"""Print docs/report/report.html to docs/fdv-adpll-report.pdf.

    python scripts/build_report.py [--browser PATH]

Uses a headless Chromium (Chrome or Edge, found on the usual install paths or
given with --browser), so the PDF is exactly what the page's print stylesheet
lays out.  The page pulls its figures from results/ and silicon/layout/, so
rebuild those first (scripts/reproduce.py does both).  Web fonts come from
Google Fonts; offline, the page falls back to local serif / sans fonts.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "docs" / "report" / "report.html"
OUT = ROOT / "docs" / "fdv-adpll-report.pdf"

CANDIDATES = [
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "google-chrome", "google-chrome-stable", "chromium", "chromium-browser",
    "microsoft-edge",
]


def find_browser() -> str:
    for c in CANDIDATES:
        path = c if os.path.isabs(c) else shutil.which(c)
        if path and os.path.exists(path):
            return path
    sys.exit("no Chrome / Edge / Chromium found; pass --browser PATH")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--browser", help="Chromium-family executable")
    ap.add_argument("-o", type=pathlib.Path, default=OUT)
    args = ap.parse_args()

    browser = args.browser or find_browser()
    # a throwaway profile, so a browser the user already has open is not reused
    with tempfile.TemporaryDirectory() as profile:
        subprocess.run([
            browser, "--headless=new", "--disable-gpu", "--no-first-run",
            f"--user-data-dir={profile}",
            "--no-pdf-header-footer",
            "--virtual-time-budget=15000",      # let the web fonts arrive
            f"--print-to-pdf={args.o}",
            SRC.as_uri(),
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not args.o.exists() or args.o.stat().st_size < 10_000:
        sys.exit(f"the browser ran but {args.o} is missing or empty")
    print(f"wrote {args.o} ({args.o.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
