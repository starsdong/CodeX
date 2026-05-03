#!/usr/bin/env python3
"""Check generated HTML links."""

from __future__ import annotations

import argparse
import sys
import urllib.error
import urllib.request
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urldefrag, urlparse

ROOT = Path(__file__).resolve().parents[1]
DOCS_DIR = ROOT / "docs"


class LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        attrs_dict = dict(attrs)
        href = attrs_dict.get("href")
        if href:
            self.links.append(href)


def internal_target_exists(page: Path, href: str) -> bool:
    href, _fragment = urldefrag(href)
    if not href or href.startswith("#"):
        return True
    parsed = urlparse(href)
    if parsed.scheme or href.startswith("//"):
        return True
    target = (page.parent / href).resolve()
    try:
        target.relative_to(DOCS_DIR.resolve())
    except ValueError:
        return False
    return target.exists()


def external_ok(href: str) -> tuple[bool, str]:
    request = urllib.request.Request(href, method="HEAD", headers={"User-Agent": "SoftHadronMirror/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=20) as response:
            return response.status < 400, str(response.status)
    except urllib.error.HTTPError as exc:
        if exc.code == 405:
            try:
                get_request = urllib.request.Request(
                    href, headers={"User-Agent": "SoftHadronMirror/1.0"}
                )
                with urllib.request.urlopen(get_request, timeout=20) as response:
                    return response.status < 400, str(response.status)
            except Exception as nested_exc:  # noqa: BLE001
                return False, str(nested_exc)
        return False, str(exc)
    except Exception as exc:  # noqa: BLE001
        return False, str(exc)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--external", action="store_true", help="Also check external URLs over the network.")
    args = parser.parse_args()

    failures: list[str] = []
    external_links: set[str] = set()
    for page in sorted(DOCS_DIR.glob("*.html")):
        link_parser = LinkParser()
        link_parser.feed(page.read_text(encoding="utf-8"))
        for href in link_parser.links:
            parsed = urlparse(href)
            if parsed.scheme in {"http", "https"}:
                external_links.add(href)
                continue
            if not internal_target_exists(page, href):
                failures.append(f"{page.relative_to(ROOT)} -> missing {href}")

    if args.external:
        for href in sorted(external_links):
            ok, reason = external_ok(href)
            if not ok:
                failures.append(f"external failed {href}: {reason}")

    if failures:
        print("Link check failed:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    if args.external:
        print(f"Link check OK: internal links and {len(external_links)} external URLs checked.")
    else:
        print("Link check OK: internal links checked. Use --external to check remote URLs.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

