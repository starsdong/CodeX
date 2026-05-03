#!/usr/bin/env python3
"""Compare local Markdown text with the published Google Site or saved fixtures."""

from __future__ import annotations

import argparse
import difflib
import re
import sys
import tomllib
import urllib.request
from html import unescape
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class GoogleSitesTextParser(HTMLParser):
    block_tags = {"h1", "h2", "h3", "h4", "p", "li", "small"}

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.in_main = False
        self.skip_depth = 0
        self.current_tag: str | None = None
        self.current_text: list[str] = []
        self.blocks: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attrs_dict = dict(attrs)
        if tag in {"script", "style", "svg"}:
            self.skip_depth += 1
            return
        if tag == "div" and attrs_dict.get("jsname") == "ZBtY8b":
            self.in_main = True
        if not self.in_main or self.skip_depth:
            return
        if tag in self.block_tags and self.current_tag is None:
            self.current_tag = tag
            self.current_text = []

    def handle_endtag(self, tag: str) -> None:
        if tag in {"script", "style", "svg"} and self.skip_depth:
            self.skip_depth -= 1
            return
        if not self.in_main or self.skip_depth:
            return
        if tag == self.current_tag:
            text = normalize_line("".join(self.current_text))
            if text and text not in {"Report abuse", "Page details", "Page updated"}:
                self.blocks.append(text)
            self.current_tag = None
            self.current_text = []

    def handle_data(self, data: str) -> None:
        if self.in_main and not self.skip_depth and self.current_tag is not None:
            self.current_text.append(data)


def normalize_line(value: str) -> str:
    value = unescape(value).replace("\xa0", " ")
    value = re.sub(r"^\s*[-*+]\s+", "", value)
    return re.sub(r"\s+", " ", value).strip()


def normalize_text(value: str) -> str:
    lines = [normalize_line(line) for line in value.splitlines()]
    lines = [line for line in lines if line]
    return "\n".join(lines)


def markdown_to_plain(markdown: str) -> str:
    markdown = re.sub(r"^!\[[^\]]*\]\([^)]+\)\s*$", "", markdown, flags=re.MULTILINE)
    markdown = re.sub(r"<div class=\"calendar-embed\">.*?</div>", "", markdown, flags=re.DOTALL)
    text = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", markdown)
    text = re.sub(r"^#{1,6}\s*", "", text, flags=re.MULTILINE)
    text = re.sub(r"<[^>]+>", "", text)
    text = re.sub(r"[*_`]", "", text)
    return normalize_text(text)


def remote_to_plain(html: str) -> str:
    parser = GoogleSitesTextParser()
    parser.feed(html)
    return normalize_text("\n".join(parser.blocks))


def fetch_remote(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "SoftHadronMirror/1.0"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8", errors="replace")


def read_site_config() -> dict:
    with (ROOT / "site.toml").open("rb") as fh:
        return tomllib.load(fh)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fixture-dir",
        type=Path,
        help="Directory with softhadron_<slug>.html files, used instead of network fetches.",
    )
    parser.add_argument("--show-diff", action="store_true")
    args = parser.parse_args()

    site = read_site_config()
    failed = False
    for item in site["nav"]:
        slug = item["slug"]
        local_md = (ROOT / "content" / "pages" / f"{slug}.md").read_text(encoding="utf-8")
        local_text = markdown_to_plain(local_md)
        if args.fixture_dir:
            html = (args.fixture_dir / f"softhadron_{slug}.html").read_text(
                encoding="utf-8", errors="replace"
            )
        else:
            html = fetch_remote(site["published_base_url"] + item["google_path"])
        remote_text = remote_to_plain(html)
        ratio = difflib.SequenceMatcher(None, local_text, remote_text).ratio()
        ok = ratio >= 0.92
        status = "OK" if ok else "DRIFT"
        print(f"{status} {slug}: similarity {ratio:.3f}")
        if not ok:
            failed = True
            if args.show_diff:
                diff = difflib.unified_diff(
                    remote_text.splitlines(),
                    local_text.splitlines(),
                    fromfile=f"remote/{slug}",
                    tofile=f"local/{slug}",
                    lineterm="",
                )
                print("\n".join(diff))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
