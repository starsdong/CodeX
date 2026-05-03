#!/usr/bin/env python3
"""Build the local static mirror and Google Sites update packet."""

from __future__ import annotations

import hashlib
import html
import json
import shutil
import subprocess
import sys
import tomllib
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTENT_DIR = ROOT / "content" / "pages"
DOCS_DIR = ROOT / "docs"
PUBLISH_DIR = ROOT / "publish"
TEMPLATE_DIR = ROOT / "templates"


def read_site_config() -> dict:
    with (ROOT / "site.toml").open("rb") as fh:
        return tomllib.load(fh)


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def markdown_to_html(markdown: str) -> str:
    try:
        result = subprocess.run(
            ["pandoc", "-f", "gfm+raw_html", "-t", "html"],
            input=markdown,
            text=True,
            capture_output=True,
            check=True,
        )
        return result.stdout
    except FileNotFoundError:
        sys.exit("pandoc is required to build this mirror.")
    except subprocess.CalledProcessError as exc:
        sys.stderr.write(exc.stderr)
        sys.exit(exc.returncode)


def page_title(markdown: str, fallback: str) -> str:
    for line in markdown.splitlines():
        if line.startswith("# "):
            return line[2:].strip()
    return fallback


def render_nav(nav_items: list[dict], current_slug: str) -> str:
    links = []
    for item in nav_items:
        attrs = ""
        if item["slug"] == current_slug:
            attrs = ' aria-current="page"'
        links.append(
            f'<a href="{html.escape(item["slug"])}.html"{attrs}>'
            f'{html.escape(item["title"])}</a>'
        )
    return "\n      ".join(links)


def render_page(template: str, site: dict, item: dict, markdown: str) -> str:
    content_html = markdown_to_html(markdown)
    replacements = {
        "content_html": content_html,
        "description": site["description"],
        "nav_html": render_nav(site["nav"], item["slug"]),
        "page_title": page_title(markdown, item["title"]),
        "site_title": site["title"],
        "source_url": site["source_url"],
    }
    page = template
    for key, value in replacements.items():
        page = page.replace("{{" + key + "}}", value)
    return page


def read_previous_manifest() -> dict:
    path = DOCS_DIR / "manifest.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_assets() -> None:
    DOCS_DIR.mkdir(exist_ok=True)
    shutil.copyfile(TEMPLATE_DIR / "styles.css", DOCS_DIR / "styles.css")
    source_assets = ROOT / "assets"
    dest_assets = DOCS_DIR / "assets"
    if dest_assets.exists():
        shutil.rmtree(dest_assets)
    if source_assets.exists():
        shutil.copytree(source_assets, dest_assets)


def write_publish_packet(site: dict, pages: list[dict], changed_slugs: set[str]) -> None:
    PUBLISH_DIR.mkdir(exist_ok=True)
    lines = [
        "# Google Sites Update Packet",
        "",
        "Generated from the local Markdown source by `python3 scripts/build.py`.",
        "",
        "## Checklist",
        "",
    ]
    if changed_slugs:
        for page in pages:
            if page["slug"] in changed_slugs:
                url = site["published_base_url"] + page["google_path"]
                lines.append(f"- [ ] Update {page['title']}: {url}")
    else:
        lines.append("- [ ] No content changes detected since the previous build.")
    lines.extend(
        [
            "- [ ] Preview in Google Sites",
            "- [ ] Publish in Google Sites",
            "- [ ] Add a row to `SYNC_LOG.md` with the git commit and publish date",
            "",
            "## Google Sites Editor Notes",
            "",
            "- Use Google Sites headings/lists/links to recreate the Markdown structure.",
            "- For the Home calendar block, use Insert > Embed with the iframe URL from `content/pages/home.md`.",
            "- Keep raw Google Takeout exports out of git until reviewed.",
            "",
            "## Changed Page Content",
            "",
        ]
    )
    for page in pages:
        if page["slug"] not in changed_slugs:
            continue
        source_path = CONTENT_DIR / f"{page['slug']}.md"
        lines.extend(
            [
                f"### {page['title']}",
                "",
                f"- Source: `content/pages/{page['slug']}.md`",
                f"- Google Sites path: `{page['google_path']}`",
                "",
                "Copy-ready Markdown/source text:",
                "",
                "```markdown",
                source_path.read_text(encoding="utf-8").rstrip(),
                "```",
                "",
            ]
        )
    (PUBLISH_DIR / "google-sites-update.md").write_text(
        "\n".join(lines).rstrip() + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--all",
        action="store_true",
        help="Include every page in the Google Sites update packet.",
    )
    args = parser.parse_args()

    site = read_site_config()
    template = (TEMPLATE_DIR / "base.html").read_text(encoding="utf-8")
    previous = read_previous_manifest()
    previous_pages = previous.get("pages", {})
    manifest = {"pages": {}}
    rendered_pages = []

    DOCS_DIR.mkdir(exist_ok=True)
    for item in site["nav"]:
        source_path = CONTENT_DIR / f"{item['slug']}.md"
        markdown = source_path.read_text(encoding="utf-8")
        source_hash = sha256_text(markdown)
        html_page = render_page(template, site, item, markdown)
        output_path = DOCS_DIR / f"{item['slug']}.html"
        output_path.write_text(html_page, encoding="utf-8")
        rendered_pages.append(item)
        manifest["pages"][item["slug"]] = {
            "source": f"content/pages/{item['slug']}.md",
            "source_sha256": source_hash,
            "output": f"docs/{item['slug']}.html",
            "google_url": site["published_base_url"] + item["google_path"],
        }
        if item["slug"] == "home":
            (DOCS_DIR / "index.html").write_text(html_page, encoding="utf-8")

    if args.all:
        changed_slugs = set(manifest["pages"])
    else:
        changed_slugs = {
            slug
            for slug, entry in manifest["pages"].items()
            if previous_pages.get(slug, {}).get("source_sha256") != entry["source_sha256"]
        }
    write_assets()
    (DOCS_DIR / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_publish_packet(site, rendered_pages, changed_slugs)
    if args.all:
        print("Built static mirror. Update packet includes all pages.")
    elif changed_slugs:
        changed = ", ".join(sorted(changed_slugs))
        print(f"Built static mirror. Changed pages: {changed}")
    else:
        print("Built static mirror. No content changes detected.")


if __name__ == "__main__":
    main()
