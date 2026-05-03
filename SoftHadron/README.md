# RNC Soft Physics Group Website Mirror

This folder is the git-first source for the public Google Site:

<https://sites.google.com/lbl.gov/rnc-softhadron/home>

Google Sites remains the official public site. Git is the editable source of truth, and updates back to Google Sites are done with a guided manual publish packet because the current Google Sites API does not support editing rebuilt/new Google Sites.

## Edit

1. Edit Markdown files in `content/pages/`.
2. Put reviewed public images/downloads in `assets/`.
3. Run:

```sh
python3 scripts/build.py
```

## Preview

```sh
python3 -m http.server -d docs 8000
```

Then open <http://localhost:8000>.

## Publish To Google Sites

1. Run `python3 scripts/build.py`.
2. Open `publish/google-sites-update.md`.
3. For each checked page, copy the updated content into the Google Sites editor.
4. Preview and publish in Google Sites.
5. Commit the local source and generated output to git.
6. Record the git commit and Google Sites publish date in `SYNC_LOG.md`.

For the first git-backed publication, generate a full packet:

```sh
python3 scripts/build.py --all
```

## Checks

```sh
python3 scripts/check_links.py
python3 scripts/check_remote.py
```

`check_remote.py` fetches the live Google Site and compares visible text with the local Markdown. If network access is unavailable, compare against saved HTML snapshots:

```sh
python3 scripts/check_remote.py --fixture-dir /private/tmp
```

External link checking is optional because it requires network access and can be slower:

```sh
python3 scripts/check_links.py --external
```

## Takeout Backups

Raw Google Takeout exports should go in `raw_takeout/`, which is ignored by git. Review any files before moving public content into tracked source files.
