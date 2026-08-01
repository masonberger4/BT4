# Packaging BT4 Studio

This directory builds a **standalone BT4 Studio** desktop bundle (no Python
install required for end users) with [PyInstaller](https://pyinstaller.org/).

- `bt4-studio.spec` — the PyInstaller spec (one-folder bundle; a macOS `.app` is
  produced on macOS). It collects the Qt runtime (via PySide6's own hook),
  pyqtgraph (skipping its GUI `examples`), and BT4's packaged codon tables.
- `bt4_studio_entry.py` — the entry script PyInstaller analyzes; it just calls
  `bt4.app.main`.

## Build locally

```bash
pip install -e '.[app,packaging]'          # from the repo root
cd packaging
pyinstaller --clean --noconfirm bt4-studio.spec
```

The bundle lands in `packaging/dist/BT4 Studio/` (run the `BT4 Studio`
executable inside it), and on macOS also as `packaging/dist/BT4 Studio.app`.

On a headless Linux box you'll need the Qt runtime libraries to launch it, e.g.
`sudo apt-get install -y libegl1 libgl1 libglib2.0-0 libxkbcommon0 libdbus-1-3`.

## Releases

Pushing a version tag (`vX.Y.Z`) runs
[`.github/workflows/release.yml`](../.github/workflows/release.yml), which builds
this bundle for Linux / macOS / Windows plus the Python sdist + wheel, and
attaches them all to a GitHub Release. `build/` and `dist/` are git-ignored. The
publish step is idempotent — if the release already exists (e.g. it was created
in the UI, which also creates the triggering tag), the assets are uploaded to it
with `--clobber` instead of failing on "release already exists". It also fails
loudly rather than publishing an empty, asset-less release.

## Repairing a release

If a tagged release ends up with **no bundles attached** — for example the tag
was created in the UI before this pipeline existed, or an early publish step
failed — rebuild it from its own source and attach the assets, without moving the
tag or deleting the release:

1. Trigger the **Release** workflow with `workflow_dispatch`, setting the `ref`
   input to the existing tag (e.g. `v0.2.0`). The workflow checks out that exact
   tag, builds the bundles + wheel/sdist from the tagged source, and idempotently
   uploads them to the existing release. Because it builds the tagged source, the
   bundles genuinely match the version they are published under.

   ```bash
   # with the GitHub CLI, from a clone with the workflow on the default branch:
   gh workflow run release.yml -f ref=v0.2.0
   ```

> **Caveat:** GitHub only exposes `workflow_dispatch` once `release.yml` is
> present on the repository's **default branch**. If the default branch does not
> yet carry the workflows, the manual dispatch is unavailable; land the workflows
> on the default branch first (a tag push always works regardless).
