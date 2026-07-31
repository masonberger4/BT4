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
attaches them all to a GitHub Release. `build/` and `dist/` are git-ignored.
