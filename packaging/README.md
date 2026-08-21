# Packaging BT4 Studio

This directory builds a **standalone, double-clickable BT4 Studio** app (no
Python install required for end users) with
[PyInstaller](https://pyinstaller.org/). The artifact shape is chosen per OS so a
non-technical user downloads one file and opens it:

| OS | `pyinstaller` output | Released asset |
|---|---|---|
| Windows | one-file `BT4 Studio.exe` | `BT4-Studio-Windows.exe` |
| macOS | one-folder `BT4 Studio.app` | `BT4-Studio-macOS.dmg` (drag-to-Applications) |
| Linux | one-file `BT4 Studio` executable | `BT4-Studio-Linux-x86_64` |

- `bt4-studio.spec` — the PyInstaller spec. It emits a **single file** on Windows
  and Linux and a **`.app` bundle** on macOS (`_onefile = not _is_mac`). It
  collects the Qt runtime (via PySide6's own hook), pyqtgraph (skipping its GUI
  `examples`), and every non-Python file BT4 ships (the codon / tRNA / enzyme
  tables, their provenance sidecars, the weight-hash pins and the committed splice
  attestations). Those patterns are pinned by `tests/test_bundle_spec.py`, because
  an enumeration of today's data files silently falls behind the tree -- and the
  resulting bundle fails only once it is launched.
- `bt4_studio_entry.py` — the entry script PyInstaller analyzes; it just calls
  `bt4.app.main`.

## Build locally

```bash
pip install -e '.[app,packaging]'          # from the repo root
cd packaging
pyinstaller --clean --noconfirm bt4-studio.spec
```

The result lands in `packaging/dist/`: a single `BT4 Studio` / `BT4 Studio.exe`
file on Linux/Windows, and a `BT4 Studio.app` folder on macOS. Just run it.

On a headless Linux box you'll need the Qt runtime libraries to launch it, e.g.
`sudo apt-get install -y libegl1 libgl1 libglib2.0-0 libxkbcommon0 libdbus-1-3`.
The CI release job installs exactly these before building the Linux one-file app.

**That set is enough to build and to run `--self-test` under `QT_QPA_PLATFORM=offscreen`
-- it is NOT enough to display the app on a real X11 desktop.** The offscreen platform
plugin needs almost none of X; `libqxcb.so` needs a dozen more `libxcb-*` libraries, and
without them the app aborts at startup with *"Could not load the Qt platform plugin
xcb"*. A normal desktop install already has them, so the gap only bites minimal,
server, WSL and container installs -- which is exactly where a maintainer tends to test.
The full list a user needs is in [`docs/INSTALL.md`](../docs/INSTALL.md#linux-it-wont-run);
it was read off the shipped binary's `NEEDED` entries and verified by launching the
packaged app on a bare system.

## Code signing (intentionally skipped)

BT4 Studio's releases are **not code-signed**, by choice. Signing would remove the
one-time OS warning first-time users see (Windows SmartScreen "Windows protected
your PC"; macOS Gatekeeper "cannot verify the developer" / "damaged"), but it needs
a paid Apple Developer ID membership and a Windows code-signing certificate — an
ongoing cost for a warning users click past once. The
[README](../README.md#install-bt4-studio-no-coding-required) and
[`docs/INSTALL.md`](../docs/INSTALL.md) walk users through that one-time step.

If you ever decide to sign the apps, add a signing step to the `bundle` job of
[`release.yml`](../.github/workflows/release.yml), keyed on GitHub Actions secrets
holding the certificates: macOS wants a Developer ID `codesign` (hardened runtime)
plus `xcrun notarytool submit --wait` and `xcrun stapler staple` on the `.dmg`;
Windows wants `signtool sign` (SHA-256 + RFC-3161 timestamp) on the `.exe`.

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
