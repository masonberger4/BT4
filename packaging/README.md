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
  `examples`), and BT4's packaged codon tables.
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

## Code signing (so the warnings go away)

Unsigned apps make first-time users see a one-time OS warning (Windows SmartScreen
"Windows protected your PC"; macOS Gatekeeper "cannot verify the developer" /
"damaged"). The [README](../README.md#install-bt4-studio-no-coding-required) and
[`docs/INSTALL.md`](../docs/INSTALL.md) tell users how to click past it.

**The signing pipeline is already wired into the `bundle` job of
[`release.yml`](../.github/workflows/release.yml).** It stays a no-op (a green,
unsigned build) until you add the certificate **secrets** below; the moment they
exist, each release is signed and — on macOS — notarized and stapled automatically,
and the warnings go away. The certificates cost money and are bound to your
identity, so they can only live as encrypted GitHub Actions secrets, never in the
repo.

### macOS — Developer ID + notarization

Requires an [Apple Developer Program](https://developer.apple.com/programs/)
membership ($99/yr). The macOS steps activate only when **both** the certificate
and the notary key are present. Add these repository secrets
(**Settings → Secrets and variables → Actions**):

| Secret | What it is / how to get it |
|---|---|
| `MACOS_CERTIFICATE_P12` | Your **Developer ID Application** cert exported from Keychain as a `.p12`, then base64-encoded: `base64 -i cert.p12 \| pbcopy`. |
| `MACOS_CERTIFICATE_PASSWORD` | The password you set when exporting the `.p12`. |
| `MACOS_SIGN_IDENTITY` | The identity string, e.g. `Developer ID Application: Your Name (TEAMID)` (see `security find-identity -v -p codesigning`). |
| `MACOS_NOTARY_API_KEY_P8` | An **App Store Connect API key** (`.p8`, Role: Developer) from App Store Connect → Users and Access → Integrations → Keys, base64-encoded. |
| `MACOS_NOTARY_KEY_ID` | The API key's Key ID. |
| `MACOS_NOTARY_ISSUER_ID` | The Issuer ID shown on the same Keys page. |

The workflow signs the `.app` with the hardened runtime, signs the `.dmg`,
submits it with `xcrun notarytool submit --wait`, and staples the ticket so
Gatekeeper opens it with no prompt, even offline.

### Windows — Authenticode

Requires a **code-signing certificate** from a CA (Sectigo, DigiCert, …). An
**OV** cert removes "unknown publisher" but SmartScreen reputation still builds
over downloads; an **EV** cert clears SmartScreen immediately. Add:

| Secret | What it is |
|---|---|
| `WINDOWS_CERTIFICATE_PFX` | The signing cert as a `.pfx`/`.p12`, base64-encoded: `base64 -w0 cert.pfx`. |
| `WINDOWS_CERTIFICATE_PASSWORD` | The `.pfx` password. |

The workflow signs `BT4-Studio-Windows.exe` with `signtool` (SHA-256, RFC-3161
timestamp) so it verifies with a trusted publisher.

> **Note:** to sign an **already-published** release (e.g. re-sign `v0.3.1` after
> adding the secrets), re-run its build — see [Repairing a release](#repairing-a-release).
> The signing steps have no local prerequisites; a plain unsigned local
> `pyinstaller` build (above) is unaffected.

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
