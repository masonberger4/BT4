# PyInstaller spec for BT4 Studio (the PySide6 desktop app).
#
# The goal is a *double-clickable* app for people who never touch a terminal, so
# the artifact shape is chosen per OS:
#
#   * Windows -> a single-file "BT4 Studio.exe" (one download, double-click).
#   * Linux   -> a single-file "BT4 Studio" executable (download, mark
#                executable, double-click). CI wraps this per release.
#   * macOS   -> a one-folder ".app" bundle (the standard Mac shape; it shows as
#                one icon). CI wraps it in a drag-to-Applications ".dmg".
#
# Build locally:
#     pip install -e '.[app,packaging]'
#     pyinstaller packaging/bt4-studio.spec
#
# PySide6's own PyInstaller hook collects the Qt runtime and platform plugins; we
# additionally pull in pyqtgraph and the packaged codon tables (loaded via
# importlib.resources at runtime).

import sys

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

APP_NAME = "BT4 Studio"
_is_mac = sys.platform == "darwin"
# Windows and Linux ship a single self-contained file; macOS ships a .app folder.
_onefile = not _is_mac

datas = []
binaries = []
hiddenimports = ["bt4.app.studio", "bt4.app.worker", "bt4.app.theme"]

# pyqtgraph does dynamic imports the static analyzer misses, so collect its
# submodules explicitly -- but SKIP `pyqtgraph.examples` (and `tests`), which
# import a live GUI and abort a headless PyInstaller analysis.
def _keep_pg(name: str) -> bool:
    return "examples" not in name and "tests" not in name


hiddenimports += collect_submodules("pyqtgraph", filter=_keep_pg)
datas += collect_data_files("pyqtgraph", excludes=["**/examples/**", "**/tests/**"])

# BT4's provenanced codon tables (TSV + provenance sidecars) and typing marker.
datas += collect_data_files(
    "bt4", includes=["**/*.tsv", "**/*.provenance.json", "py.typed"]
)

a = Analysis(
    ["bt4_studio_entry.py"],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["tkinter", "matplotlib", "torch", "numpy.distutils"],
    noarchive=False,
)

pyz = PYZ(a.pure)

if _onefile:
    # A single self-contained executable: everything (Qt, pyqtgraph, tables) is
    # embedded and unpacked to a temp dir at launch. One file to download and run.
    exe = EXE(
        pyz,
        a.scripts,
        a.binaries,
        a.datas,
        [],
        name=APP_NAME,
        debug=False,
        bootloader_ignore_signals=False,
        strip=False,
        upx=False,
        console=False,
        disable_windowed_traceback=False,
    )
else:
    # macOS: one-folder build wrapped in a .app bundle (shown as a single icon).
    exe = EXE(
        pyz,
        a.scripts,
        [],
        exclude_binaries=True,
        name=APP_NAME,
        debug=False,
        bootloader_ignore_signals=False,
        strip=False,
        upx=False,
        console=False,
        disable_windowed_traceback=False,
    )
    coll = COLLECT(
        exe,
        a.binaries,
        a.datas,
        strip=False,
        upx=False,
        name=APP_NAME,
    )
    app = BUNDLE(
        coll,
        name=f"{APP_NAME}.app",
        icon=None,
        bundle_identifier="com.bt4.studio",
    )
