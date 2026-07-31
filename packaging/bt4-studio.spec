# PyInstaller spec for BT4 Studio (the PySide6 desktop app).
#
# Build a one-folder bundle:
#     pip install -e '.[app]' pyinstaller
#     pyinstaller packaging/bt4-studio.spec
#
# Output: dist/BT4 Studio/  (a self-contained folder; run the "BT4 Studio"
# executable inside it). PySide6's own PyInstaller hook collects the Qt runtime
# and platform plugins; we additionally pull in pyqtgraph and the packaged codon
# tables (loaded via importlib.resources at runtime).

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

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

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="BT4 Studio",
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
    name="BT4 Studio",
)

app = BUNDLE(
    coll,
    name="BT4 Studio.app",
    icon=None,
    bundle_identifier="com.bt4.studio",
)
