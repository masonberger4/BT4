"""Light/dark-aware styling for BT4 Studio.

This module is intentionally free of Qt imports: it only builds Qt Style Sheet
(QSS) strings from a small neutral palette. ``studio.py`` decides light vs dark
from the running application's palette and asks this module for the matching
sheet. The certificate badge gets its own colour helper so the solve quality is
shown honestly (green proven-optimal, amber truncated/heuristic, red relaxed).
"""

from __future__ import annotations

__all__ = ["badge_category", "badge_qss", "stylesheet"]

# Neutral palettes. No branding, just a legible surface in each theme.
_LIGHT = {
    "bg": "#f4f5f7",
    "panel": "#ffffff",
    "text": "#1b1f24",
    "muted": "#5b6169",
    "border": "#ccd1d8",
    "base": "#ffffff",
    "accent": "#2d6cdf",
    "accent_text": "#ffffff",
    "hover": "#3b78e6",
    "disabled": "#9aa2ac",
}
_DARK = {
    "bg": "#1c1f26",
    "panel": "#242832",
    "text": "#e6e9ed",
    "muted": "#9aa2ac",
    "border": "#39414d",
    "base": "#2b303b",
    "accent": "#4a83e0",
    "accent_text": "#ffffff",
    "hover": "#5a90e6",
    "disabled": "#5b6169",
}

# QSS template. Uses %(name)s substitution (there are no literal percent signs in
# the sheet), so palette values drop in without brace-escaping an f-string.
_QSS_TEMPLATE = """
QWidget { background-color: %(bg)s; color: %(text)s; }
QGroupBox {
    background-color: %(panel)s;
    border: 1px solid %(border)s;
    border-radius: 8px;
    margin-top: 14px;
    padding: 8px;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 4px;
    color: %(muted)s;
}
QLabel { background-color: transparent; }
QPlainTextEdit, QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QTableWidget {
    background-color: %(base)s;
    color: %(text)s;
    border: 1px solid %(border)s;
    border-radius: 6px;
    padding: 4px;
    selection-background-color: %(accent)s;
    selection-color: %(accent_text)s;
}
QHeaderView::section {
    background-color: %(panel)s;
    color: %(muted)s;
    border: none;
    border-bottom: 1px solid %(border)s;
    padding: 4px;
}
QTableWidget { gridline-color: %(border)s; }
QPushButton {
    background-color: %(accent)s;
    color: %(accent_text)s;
    border: none;
    border-radius: 6px;
    padding: 7px 14px;
    font-weight: 600;
}
QPushButton:hover { background-color: %(hover)s; }
QPushButton:disabled { background-color: %(disabled)s; color: %(panel)s; }
QStatusBar { color: %(muted)s; }
"""


def stylesheet(dark: bool) -> str:
    """Return the application QSS for a light or dark theme.

    Args:
        dark: True for the dark palette, False for the light palette.

    Returns:
        A QSS string suitable for ``QApplication.setStyleSheet``.
    """
    palette = _DARK if dark else _LIGHT
    return _QSS_TEMPLATE % palette


def badge_category(status_value: str) -> str:
    """Map an optimality-status value to a badge colour category.

    Args:
        status_value: ``OptimalityStatus.value`` such as ``"proven_optimal"``.

    Returns:
        One of ``"green"``, ``"amber"``, ``"red"``, or ``"neutral"``.
    """
    if status_value == "proven_optimal":
        return "green"
    if status_value in {"beam_truncated", "heuristic"}:
        return "amber"
    if status_value in {"relaxed", "context_capped", "gap_bounded"}:
        return "red"
    if status_value == "sampled":
        # A sampler makes no optimality claim at all, so it gets the neutral
        # colour rather than a warning shade: there is nothing degraded here to
        # warn about, and colouring "no claim" as "worse" would be its own lie.
        return "neutral"
    return "neutral"


def badge_qss(status_value: str) -> str:
    """Return the QSS for the certificate badge given an optimality status.

    Args:
        status_value: ``OptimalityStatus.value`` (empty for "no result yet").

    Returns:
        A QSS string targeting ``QLabel#certBadge``.
    """
    colours = {
        "green": "#1e8b4d",
        "amber": "#b9770e",
        "red": "#c0392b",
        "neutral": "#5b6169",
    }
    background = colours[badge_category(status_value)]
    return (
        "QLabel#certBadge {"
        f" background-color: {background};"
        " color: #ffffff;"
        " border-radius: 8px;"
        " padding: 10px 14px;"
        " font-weight: 700;"
        " }"
    )
