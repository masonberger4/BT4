# Installing BT4 Studio — a friendly step-by-step guide

Welcome! This guide walks you through downloading and opening **BT4 Studio**, the
desktop app for BT4. You do **not** need to know anything about coding or the
command line — you download one file and open it, just like any other app.

A few reassurances before we start:

- **Everything runs on your own computer.** BT4 Studio is 100% local and offline.
  Nothing you paste in is ever uploaded anywhere.
- **The first-time security warning is normal.** The app isn't code-signed, so the
  very first time you open it, your operating system shows a one-time "we don't
  recognize this app" message. This is expected and safe, and below we show you
  exactly, click by click, how to get past it.
- **You only do the warning step once.** After the first launch, the app opens
  normally like anything else.

> This is the same information the [README summarizes](../README.md#install-bt4-studio-no-coding-required)
> — just with more detail and more hand-holding.

---

## Step 1 — Download the app

1. Open the [**Releases page**](https://github.com/masonberger4/BT4/releases).
2. Click the **latest release** at the top.
3. Scroll down to the **Assets** section (a list of downloadable files).
4. Download the **one** file that matches your computer:

| Your computer | Download this file |
|---|---|
| **Windows** | `BT4-Studio-Windows.exe` |
| **Mac** | `BT4-Studio-macOS.dmg` |
| **Linux** | `BT4-Studio-Linux-x86_64` |

Pick only the file for the kind of computer you have. If you're not sure, you're
almost certainly on **Windows** or **Mac**.

> **Mac — Apple Silicon vs Intel:** the `.dmg` is built for **Apple Silicon** Macs
> (M1/M2/M3/M4 — 2020 or later). To check, click the **Apple menu → About This
> Mac**: if it says **Chip** you're on Apple Silicon (use the `.dmg`); if it says
> **Processor / Intel**, use the "from source" developer install in the
> [README](../README.md#install-for-developers-from-source) for now.

> **Don't see any of these files under Assets?** The download for that release
> hasn't been built yet. Check back shortly, or ask a maintainer — the
> [packaging guide](../packaging/README.md#repairing-a-release) explains how they
> produce them.

---

## Step 2 — Open the app (per computer)

Jump to the section for your computer.

### Windows

1. Find the downloaded `BT4-Studio-Windows.exe` (usually in your **Downloads**
   folder).
2. **Double-click it.**

That's it — there's nothing to install. The first time, you'll likely see a blue
box (see [Getting past the first-time warning](#getting-past-the-first-time-warning)
below) — that's normal.

### Mac

1. **Double-click** the downloaded `BT4-Studio-macOS.dmg`. A small window opens.
2. **Drag the `BT4 Studio` icon into the `Applications` folder** shown in that
   window.
3. Open **Applications** (in Finder, click **Go → Applications**), find
   **BT4 Studio**, and open it from there.

The first time you open it, macOS shows a warning (see
[Getting past the first-time warning](#getting-past-the-first-time-warning)) —
this is expected.

### Linux

1. Find the downloaded `BT4-Studio-Linux-x86_64`.
2. Linux won't run a downloaded file until you mark it as a program. In your file
   manager, **right-click the file → Properties → Permissions**, and turn on
   **"Allow executing as program"** (the exact wording varies by desktop).
3. **Double-click** the file to open it.

If your file manager doesn't have that option, see the
[Linux terminal fallback](#linux-it-wont-run) below.

---

## Getting past the first-time warning

Because the app isn't code-signed, the **first** time you open it your computer
warns you that it doesn't recognize the developer. **This is expected and the app
is safe.** Here's exactly how to open it anyway.

### Windows — "Windows protected your PC"

You'll see a blue box titled **"Windows protected your PC"**.

1. Click **More info** (a small link in that box — it may not be obvious at
   first).
2. A **Run anyway** button appears. Click **Run anyway**.

The app opens. You won't need to do this again.

### Mac — "cannot be opened because Apple cannot check it..."

You may see a message like *"'BT4 Studio' cannot be opened because Apple cannot
check it for malicious software."* To open it:

1. In **Applications**, **right-click** (or hold **Control** and click) the
   **BT4 Studio** app.
2. Choose **Open** from the menu.
3. In the dialog that appears, click **Open** again.

You only need this right-click-Open trick the **first** time. After that, a normal
double-click works.

On newer versions of macOS the wording differs slightly — instead, open
**System Settings → Privacy & Security**, scroll down, and click the **Open
Anyway** button that appears there after your first attempt to open the app.

### Mac — "is damaged and can't be opened. You should move it to the Trash"

Don't worry — the app almost certainly **isn't** damaged. macOS adds a hidden
"downloaded from the internet" quarantine tag to unsigned apps, and sometimes
shows this misleading message instead of the normal warning.

**First, try this Finder-only step:** make sure you actually **dragged BT4 Studio
into your Applications folder** (Step 2 above) and are opening it from there — not
running it from inside the `.dmg` window or the Downloads folder. Then try the
right-click **→ Open → Open** trick above.

**If it still says "damaged,"** here's the reliable fix. Open the **Terminal**
app (find it via **Applications → Utilities → Terminal**, or press **Cmd + Space**
and type `Terminal`), then copy and paste this line and press **Return**:

```
xattr -dr com.apple.quarantine "/Applications/BT4 Studio.app"
```

That command simply removes the hidden "downloaded from the internet" quarantine
flag from the app. You only need to run it once. Now open **BT4 Studio** normally
from Applications.

### Linux: it won't run

If double-clicking does nothing, the file probably isn't marked as a program yet.
The easiest fix is in your file manager: **right-click → Properties →
Permissions → Allow executing as program**, then double-click.

If that isn't available, open a **terminal** in the folder where you downloaded
the file and run these two lines (press **Enter** after each):

```bash
chmod +x "BT4-Studio-Linux-x86_64"
./"BT4-Studio-Linux-x86_64"
```

The first line marks the file as runnable; the second one starts it.

---

## "My antivirus flagged it"

Occasionally an antivirus program flags BT4 Studio as suspicious. **This is a
false alarm.** BT4 Studio is bundled into a single file with a common tool
(PyInstaller), and because *many* apps are packaged the exact same way, some
antivirus software flags the whole pattern as a precaution rather than because
anything is actually wrong.

The app is safe. You can tell your antivirus to **allow** or **whitelist** it (the
option is usually right there in the warning, or under your antivirus's
"quarantine" / "allowed apps" settings).

---

## It still won't open

If you're stuck, run through these:

- **Re-download the file.** Downloads sometimes finish incompletely, which can
  cause "damaged" or "won't open" errors. Downloading it fresh often fixes it.
- **Double-check you picked the right file** for your computer (Windows `.exe`,
  Mac `.dmg`, or the Linux file) from Step 1.
- **On a Mac,** make sure you **dragged BT4 Studio into Applications** and are
  opening it from there — not from inside the `.dmg` window.
- **As a last resort,** the app can also be installed from source (the "from
  source" / developer install in the [README](../README.md#install-for-developers-from-source)).
  It needs a bit of setup, but it always works.

---

## Using BT4 Studio

Once it opens, the basic flow is simple:

1. **Paste a protein** sequence (or upload a FASTA file).
2. **Pick a target organism.**
3. **Set your constraints** (GC target, max homopolymer, forbidden motifs, and
   more).
4. Click **Optimize**.

You'll get an optimized coding sequence, an honest optimality badge, a metrics
table, and an interactive trade-off frontier — all computed right on your machine.

---

## Uninstall

Removing BT4 Studio is just deleting the file you downloaded — there's no
installer to run.

- **Windows:** delete `BT4-Studio-Windows.exe`.
- **Mac:** open **Applications** and drag **BT4 Studio** to the **Trash**.
- **Linux:** delete the `BT4-Studio-Linux-x86_64` file.

---

## A note on the security warnings

The one-time warning appears because the app isn't **code-signed**. Removing it
would mean paying for certificates (an Apple Developer ID for Mac, an Authenticode
certificate for Windows) just to silence a prompt you click past once — so BT4
Studio is distributed unsigned on purpose. As shown above, getting past the warning
takes a couple of clicks, and the app is safe.
