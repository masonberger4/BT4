# Security Policy

## Supported versions

BT4 is under active, single-trunk development. Security fixes are made against
the `main` branch and released in the latest release. Please make sure you're on
the latest release (or current `main`) before reporting an issue.

| Version           | Supported          |
| ----------------- | ------------------ |
| `main` / latest release | ✅            |
| older releases    | ❌                 |

## Reporting a vulnerability

Please report security vulnerabilities **privately** — do not open a public
issue, pull request, or discussion for a security problem.

The preferred channel is GitHub's **private security advisory**: go to the
repository's **Security** tab and choose **"Report a vulnerability"** to open a
private advisory that only the maintainers can see.

If you cannot use GitHub security advisories, contact the maintainers privately
at **[INSERT CONTACT METHOD]**.

> **Maintainers:** replace `[INSERT CONTACT METHOD]` above with a real private
> contact address before publishing this project.

When reporting, please include as much as you can:

- a description of the vulnerability and its impact,
- steps to reproduce (a minimal protein + config, or sequence, if relevant),
- the BT4 version (`bt4 --version`) and how you installed it, and
- any relevant logs or proof-of-concept.

We will acknowledge your report as soon as we can, keep you informed of
progress, and coordinate a disclosure timeline with you. Please give us a
reasonable opportunity to address the issue before any public disclosure.

## Scope and threat model

BT4 is a **local, offline tool**. The library, the CLI, and the BT4 Studio
desktop app run entirely on your machine and compute results locally — **nothing
leaves the machine** during normal use. There is no telemetry, and the app does
not phone home.

The one designed exception is the opt-in **ASSP splice cross-check**, which is a
network call to an external service. It is described in `CLAUDE.md` as
**roadmap** and is **not shipped today**; when it lands it will be strictly
opt-in, gated behind a flag and an optional extra, out of the optimization
loop, and clearly labeled as network-derived. Until then, no BT4 code contacts
the network as part of a run.

The optional `service/` HTTP API (also roadmap) is intended for automation and
would expose BT4 over the network; when it ships it is expected to be
authenticated and resource-bounded. If you deploy any future networked surface,
treat it as internet-facing and secure it accordingly.
