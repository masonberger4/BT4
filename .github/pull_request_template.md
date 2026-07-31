<!--
Thanks for contributing to BT4! Please read CONTRIBUTING.md and CLAUDE.md (the
design constitution) first. Fill in the sections below and keep the PR focused.
-->

## Summary

<!-- What does this PR do, and why? One or two sentences. -->

## Changes

<!-- Bullet the notable changes. Note any new constraints/objectives, contract
     changes, or user-facing changes. -->

-
-

## Checklist

- [ ] `ruff check src tests` is clean
- [ ] `mypy` is clean
- [ ] `lint-imports` is clean (strict layering contract preserved)
- [ ] `pytest` is green
- [ ] The layering is preserved (no new cross-layer or private-symbol imports)
- [ ] If this adds a `Constraint` or `ObjectiveTerm`: it's a new file + a public
      export, **and** its property test is added
      (`ok_suffix` ⇔ `validate` / sufficient `context_len`, or `delta` == `score`)
- [ ] Docs updated if the architecture or contracts changed
      (`CLAUDE.md`, and `README.md` if user-facing)
- [ ] No roadmap-only feature is described or presented as if it already ships

## Notes for reviewers

<!-- Anything reviewers should pay special attention to, open questions, or
     follow-ups. Optional. -->
