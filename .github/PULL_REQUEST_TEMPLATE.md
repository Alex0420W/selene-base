## Summary

What this PR does, in 1–3 sentences.

## Motivation

Why this change is worth making (link the issue if there is one).

## Checklist

- [ ] `ruff check .` passes
- [ ] `ruff format --check .` passes
- [ ] `pytest` passes locally
- [ ] New behaviour has a test (or an explicit `# noqa` reason if not)
- [ ] Public-facing changes (CLI flags, criterion math, file layout) are
      reflected in the README or the relevant module docstring
- [ ] No real-data fixtures were committed; new tests run in CI without
      downloaded data (or are guarded with `pytest.mark.skipif`)

## Scope notes

- Things this PR explicitly does **not** address (so reviewers know what
  to expect in follow-ups).
