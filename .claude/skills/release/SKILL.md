---
name: release
description: Run the cnr-ibba/nf-blobpurge release procedure (dev -> master) — version bump, validation, docs/ro-crate, PR/merge, git tag, and gh release create.
argument-hint: "<version> (e.g. 1.0.0)"
allowed-tools:
  - Read
  - Edit
  - Bash(git status*)
  - Bash(git diff*)
  - Bash(git log*)
  - Bash(git branch*)
  - Bash(nf-core pipelines *)
  - Bash(nf-test *)
  - Bash(gh pr list *)
  - Bash(gh pr view *)
  - Bash(gh pr checks *)
  - Bash(gh issue view *)
  - Bash(gh issue list *)
  - Bash(gh release list *)
disable-model-invocation: true
---

## Arguments

`$ARGUMENTS` is the release version, e.g. `1.0.0`. Required.

## Implementation note on confirmation gating

`gh pr create`, `gh pr merge`, and `gh release create` are deliberately
**not** in `allowed-tools` above, so they always trigger a permission
prompt — those are the genuinely irreversible/public actions.

`git tag`, `git push`, and `git commit`, however, are commonly covered by a
blanket `Bash(git *)` grant in `.claude/settings.local.json` — as they are
in this repo's own checkout. That file is gitignored and machine-specific,
so its exact contents aren't visible in this skill's diff and won't
necessarily match on another clone. Tool permissions alone will **not**
reliably gate these commands either way. This skill's own steps must
therefore pause and ask for explicit confirmation before running any `git
push` or `git tag`, regardless of what `settings.local.json` happens to
allow locally — the explicit-confirmation steps below are what actually
protect these actions, not the permission config. Treat this as a known
limitation, not an oversight.

## Resume logic

This procedure spans real calendar time (PR review, CI, merge can take
hours or days), so don't assume you're starting at step 1. At the start of
every invocation, inspect actual repo state and figure out what's already
done:

- Current branch (`git branch --show-current`).
- Whether `manifest.version` in `nextflow.config` still ends in `dev`.
- Whether tag `<version>` already exists (`git tag -l <version>`, also
  check `origin`).
- Whether a `dev -> master` PR exists and its state:
  `gh pr list --base master --head dev --state all --json number,state,url`.
- Whether a GitHub Release for `<version>` already exists:
  `gh release list --limit 5`.

Report a short status summary, then proceed from the next undone step.

## Phase 1 — prep on `dev`

1. Replace the `[date]` placeholder in the `## v<version>dev - [date]`
   heading of `CHANGELOG.md` with today's date, and finalize its entries.
2. Bump `manifest.version` in `nextflow.config` from `<version>dev` to
   `<version>`.
3. Align `template.version` in `.nf-core.yml` to `<version>`.

Show the diff of these edits before committing anything.

## Phase 2 — validate (safe to run automatically; local and reversible)

```bash
nf-core pipelines schema build
nf-core pipelines lint --dir .
nf-test test --tag pipeline --profile +docker --verbose
```

All three must pass clean before moving on.

## Phase 3 — docs & ro-crate

- Check whether `docs/usage.md`, `docs/output.md`, or `CITATIONS.md` need
  updates for anything added since the last release. Draft changes for
  review rather than auto-applying prose.
- Regenerate the RO-Crate metadata:
  ```bash
  nf-core pipelines rocrate .
  ```

## Phase 4 — Zenodo checkpoint (manual decision, no automation)

Surface this reminder before moving to Phase 7 (publish); do not act on it
automatically:

> No Zenodo <-> GitHub integration is configured yet — `README.md` still has
> `zenodo.XXXXXXX` placeholders and `manifest.doi` is empty in
> `nextflow.config`. If a DOI is wanted for this release, enable the Zenodo
> GitHub integration **now**, before `gh release create` — it only mints a
> DOI for releases created after the webhook is enabled. Otherwise, proceed
> and leave the placeholders as-is.

## Phase 5 — commit & PR

1. **Confirm explicitly** before `git push` of the `dev` commits from
   Phases 1–3.
2. Find the aggregator issue from the `prepare-release-issue` skill:
   ```bash
   gh issue list --search "in:title Release prep"
   ```
   Pull its `Closes #N` block.
3. **Confirm explicitly**, then open the release PR:
   ```bash
   gh pr create --base master --head dev \
     --title "Release v<version>" \
     --body "<Closes #N block from the aggregator issue> ...changelog summary..."
   ```

## Phase 6 — merge & tag

1. Wait for CI: `gh pr checks <pr-number>`. Only proceed once green.
2. **Confirm explicitly**, then: `gh pr merge <pr-number>`.
3. **Confirm explicitly**, then on `master`:
   ```bash
   git checkout master && git pull --ff-only
   git tag -a <version> -m "Release v<version>"
   git push origin <version>
   ```

## Phase 7 — publish

**Confirm explicitly**, then:

```bash
gh release create <version> --notes-file <path-to-changelog-section>
```

(Extract the `## v<version>` section of `CHANGELOG.md` to a temp file first,
or pass `--notes` inline.)

## Phase 8 — post-release `dev` bump

1. Back on `dev`: bump `manifest.version`/`.nf-core.yml` `template.version`
   to the next `<next>dev` (e.g. `1.1.0dev`).
2. Add a new empty "in progress" section at the top of `CHANGELOG.md`.
3. **Confirm explicitly**, then commit and push.

## Out of scope

No Zenodo API automation — Phase 4 is a manual checkpoint only.
