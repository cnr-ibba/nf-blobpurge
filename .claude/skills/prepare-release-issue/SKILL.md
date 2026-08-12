---
name: prepare-release-issue
description: Draft (and, after explicit confirmation, create) a GitHub issue aggregating the issues/PRs that will be closed by the upcoming dev -> master release merge, as a base for CHANGELOG/release notes.
argument-hint: "[target-version]"
allowed-tools:
  - Bash(git log *)
  - Bash(git describe *)
  - Bash(gh pr list *)
  - Bash(gh pr view *)
  - Bash(gh issue list *)
  - Bash(gh issue view *)
disable-model-invocation: true
---

## Arguments

`$ARGUMENTS` is an optional target version (e.g. `1.0.0`). If omitted, infer
it by reading `manifest.version` in `nextflow.config` and dropping the
trailing `dev` suffix.

## What this does

Surveys everything that has landed on `dev` since the last release and
drafts a single GitHub issue listing:
- which existing issues will actually close on the `dev -> master` merge,
- which merged PRs are included,
- which PRs are still open/in-flight and therefore excluded unless merged
  before the cut,
- a ready-to-paste `Closes #N` block,
- a best-effort CHANGELOG draft grouped by category.

This is deliberately read-mostly. The only mutating action is a single,
explicitly-confirmed `gh issue create` at the end.

**Why this matters, not just documents**: GitHub only auto-closes an issue
when a `Closes/Fixes/Resolves #N` keyword lands in something merged into the
repo's **default branch** (`master` here), not `dev`. A PR merged into `dev`
with `Closes #N` in its body does **not** close that issue — confirmed in
this repo: issue #1 is still open even though PR #2 (merged into `dev`)
contains `Closes #1`. So the `Closes #N` block this skill produces is a
functional input to the `release` skill's `dev -> master` PR body, not just
release-notes material.

## Steps

1. **Find the last release tag** (read-only):
   ```bash
   git describe --tags --abbrev=0 origin/master 2>/dev/null || true
   ```
   If empty, there is no prior release — treat "since last release" as
   "since the start of the repo" (this is the case for the first release).

2. **List merged PRs targeting `dev`** since that tag's date (or all, if no
   tag exists):
   ```bash
   gh pr list --state merged --base dev \
     --json number,title,mergedAt,body,url,author --limit 200
   # if a last-release date is known, narrow with:
   #   --search "merged:>=<YYYY-MM-DD>"
   ```

3. **Cross-check against the actual commit range** reachable on `dev` but
   not `master`, to catch anything a search might miss:
   ```bash
   git log --oneline master..dev
   ```

4. **Extract closing-keyword issue references** from each merged PR body
   (case-insensitive, singular/plural/past tense, optional colon):
   ```
   \b(clos(e|es|ed)|fix(es|ed)?|resolv(e|es|ed))\s*:?\s*#[0-9]+\b
   ```

5. **Pull metadata for each referenced issue**, for categorization:
   ```bash
   gh issue view <N> --json number,title,state,labels,url
   ```

6. **List PRs open against `dev`** that are not yet merged, separately —
   these are "in-flight, excluded unless merged before the cut":
   ```bash
   gh pr list --state open --base dev --json number,title,body,url
   ```

## Categorization heuristic (best-effort — flag for human review)

Map existing repo labels to the CHANGELOG.md headings:

| Label | CHANGELOG bucket |
|---|---|
| `enhancement` | `### Added` |
| `bug` | `### Fixed` |
| `dependencies` | `### Dependencies` |
| `performance` | `### Fixed` (closest existing bucket — flag for review) |
| anything else / unlabeled | "Other / needs manual placement" |

State explicitly in the drafted issue that this categorization is a
heuristic based on label text, not reliable enough to paste verbatim into
`CHANGELOG.md` — it always needs a human pass.

## Draft issue body

```markdown
## Summary

N issue(s) will be closeable and M PR(s) merged into `dev` since <last tag | project start>.

## Issues that will close on the dev -> master merge

- [ ] #N <title> (<label>)
...

## Merged PRs included

- #N <title> (@author, merged <date>)
...

## Not yet merged (excluded from this release unless merged before the cut)

- #N <title> (open, targets dev, closes #M)
...

## Closes block for the release PR body

<!-- copy-paste ready: this is what the `release` skill's dev -> master PR body needs -->
Closes #N
Closes #M

## CHANGELOG draft (heuristic — needs manual cleanup)

### `Added`
- ...

### `Fixed`
- ...

### `Dependencies`

### `Deprecated`
```

## Confirm before creating

Never create the issue automatically. Print the full draft in chat, ask the
user to confirm, edit, or cancel. Only after explicit confirmation, run:

```bash
gh issue create --title "Release prep: v<version>" --body "$(cat <<'EOF'
<confirmed draft body>
EOF
)"
```

Use the `Release prep: v<version>` title prefix so the `release` skill can
find it later without needing a dedicated label:

```bash
gh issue list --search "in:title Release prep"
```

Report back the created issue's URL and number.
