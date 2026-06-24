# Wiki schema — GUP integral research

This is a persistent, LLM-maintained wiki for the GUP (Generalized Uncertainty
Principle) integral project, following the pattern described at
https://gist.github.com/karpathy/442a6bf555914893e9891c11519de94f: rather than
re-deriving context from the raw scripts on every question, the wiki compiles
findings once and keeps them current.

## Layers

1. **Raw sources** — the `.py` scripts and plots at the repo root. Immutable;
   never edit them to "fix" the wiki.
2. **This wiki** (`wiki/`) — markdown pages summarizing what each script does,
   what was found, and how concepts relate.
3. **This file** — the schema governing how the wiki is structured.

## Structure

- `wiki/index.md` — catalog of all pages, grouped by category, one-line summary each.
- `wiki/log.md` — append-only chronological record of ingests and findings.
- `wiki/pages/*.md` — one page per script, concept, or open question.

## Page conventions

Each page in `wiki/pages/` starts with:

```markdown
# Title

**Source:** path/to/file.py (or "concept")
**Status:** active | superseded | open-question

One-paragraph summary.
```

Followed by free-form sections as needed (Key results, Open questions, Related pages).
Cross-link with relative markdown links, e.g. `[Stable Integrand](StableIntegrand.md)`.

## Operations

**Ingest** — when a script changes or a new one is added: re-read it, update
its page, update `index.md`'s one-line summary, append an `log.md` entry.

**Query** — when asked a question about the project, search `wiki/pages/`
first instead of re-reading every script. Cite the page. If the answer isn't
there, read the source, answer, then file the finding back as a page update.

**Lint** — periodically check for: stale claims (page says something the
current script no longer does), orphaned pages (not linked from `index.md`),
and contradictions between pages (e.g. two pages disagreeing about which N is
"recommended" for the same case).
