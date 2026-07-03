# Purification Agent string format

Flow field: **Purification Agent** → API `purification_agent` (GraphQL `purificationAgent`).

## Preferred format (new uploads and post-upload patches)

```
{Species} Anti-{TARGET} ({Vendor} {Catalog})
```

Examples:

- `Mouse Anti-FLAG (Sigma F1804-200UG)`
- `Rabbit Anti-GRB2 (Bethyl A302-234A)`
- `Mouse Anti-V5 (Sigma V8012)`

Rules:

1. **Species** — `Mouse`, `Rabbit`, `Goat`, etc., when stated in Methods.
2. **Anti-{TARGET}** — gene/protein immunogen as in the paper (may differ from purification target symbol, e.g. CPSF5 / NUDT21).
3. **Vendor + catalog** in parentheses — use the exact catalog from Methods; include clone name in the catalog segment if that is how the paper lists it.
4. One primary antibody per IP sample. If Methods list two vendors for the same target, prefer the one used for eCLIP / IP in that sentence.
5. **INPUT / SMInput** — leave `purification_agent` **empty** (no IP antibody).

## What not to use

| Avoid | Use instead |
|-------|-------------|
| `V5-antibody` | Full vendor string from Methods, or anti-V5 catalog if that was the actual IP |
| `CPSF5 antibody` | `Mouse Anti-CPSF5 (Vendor Cat#)` from Methods |
| GEO boilerplate only | Paper Methods table |

## Agent extraction (Methods)

1. GEO matrix often says *“immunoprecipitated with RBP specific or V5-antibody”* — insufficient for Flow.
2. Read **eCLIP / CLIP Methods** (antibody table, “Antibodies”, or IP paragraph).
3. Build one string per **purification target**; apply to all IP replicates.
4. Record `evidence_quote` in proposals TSV for human review.

## Relation to Purification Target Annotation

| Sample type | `purification_target` | `purification_target__annotation` | `purification_agent` |
|-------------|----------------------|-----------------------------------|----------------------|
| eCLIP IP | Gene symbol | Tag if any (`cV5`, `FLAG`) | Antibody from Methods |
| Size-matched input | `SMInput` | empty | empty |

Tethered V5 screens (e.g. GSE290281): annotation may be `cV5` while **IP eCLIP** still uses the **RBP-specific** antibody listed in Methods for that protein.
