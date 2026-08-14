# data/

Do not commit raw or large processed data.

Recommended layout:

```text
data/
  inventory/
  raw/          # gitignored
  processed/    # gitignored unless very small
  external/     # gitignored
  metadata/     # small curated metadata can be tracked
```

Track small metadata tables only when licenses allow it.

