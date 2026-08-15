# Release review checklist

## Completed for the public GitHub review snapshot

- compact code and selected-result export;
- server/private path removal;
- raw data, environments, caches, TIFFs, logs and internal author documents excluded;
- deterministic manifest and integrity gate;
- public GitHub visibility for author and reviewer access.

## Corresponding-author decisions required before a versioned archive release

- final author order, affiliations, corresponding-author details and ORCIDs;
- final code licence and derived-output licence;
- whether the manuscript drafts belong in the public repository;
- release-tag timing relative to preprint or journal submission;
- release version (`v0.1.0` or `v1.0.0`);
- Zenodo archive and DOI authorization;
- whether an independent off-xulab clean execution or container is required.

## Proposed defaults

- code: Apache-2.0;
- original figures and derived summary tables: CC BY 4.0;
- third-party data: original repository terms;
- initial public tag: `v0.1.0` for the bounded pilot;
- archive: GitHub Release connected to Zenodo.
