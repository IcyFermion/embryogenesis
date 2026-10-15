# Track C: outlier subtrees

Collaborators asked whether the subtree outliers carry a biological story.
Two steps: a data dossier from existing results, then a literature research
session that uses the dossier as input. Start after B1, which provides the
edge-change tables.

## Separate numerical artifacts first

Primary cohort (Table 1):

| Subtree | n | On front | u_L | d_LP | d_CP | Max ER | Reading |
|---|---|---|---|---|---|---|---|
| ABplpa | 17 | yes | 1.000 | 0.000 | 64.0 | 1.000 | Natural lineage on its front; the large d_CP and u = 1 come from a near-zero endpoint span. Probably numerical, like C. briggsae P2 in the cross-species sketches. Confirm from raw costs. |
| ABarp | 15 | no | 0.000 | 0.165 | 0.386 | 0.800 | Genuinely off front. |
| ABplppp | 12 | no | 0.000 | 0.169 | 0.450 | 0.750 | Genuinely off front. |

Also review other off-front subtrees (ABplpp d_LP 0.119; ABplp 0.061) and any
with tiny endpoint spans (ABalpa d_CP 18.8; ABprpa 3.4).

## Dossier, per subtree

- Edges that differ at P* and at maximum retention: child, natural and
  assigned parents, tree distance, raw travel and state changes.
- Fates, positions and division timing of those cells; whether they end at the
  tracking cutoff in any embryo.
- Robustness: per tracking embryo (S3 caches), 3D vs XY, and in the 187-edge
  cross-species cohort (CE protein, CE RNA, CB RNA) where the subtree survives
  (ABarp and ABplpa have 12 matched edges each; ABplppp falls below the
  12-edge threshold there).
- Raw cost differences, not only endpoint-normalized distances.

## Research session

A separate session with web search, given the dossier: lineage, fates,
migrations, cell deaths and known peculiarities of these branches (Sulston et
al. 1983; WormAtlas, WormBase). The goal is candidate explanations for the
collaborators, each tied to specific edges in the dossier.
