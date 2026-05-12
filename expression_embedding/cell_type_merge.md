# C. elegans Cell Type Merger Schemes

## Original Distribution (reference)

| Cell type | Acc | n |
|---|---|---|
| marginal | 1.000 | 9 |
| hypoderm | 1.000 | 51 |
| excretory | 1.000 | 4 |
| rectal | 1.000 | 2 |
| tail | 1.000 | 2 |
| intestine | 1.000 | 2 |
| mesoderm | 1.000 | 1 |
| repro | 1.000 | 2 |
| muscle | 0.982 | 110 |
| neuron | 0.964 | 167 |
| socket | 0.944 | 18 |
| other | 0.889 | 9 |
| valve | 0.875 | 8 |
| programmed_death | 0.865 | 74 |
| sheath | 0.864 | 22 |
| epithelium | 0.538 | 13 |

---

## Scheme A — Mirror Large et al. (~7–8 classes, recommended starting point)

| Merged class | Components | n | Notes |
|---|---|---|---|
| Neuron | neuron | 167 | |
| Muscle | muscle | 110 | |
| Programmed death | programmed_death | 74 | Keep separate — fate decision, not tissue identity |
| Skin (hypodermis + seam) | hypoderm + tail + epithelium* | 66 | *Inspect epithelium lineages before merging |
| Glia + excretory | sheath + socket + excretory | 44 | Convention from Large et al. and Packer et al. |
| Pharynx + rectal | marginal + valve + rectal | 19 | Both ends of gut share PHA-4-dependent programs |
| Mesoderm + reproductive | mesoderm + repro | 3 | Too small to stand alone |
| Intestine | intestine | 2 | Still very small — consider dropping or folding into Pharynx+rectal |
| *(drop)* | other | — | Unspecified — exclude from supervised tasks |

---

## Scheme B — Strict Ma et al. style (5 classes + drop)

| Merged class | Components | n | Notes |
|---|---|---|---|
| Neuron | neuron | 167 | |
| Muscle | muscle | 110 | |
| Skin | hypoderm + tail | 53 | |
| Pharynx | marginal + valve | 17 | Matches Ma et al. pharynx-specific TF set (20 TFs) |
| Intestine | intestine + rectal | 4 | Still small |
| *(drop)* | sheath + socket + excretory + mesoderm + repro + epithelium + programmed_death + other | ~150 | Cleanest for using Ma et al. TF sets without circularity |

---

## Scheme C — Pragmatic 5-class keeping apoptotic

| Merged class | Components | n | Notes |
|---|---|---|---|
| Neuron | neuron | 167 | |
| Muscle | muscle | 110 | |
| Programmed death | programmed_death | 74 | Tests whether TF profile predicts apoptosis fate prospectively |
| Epithelial / skin | hypoderm + tail + epithelium + rectal | 68 | *Inspect epithelium + rectal cell identities before merging |
| Glia + excretory | sheath + socket + excretory | 44 | |
| Pharynx | marginal + valve | 17 | |
| *(drop or other)* | intestine + mesoderm + repro + other | 14 | |

# Schema D - feedback from colleague
neuron -> neuron
muscle -> muscle
reproduction -> reproduction
other-arc, hypoderm, epithelium + seam -> epithelium
sheath + socket -> glial
coelomocyte -> coelomocyte
excretory -> excretory
mesoderm -> mesoderm
other -> other
intestine, valve, marginal, gland, rectal  ->  alimentary
tail spike -> die anyways so can be bundled with programmed_death