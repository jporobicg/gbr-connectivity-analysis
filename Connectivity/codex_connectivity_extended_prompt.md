# Codex Prompt Add-on: Extended Connectivity Comparison and Uncertainty Analysis

Use this prompt **in addition to** the existing connectivity-comparison prompt. The goal is to extend the workflow so it includes stronger cross-disciplinary analyses of:

1. **system structure**
2. **strength of connection**
3. **overrepresentation of connections**
4. **uncertainty representation**

The core scientific question is:

> **How much better does the bootstrap analysis represent uncertainty compared with a single connectivity estimate?**

I want you to update the existing project so it includes the new analyses below, and also create a **very simple markdown summary file** that documents both the **old analyses** and the **new analyses**, with:
- the name of each analysis
- a **very short and simple explanation**
- the figure produced for that analysis
- a **short interpretation guide** for that figure

---

## 1. Main implementation request

Please extend the existing connectivity comparison workflow so it includes the following new analyses.

Keep the code:
- in **Python**
- simple
- modular
- easy to debug
- easy to inspect
- memory-aware

Do not overengineer the project.

Use standard scientific packages only, such as:
- `xarray`
- `numpy`
- `pandas`
- `matplotlib`
- `scipy`
- optionally `networkx`
- optionally `scikit-learn` for simple ordination if needed

Avoid complicated frameworks.

---

## 2. Important scientific framing

The analysis must distinguish between:

- the **single connectivity estimate**
- the **bootstrap mean connectivity**
- the **bootstrap distribution / uncertainty**

The goal is **not only** to compare matrices edge by edge.

The goal is also to assess whether the **ecological conclusions** change when using bootstrap-based connectivity instead of a single estimate.

In particular, test whether bootstraping changes:

- overall system structure
- identity of important source reefs
- identity of important sink reefs
- dominant pathways
- concentration of connectivity into a few strong links
- uncertainty in conclusions

---

## 3. New analyses to implement

# A. System structure analyses

## A1. Source and sink rank stability
For each time step:

- compute outgoing strength per source
- compute incoming strength per sink
- rank reefs by outgoing strength
- rank reefs by incoming strength

Compare:
- single matrix vs bootstrap mean
- single matrix vs bootstrap distribution
- bootstrap mean vs bootstrap samples

Metrics:
- Spearman rank correlation
- top-k overlap for sources and sinks
- absolute rank change
- largest movers up/down in ranking

Use several `k` values, for example:
- 10
- 50
- 100

Outputs:
- CSV table of source rank comparisons by time
- CSV table of sink rank comparisons by time
- tables of largest rank changes
- figure: time series of rank correlation
- figure: bar plot of largest source/sink rank changes for selected times

Interpretation goal:
- are the important sources/sinks stable?
- does bootstraping change which reefs matter most?

---

## A2. Destination-profile and source-profile change
For each source row:
- normalize outgoing connectivity across sinks
- compare the row-profile of:
  - single matrix
  - bootstrap mean

For each sink column:
- normalize incoming connectivity across sources
- compare the column-profile of:
  - single matrix
  - bootstrap mean

Use profile distances such as:
- Jensen-Shannon divergence
- Hellinger distance
- optionally Wasserstein distance if easy and justified

Outputs:
- per-time summary of mean/median row-profile divergence
- list of sources with biggest destination-profile change
- list of sinks with biggest origin-profile change
- figure: histogram of row-profile divergences
- figure: top changed sources by profile divergence

Interpretation goal:
- even if total connectivity is similar, do reefs connect to different places?

---

## A3. Higher-level structural grouping
Implement a simple and robust analysis of broader system organization.

Possible approaches:
- clustering using row sums / column sums / connectivity profiles
- thresholded network connected components
- simple community-like grouping if computationally manageable

Compare:
- grouping from single matrix
- grouping from bootstrap mean
- grouping stability across bootstrap samples

Outputs:
- summary table of group similarity across methods/times
- simple figure showing group-size summaries or similarity scores across time

Important:
Do **not** make this too complicated.
If full community detection is too heavy or unstable, use a simpler grouping approach and explain that choice in comments.

Interpretation goal:
- does the system break into similar large-scale patterns under both approaches?

---

# B. Strength of connection analyses

## B1. Concentration / inequality of edge strengths
Borrow ideas from economics.

For each time step, compute for:
- single matrix
- bootstrap mean
- bootstrap samples where feasible

Metrics:
- Lorenz curve
- Gini coefficient
- share of total connectivity held by:
  - top 1% edges
  - top 5% edges
  - top 10% edges
- 90th / 95th / 99th percentile edge values

Outputs:
- CSV table of concentration metrics by time
- figure: Lorenz curves for selected times
- figure: time series of Gini coefficient
- figure: top-share comparison over time

Interpretation goal:
- does the single estimate concentrate connectivity in a few strong links more than the bootstrap mean?
- does bootstrap averaging smooth extreme values?

---

## B2. Tail-behavior analysis
Compare the upper tail of edge weights between:
- single matrix
- bootstrap mean
- bootstrap sample distribution

Metrics:
- upper quantiles
- max edge value
- mean of top-k edges
- ratio of strongest edge to median positive edge
- ratio of top 1% mean to overall mean

Outputs:
- time summary table
- figure: tail metric time series
- figure: comparison of top-k edge strength distributions

Interpretation goal:
- are differences mainly in the strongest pathways?

---

## B3. Node-level strength uncertainty
For each node and time:
- outgoing strength
- incoming strength

Across bootstrap samples compute:
- mean
- SD
- CV
- quantiles / CI

Then compare the single value to bootstrap uncertainty.

Outputs:
- node-level uncertainty tables
- list of nodes with high uncertainty
- list of nodes where single estimate is far from bootstrap mean
- figure: selected node uncertainty intervals
- figure: histogram/distribution of node CV

Interpretation goal:
- are some reefs consistently important while others are highly uncertain?

---

# C. Overrepresentation / backbone analyses

## C1. Top-k edge overlap and churn
For each time:
- identify top-k strongest edges in the single matrix
- identify top-k strongest edges in the bootstrap mean

Compute:
- overlap count
- overlap proportion
- Jaccard similarity
- edges entering/leaving the top-k set

Use several k values:
- 100
- 500
- 1000

Outputs:
- CSV table by time and k
- figure: overlap vs time
- figure: selected-time churn plot

Interpretation goal:
- do the strongest pathways stay the same or change?

---

## C2. Stable-edge frequency across bootstrap samples
For each edge and time:
- determine how often the edge appears in:
  - top-k set
  - above-threshold set
  - optionally a local-significance backbone if implemented

Then define stable edges, for example:
- present in at least 50% of bootstrap samples
- present in at least 80% of bootstrap samples
- present in at least 95% of bootstrap samples

Compare these stable edges with:
- top edges in the single matrix
- top edges in the bootstrap mean

Outputs:
- frequency tables
- stable-edge summary by time
- figure: stable-edge counts vs threshold/frequency rule
- figure: overlap of single matrix with stable bootstrap edges

Interpretation goal:
- which pathways are robust to sampling uncertainty?

---

## C3. Local dominance / entropy of outgoing and incoming patterns
For each source and time:
- normalize outgoing edge weights into proportions
- compute:
  - entropy
  - effective number of destinations
  - dominance ratio (largest outgoing edge / total outgoing strength)

For each sink and time:
- normalize incoming edge weights into proportions
- compute:
  - entropy
  - effective number of sources
  - dominance ratio

Compare:
- single matrix
- bootstrap mean
- bootstrap uncertainty where feasible

Outputs:
- node summary tables
- figure: entropy vs total strength
- figure: dominance-ratio comparison
- figure: distribution of effective number of destinations/sources

Interpretation goal:
- does bootstraping make connectivity patterns look more diffuse or more concentrated?

---

## C4. Optional backbone extraction
If computationally manageable, implement a simple backbone extraction.

Options:
- threshold backbone
- local-significance / disparity-style backbone
- row-share dominance backbone

Important:
Only use a more advanced backbone method if it can be implemented clearly and robustly.
Otherwise use a transparent simpler proxy.

Compare:
- backbone of single matrix
- backbone of bootstrap mean
- stable backbone across bootstrap samples

Outputs:
- backbone edge counts
- overlap tables
- figure: backbone overlap by time

Interpretation goal:
- do both approaches identify the same important structural links?

---

# D. Uncertainty and robustness analyses

## D1. Coverage-style comparison
For each time, compare the single estimate with the bootstrap uncertainty.

At minimum do this for:
- edges (possibly thresholded or sampled if full matrix is too large)
- outgoing node strength
- incoming node strength
- Gini coefficient
- connectance
- top-k overlap metrics if feasible

Compute:
- bootstrap CI / quantile interval
- whether the single value falls inside or outside the interval

Outputs:
- coverage summary table
- proportion outside CI by metric and time
- figure: coverage rates across time
- figure: selected example intervals

Interpretation goal:
- is the single estimate consistent with the uncertainty implied by bootstrap samples?

---

## D2. Compare single-to-bootstrap discrepancy with within-bootstrap variability
For each time:

- compute similarity/distance among bootstrap samples
- compute similarity/distance between the single matrix and the bootstrap mean
- compare these two scales

Possible similarity metrics:
- Pearson correlation
- Spearman correlation
- Frobenius norm
- RMSE
- row-sum rank correlation
- top-k overlap

Outputs:
- summary table by time
- figure: single-vs-bootstrap difference compared with within-bootstrap variability

Interpretation goal:
- is the single estimate just another plausible realization?
- or is it outside the normal variability of the bootstrap ensemble?

This is a key analysis.

---

## D3. Bootstrap rank uncertainty
For each time:
- for each node, compute bootstrap distribution of:
  - outgoing rank
  - incoming rank

Then calculate:
- rank median
- rank interval
- rank SD
- probability of being in top-k

Compare with the single estimate rank.

Outputs:
- rank-uncertainty tables
- top-k probability tables
- figure: probability of being top-k for selected nodes
- figure: rank uncertainty intervals for important nodes

Interpretation goal:
- are “important” reefs truly stable, or only apparently important in one estimate?

---

## D4. Ordination / low-dimensional uncertainty view
Implement a simple ordination of bootstrap replicates if computationally feasible.

Possible inputs:
- row sums
- column sums
- top-link summary vectors
- thresholded edge summaries

Possible methods:
- PCA
- classical MDS
- non-metric MDS only if easy

For each time:
- represent bootstrap samples in low-dimensional space
- place the single estimate and bootstrap mean in the same space

Outputs:
- figure: bootstrap cloud with single estimate and bootstrap mean overlaid

Interpretation goal:
- does the single estimate sit inside the bootstrap cloud or outside it?

Important:
Do this only if manageable.
Use a simple feature representation.
Do not flatten the entire matrix if that is too expensive.

---

# E. Threshold sensitivity

Repeat selected important metrics under several thresholds for defining an active edge.

Use thresholds such as:
- 0
- a very small positive value
- several progressively stronger values
- optionally data-informed quantiles

Apply threshold sensitivity to:
- connectance
- top-k edge overlap
- backbone edge counts
- stable-edge frequency
- concentration metrics where relevant

Outputs:
- threshold sensitivity tables
- figure: threshold sensitivity curves

Interpretation goal:
- are conclusions robust to what counts as a meaningful connection?

---

## 4. Required summary markdown file

In addition to the code and figures, create a markdown file called something like:

```text
outputs/analysis_catalog.md
```

This file must summarize **both the original analyses and the new analyses**.

For each analysis include a very short block with this structure:

```md
## Analysis name

**What it does:** one or two very short sentences in simple language.

**Figure:** filename of the figure produced.

**How to read it:** one or two very short sentences explaining what a high/low value or pattern means.
```

Keep the wording:
- very short
- very simple
- easy to scan

Do not write a long report.
This is a compact catalog of analyses and figures.

Also include:
- a short introduction at the top
- a very short “recommended figures for a paper/report” section at the end

---

## 5. Also update the project report

Please update the existing reporting workflow so the final report includes:
- old analyses
- new analyses
- a short paragraph on what each analysis adds
- a final synthesis focused on uncertainty

The report should explicitly address:

1. Are the single and bootstrap-average matrices broadly similar?
2. Where do they differ most?
3. Do bootstraps change the identity of key sources, sinks, or pathways?
4. Are differences concentrated in weak links or strong links?
5. Does bootstraping materially improve representation of uncertainty?

---

## 6. Output requirements

Please save outputs in a clear structure such as:

```text
outputs/
├── tables/
├── figures/
├── intermediate/
└── analysis_catalog.md
```

At minimum produce:
- CSV summary tables
- figure files
- analysis catalog markdown file
- updated report or summary markdown

---

## 7. Efficiency and fallback rules

These matrices are large, so follow these rules:

- process one time step at a time when needed
- avoid unnecessary full dense copies
- use chunking/lazy loading where helpful
- save intermediate summaries

If a full edge-level analysis is too expensive:
- use thresholding
- use top-k edges
- use random subsampling for visualization only
- use node-level summaries
- explain the fallback clearly in comments

Always prefer:
- interpretable
- robust
- biologically meaningful
over computationally flashy.

---

## 8. Coding style

Use:
- small functions
- clear names
- docstrings
- comments
- no unnecessary classes
- no overcomplicated abstractions

The project should still feel like a clean scientific analysis workflow.

---

## 9. Final deliverables

After implementing the update, provide:

1. updated code
2. updated README
3. `analysis_catalog.md`
4. list of new figures produced
5. list of new CSV outputs
6. short note on which analyses are most informative for:
   - structure
   - strength
   - overrepresentation
   - uncertainty

---

## 10. Important scientific guidance for interpretation

Please keep these points in mind while implementing and summarizing results:

- the bootstrap mean may smooth stochastic variation
- the single estimate may exaggerate sharp pathways
- differences in weak edges may matter less than differences in dominant exporters/importers
- uncertainty should be interpreted at:
  - edge level
  - node level
  - whole-system level
- the main question is whether ecological conclusions are robust, not whether every edge is identical

---

## 11. Do not stop for minor decisions

Make sensible defaults.
Only stop if a decision would fundamentally change the scientific interpretation or make the analysis invalid.

---

## 12. Final note

The central idea is:

> Test whether a single connectivity estimate gives the same ecological interpretation as the bootstrap ensemble, and quantify how much extra uncertainty information the bootstrap approach provides.

Please keep that central goal visible throughout the implementation, figures, and markdown summaries.
