# Decay Analysis: Monegetti (Metamorphosis) vs. Randal (Settlement)

## Key Clarification: What Each Dataset Measures

### Monegetti Data: **Metamorphosis** (Spontaneous, No Cues)

**What is measured:**
- **Metamorphosis**: The biological transformation from larval to juvenile form
- **No settlement cues provided**: Larvae metamorphose spontaneously in control conditions
- **Data structure**: 
  - `larvae` = total larvae at start (20 per replicate, constant)
  - `meta` = number metamorphosed
  - `swimming` = number still swimming (not metamorphosed)
  - **No mortality tracking**: `larvae = meta + swimming` (all larvae accounted for)

**Biological process:**
- Metamorphosis is an **internal developmental process**
- Occurs without external settlement cues
- Represents the larva's **intrinsic competency** to transform
- Once metamorphosed, the larva has completed development

### Randal Data: **Settlement** (With Cues)

**What is measured:**
- **Settlement**: The behavioral act of attaching to a substrate
- **Settlement cues provided**: Multiple treatments (rubble, CCA, disc, etc.)
- **Data structure**:
  - `NoAlive` = number of larvae alive at observation
  - `NoSet` = number settled
  - `NoNotSet` = number not settled (but still alive)
  - **No separate mortality column**: `NoAlive = NoSet + NoNotSet` (mortality not explicitly tracked)

**Biological process:**
- Settlement is a **behavioral response** to environmental cues
- Requires appropriate settlement cues to trigger
- Represents the larva's **ability to respond** to settlement opportunities
- Larvae can remain alive but not settle if cues are absent or suboptimal

**Important Note**: While Randal's data structure shows `NoAlive = NoSet + NoNotSet`, this doesn't mean there's no mortality. It means:
- Mortality may be implicit (dead larvae are simply not counted in `NoAlive`)
- Or mortality is tracked separately and not included in these columns
- The data focuses on **settlement behavior** of **living larvae**

---

## Decay Parameter Comparison

### Late-Phase Exponential Decay (b2)

| Dataset | b2 Value | Interpretation | Biological Meaning |
|---------|----------|----------------|-------------------|
| **Monegetti (A. tenuis)** | 0.3969 | **Rapid decline** | Metamorphosis competency declines rapidly after change point |
| **Randal (Acroporidae)** | 0.001 | **Essentially no decline** | Settlement competency maintained after change point |

**Difference**: Monegetti shows **397× faster** late-phase decay than Randal

### Why the Dramatic Difference?

#### 1. **Biological Process Difference**

**Metamorphosis (Monegetti)**:
- Once a larva metamorphoses, it **cannot revert**
- Metamorphosis is a **one-way developmental process**
- After the optimal window, larvae that haven't metamorphosed may:
  - Lose the ability to metamorphose (developmental window closes)
  - Experience physiological decline
  - Die without completing development
- **Result**: Clear late-phase decline in metamorphosis competency

**Settlement (Randal)**:
- Settlement is a **behavioral choice** that can be delayed
- Larvae can remain **competent but not settled** for extended periods
- Settlement competency may be **maintained** if larvae remain healthy
- Larvae can settle later if appropriate cues become available
- **Result**: Maintained settlement competency with minimal decline

#### 2. **Cue Dependency**

**Monegetti (No Cues)**:
- Metamorphosis occurs spontaneously based on internal development
- If the developmental window passes, competency is lost
- No external cues can "rescue" larvae that missed the window
- **Decay reflects biological constraint**

**Randal (With Cues)**:
- Settlement depends on cue availability
- Larvae may remain competent but not settle due to:
  - Absence of appropriate cues
  - Suboptimal cue conditions
  - Waiting for better settlement opportunities
- Competency is maintained as long as larvae are alive and healthy
- **Decay reflects behavioral/ecological constraint, not biological loss**

#### 3. **Data Structure Effects**

**Monegetti**:
- Constant cohort size (20 larvae per replicate)
- Tracks metamorphosis directly
- Clear decline as larvae that haven't metamorphosed lose ability
- **Direct measurement of competency loss**

**Randal**:
- Variable cohort sizes across replicates
- Tracks settlement behavior, not intrinsic competency
- Decline may be masked by:
  - Multiple treatments averaging
  - Species aggregation
  - Behavioral variability
- **Indirect measurement of competency expression**

---

## Early-Phase Decay (b1, v1) Comparison

### Weibull Scale (b1)

| Dataset | b1 Value | Difference |
|---------|----------|------------|
| **Monegetti** | 0.001878 | Baseline |
| **Randal** | 0.0017 | -9.5% (very similar) |

**Interpretation**: Early-phase loss rates are **very similar**, suggesting comparable early competency dynamics regardless of metamorphosis vs. settlement.

### Weibull Shape (v1)

| Dataset | v1 Value | Difference |
|---------|----------|------------|
| **Monegetti** | 0.3645 | Baseline |
| **Randal** | 0.100 | -72.6% (much slower decay) |

**Interpretation**: Randal shows **much slower early-phase decay shape**, suggesting:
- More gradual competency loss in early phase
- Settlement competency may be more "flexible" than metamorphosis competency
- Behavioral competency (settlement) may be maintained longer than developmental competency (metamorphosis)

---

## Complete Decay Profile Comparison

### Early Phase (tc to Tcp: ~3.3-70 days)

**Monegetti (Metamorphosis)**:
- Moderate decay rate (b1 = 0.001878)
- Moderate decay shape (v1 = 0.3645)
- Gradual loss of metamorphosis competency
- **Biological constraint**: Developmental window closing

**Randal (Settlement)**:
- Similar decay rate (b1 = 0.0017)
- Much slower decay shape (v1 = 0.100)
- Very gradual loss of settlement competency
- **Behavioral/ecological constraint**: Cue availability, not biological loss

### Late Phase (After Tcp: >70 days)

**Monegetti (Metamorphosis)**:
- **Rapid exponential decay** (b2 = 0.3969)
- Clear decline from ~70% to near zero by day 80
- **Biological loss**: Larvae that haven't metamorphosed lose ability
- One-way process with clear endpoint

**Randal (Settlement)**:
- **Minimal decay** (b2 = 0.001)
- Maintained competency (~40-50%) through day 76
- **Behavioral maintenance**: Larvae remain competent but may not settle
- Reversible process with maintained potential

---

## Biological Interpretation

### Metamorphosis Decay (Monegetti)

**Mechanism**: Developmental window closing
- Metamorphosis is a **developmental milestone**
- Once the optimal window passes, physiological changes make metamorphosis less likely
- Larvae that don't metamorphose may:
  - Experience physiological decline
  - Lose developmental competency
  - Die without completing development
- **Result**: Clear late-phase decline in metamorphosis competency

**Ecological Implication**: 
- Larvae have a **limited window** for metamorphosis
- After ~70 days, metamorphosis becomes unlikely
- Late-phase decline reflects biological constraint

### Settlement Decay (Randal)

**Mechanism**: Behavioral/ecological constraint
- Settlement is a **behavioral response** to cues
- Competency may be maintained as long as larvae are alive
- Larvae that don't settle may:
  - Remain competent but waiting for cues
  - Experience reduced settlement response (not loss of competency)
  - Maintain settlement ability but choose not to settle
- **Result**: Maintained settlement competency with minimal decline

**Ecological Implication**:
- Larvae can **maintain settlement competency** for extended periods
- Settlement depends on cue availability, not just competency
- Late-phase maintenance reflects behavioral flexibility

---

## Implications for Connectivity Modeling

### 1. **Model Selection**

**For Metamorphosis-Based Models** (like Monegetti):
- Use Monegetti model with **high b2** (rapid late-phase decay)
- Captures developmental window closing
- Appropriate for models focusing on metamorphosis competency

**For Settlement-Based Models** (like Randal):
- Use simpler models (Weibull) or Monegetti with **low b2** (maintained competency)
- Captures behavioral settlement patterns
- Appropriate for models focusing on settlement behavior

### 2. **Parameter Interpretation**

**b2 Parameter Meaning**:
- **High b2 (Monegetti)**: Rapid loss of developmental competency
- **Low b2 (Randal)**: Maintained behavioral competency

**Model Selection Result**:
- Randal's Acroporidae: Weibull preferred (Δ AICc = +201.9)
- Suggests two-phase structure less important for settlement behavior
- Simpler models capture settlement patterns adequately

### 3. **Connectivity Modeling Applications**

**Metamorphosis Competency** (Monegetti approach):
- Use when modeling **developmental competency**
- Appropriate for understanding **intrinsic larval development**
- Captures **biological constraints** on metamorphosis
- **High b2** reflects developmental window closing

**Settlement Competency** (Randal approach):
- Use when modeling **settlement behavior**
- Appropriate for understanding **settlement responses**
- Captures **ecological constraints** on settlement
- **Low b2** reflects maintained behavioral competency

---

## Summary

### Key Differences in Decay

1. **Late-Phase Decay (b2)**:
   - **Monegetti**: 0.3969 (rapid decline) - **397× faster**
   - **Randal**: 0.001 (minimal decline) - essentially flat
   - **Difference**: Metamorphosis shows clear biological decline; settlement shows maintained competency

2. **Early-Phase Decay Shape (v1)**:
   - **Monegetti**: 0.3645 (moderate decay shape)
   - **Randal**: 0.100 (very slow decay shape) - **3.6× slower**
   - **Difference**: Settlement competency decays more gradually than metamorphosis competency

3. **Early-Phase Decay Rate (b1)**:
   - **Monegetti**: 0.001878
   - **Randal**: 0.0017
   - **Difference**: Very similar (only 9.5% difference)

### Biological Interpretation

**You are correct**:
- **Monegetti**: Measures **metamorphosis** (spontaneous, no cues)
- **Randal**: Measures **settlement** (with cues, includes behavioral response)

**About Mortality**:
- **Monegetti**: No mortality tracking - `larvae = meta + swimming` (all accounted for)
- **Randal**: No explicit mortality column - `NoAlive = NoSet + NoNotSet` (mortality not separately tracked, but may be implicit in `NoAlive`)

**Decay Differences Reflect**:
1. **Biological process**: Metamorphosis (developmental) vs. Settlement (behavioral)
2. **Cue dependency**: No cues (Monegetti) vs. With cues (Randal)
3. **Competency type**: Developmental competency vs. Behavioral competency
4. **Temporal dynamics**: Developmental window closing vs. Maintained behavioral potential

The dramatic difference in late-phase decay (b2) is the key finding: **metamorphosis competency declines rapidly after the optimal window, while settlement competency is maintained as long as larvae remain alive and healthy.**

---

## References

1. **Moneghetti et al. (2019)**: "High-frequency sampling and piecewise models reshape dispersal kernels of a common reef coral." *Ecology Letters*. [Metamorphosis data, no cues]

2. **Randal et al. (2024)**: "Age to settlement competency and settlement cue preference in coral larvae." *Communications Biology*. [Settlement data, with cues]

3. **This Analysis (2025)**: Direct comparison of decay parameters between metamorphosis and settlement datasets, revealing fundamental differences in competency dynamics.

---

**Report Date**: December 2, 2025  
**Key Finding**: Metamorphosis (Monegetti) shows rapid late-phase decay (b2 = 0.397), while settlement (Randal) shows maintained competency (b2 = 0.001), reflecting fundamental differences between developmental and behavioral processes.

