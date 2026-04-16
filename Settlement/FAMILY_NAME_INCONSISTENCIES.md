# Family Name Inconsistencies Report

## Summary

**The Species abbreviations.csv file was NOT used** to classify families. The analysis used family names directly from `Settlement.csv`, which has some inconsistencies with the Species abbreviations file.

## Inconsistencies Found

### 1. **Diploastraeidae vs Diploastreidae** ⚠️

- **Species abbreviations.csv**: `Diploastraeidae`
- **Settlement.csv**: `Diploastreidae`
- **Status**: Spelling difference (missing 'a' in Settlement.csv)
- **Impact**: Same family, different spelling

### 2. **Lobophyllidae vs Lobophylliidae** ⚠️

- **Species abbreviations.csv**: `Lobophyllidae`
- **Settlement.csv**: `Lobophylliidae`
- **Status**: Spelling difference (extra 'i' in Settlement.csv)
- **Impact**: Same family, different spelling

### 3. **Agariciidae** ⚠️

- **Species abbreviations.csv**: NOT PRESENT
- **Settlement.csv**: PRESENT
- **Status**: Family in Settlement data but not in Species abbreviations file
- **Impact**: Missing from reference file

### 4. **Pachyseridae** ℹ️

- **Species abbreviations.csv**: PRESENT
- **Settlement.csv**: NOT PRESENT
- **Status**: Family in Species abbreviations but no data in Settlement.csv
- **Impact**: No data to analyze, not a problem

## Current Analysis Status

**Families analyzed in our results:**
1. Acroporidae ✅ (matches both files)
2. Agariciidae ⚠️ (in Settlement.csv, not in Species abbreviations)
3. Diploastreidae ⚠️ (spelling: should be Diploastraeidae)
4. Euphylliidae ✅ (matches both files)
5. Lobophylliidae ⚠️ (spelling: should be Lobophyllidae)
6. Merulinidae ✅ (matches both files)
7. Poritidae ✅ (matches both files)

## Recommendations

### Option 1: Standardize to Species abbreviations.csv (Recommended)

Update `Settlement.csv` or create a mapping to standardize:
- `Diploastreidae` → `Diploastraeidae`
- `Lobophylliidae` → `Lobophyllidae`
- Keep `Agariciidae` (add to Species abbreviations if needed)

### Option 2: Use Settlement.csv as source of truth

Accept Settlement.csv names as correct and update Species abbreviations.csv:
- `Diploastraeidae` → `Diploastreidae`
- `Lobophyllidae` → `Lobophylliidae`
- Add `Agariciidae` to Species abbreviations

### Option 3: Create a mapping function

Create a family name mapping function in `utility_tools.py` to standardize names during data loading.

## Impact on Analysis

**Current impact**: **MINIMAL**
- The analysis is internally consistent (uses Settlement.csv names throughout)
- Results are valid for the families as named in Settlement.csv
- No data loss or incorrect grouping

**Potential issues**:
- If comparing with other datasets using Species abbreviations names, there will be mismatches
- Publication should clarify which naming convention is used
- May cause confusion when referencing species-level data

## Next Steps

1. **Decide on naming convention**: Which file is the source of truth?
2. **Create mapping function**: If standardizing, add to `utility_tools.py`
3. **Update analysis**: Re-run analysis with standardized names if needed
4. **Document decision**: Add to analysis documentation

## Code Changes Needed (if standardizing)

If we standardize to Species abbreviations.csv names, we would need to:

1. Add mapping function to `utility_tools.py`:
```python
def standardize_family_names(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize family names to match Species abbreviations.csv"""
    mapping = {
        'Diploastreidae': 'Diploastraeidae',
        'Lobophylliidae': 'Lobophyllidae'
    }
    df['Family'] = df['Family'].replace(mapping)
    return df
```

2. Update `load_settlement_data()` to use this mapping

3. Re-run analysis with standardized names

---

**Status**: ✅ **STANDARDIZED** - Family names are now automatically standardized to match Species abbreviations.csv during data loading.

## Implementation

Family name standardization has been implemented in `utility_tools.py`:

1. **`get_family_name_mapping()`**: Creates mapping from Settlement.csv names to Species abbreviations.csv names
2. **`standardize_family_names()`**: Applies the mapping to a DataFrame
3. **`load_settlement_data()`**: Now automatically standardizes family names by default

### Standardization Applied

- `Diploastreidae` → `Diploastraeidae` (539 rows)
- `Lobophylliidae` → `Lobophyllidae` (760 rows)

### Note on Agariciidae

`Agariciidae` remains unchanged as it is present in Settlement.csv but not in Species abbreviations.csv. This is acceptable as it represents a family with data but not in the reference file.

## Next Steps

**⚠️ IMPORTANT**: The analysis should be re-run with standardized family names to ensure all results use consistent naming.


