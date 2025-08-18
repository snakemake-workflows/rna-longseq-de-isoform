# Multi-Factor Differential Expression Analysis

This document explains how to use the enhanced differential expression analysis script that supports multiple categorical variables.

## Overview

The enhanced script now supports four analysis modes:

1. **Traditional Binary Analysis** (default): Compares exactly two conditions in a single factor
2. **Multi-Factor Analysis**: Performs pairwise comparisons within each categorical factor independently
3. **Cross-Factor Analysis**: Performs pairwise comparisons between combinations of factor levels across different factors
4. **Unknown Value Handling** (new): Handles unknown/missing values with multiple analysis scenarios

## Unknown Value Handling

### Supported Unknown Indicators
The script recognizes these values as "unknown":
- `"unknown"`
- `"not available"`
- `"na"` or `"n/a"` 
- `"NA"`
- Empty strings `""`

### Analysis Scenarios
For each factor containing unknown values, the script runs **three separate analyses**:

1. **Unknown as Factor Level 1**: Unknown samples are grouped with the first known factor level
2. **Unknown as Factor Level 2**: Unknown samples are grouped with the second known factor level  
3. **Exclude Unknown**: Unknown samples are completely excluded from analysis

### Example
Given samples with `condition` = `[male, female, unknown, male, female, unknown]`:

- **Scenario 1**: `unknown` → `male` → `[male, female, male, male, female, male]`
- **Scenario 2**: `unknown` → `female` → `[male, female, female, male, female, female]`
- **Scenario 3**: Exclude unknowns → `[male, female, male, female]`

Each scenario generates complete DE analysis results, allowing you to assess how sensitive your findings are to the treatment of unknown values.

## Types of Comparisons

### Within-Factor Analysis (multi_factor_analysis: true)
- **condition**: control vs case
- **treatment**: treated vs untreated  
- **genotype**: WT vs KO

### Cross-Factor Analysis (cross_factor_analysis: true)
- **condition × treatment**: 
  - control+treated vs control+untreated
  - control+treated vs case+treated
  - control+treated vs case+untreated
  - control+untreated vs case+treated
  - control+untreated vs case+untreated
  - case+treated vs case+untreated

- **condition × genotype**:
  - control+WT vs control+KO
  - control+WT vs case+WT
  - control+WT vs case+KO
  - control+KO vs case+WT
  - control+KO vs case+KO
  - case+WT vs case+KO

- **treatment × genotype**:
  - treated+WT vs treated+KO
  - treated+WT vs untreated+WT
  - treated+WT vs untreated+KO
  - treated+KO vs untreated+WT
  - treated+KO vs untreated+KO
  - untreated+WT vs untreated+KO

## Configuration

### Enable All Features

Add the following to your `config.yml` under the `deseq2` section:

```yaml
deseq2:
    # ... other parameters ...
    
    # Enable within-factor analysis
    multi_factor_analysis: true
    
    # Enable unknown value handling (NEW!)
    handle_unknown_values: true
    
    # Enable cross-factor analysis (combinations between factors)
    cross_factor_analysis: true
    
    # Maximum number of factors to combine (default: 2)
    max_cross_factors: 2
    
    # Optional: specify which factors to analyze
    factors_to_analyze:
        - "condition"
        - "condition2"
    
    # Include all factors in the design
    design_factors:
        - "condition"
        - "condition2"
        - "batch_effect"  # control variable
```

**Important:** 
- Unknown handling creates 3x more comparisons per factor with unknowns
- Cross-factor analysis uses "exclude unknown" scenario for cleaner combinations
- Total output files can grow very quickly with multiple features enabled

### Automatic Factor Detection

If you don't specify `factors_to_analyze`, the script will automatically analyze all categorical factors (excluding those listed in `continuous_factors`):

```yaml
deseq2:
    multi_factor_analysis: true
    design_factors:
        - "condition"
        - "treatment"
        - "batch"
    continuous_factors:
        - "age"
        - "weight"
    # factors_to_analyze not specified = analyze condition, treatment, batch
```

## Sample File Requirements

Your `samples.csv` should include all the factors you want to analyze, and can contain unknown values:

```csv
sample,condition,condition2,batch_effect,platform,purity
A,male,young,batch1,NANOPORE,0.98
B,female,young,batch1,NANOPORE,0.97
C,male,old,batch1,NANOPORE,0.96
D,female,old,batch1,NANOPORE,0.99
E,unknown,middle,batch2,NANOPORE,0.95
F,male,unknown,batch2,NANOPORE,0.98
G,female,na,batch2,NANOPORE,0.97
H,not available,young,batch2,NANOPORE,0.96
```

**Supported Unknown Values:**
- `unknown`
- `not available` 
- `na`, `n/a`, `NA`
- Empty strings

**Requirements:**
- Each factor must have at least 2 known (non-unknown) levels
- At least 4 samples total for meaningful analysis
- Unknown values are handled automatically when `handle_unknown_values: true`

## Output Files

### Multi-Factor Analysis Outputs

When `multi_factor_analysis: true`, the script generates:

1. **Within-Factor Comparison Files**: For each factor analyzed independently
   - Standard: `de_analysis_{factor}_{level1}_vs_{level2}_lfc_analysis.csv`
   - With unknowns: `de_analysis_{factor}_unknown_as_{level1}_{level1}_vs_{level2}_lfc_analysis.csv`
   - Excluding unknowns: `de_analysis_{factor}_exclude_unknown_{level1}_vs_{level2}_lfc_analysis.csv`

2. **Cross-Factor Comparison Files**: (when `cross_factor_analysis: true`)
   - `de_analysis_cross_{factor1}_{factor2}_{combo1}_vs_{combo2}_lfc_analysis.csv`
   - Uses "exclude unknown" scenario for cleaner combinations

3. **Unknown Handling Outputs**: (when `handle_unknown_values: true`)
   - **Three scenarios per factor** with unknown values
   - **Sensitivity analysis** across different treatments of unknowns

4. **Summary File**: `de_analysis_summary.csv`
   - Overview of all comparisons with counts of significant genes

5. **Main Output Files**: (using the comparison with most significant genes)
   - Standard output files for compatibility with downstream analysis

### Example Output Structure

```
results/
# Standard within-factor comparisons
├── de_analysis_condition_male_vs_female_lfc_analysis.csv
├── de_analysis_condition2_young_vs_old_lfc_analysis.csv

# Unknown handling scenarios
├── de_analysis_condition_unknown_as_male_male_vs_female_lfc_analysis.csv
├── de_analysis_condition_unknown_as_female_male_vs_female_lfc_analysis.csv
├── de_analysis_condition_exclude_unknown_male_vs_female_lfc_analysis.csv
├── de_analysis_condition2_unknown_as_young_young_vs_old_lfc_analysis.csv
├── de_analysis_condition2_unknown_as_old_young_vs_old_lfc_analysis.csv
├── de_analysis_condition2_exclude_unknown_young_vs_old_lfc_analysis.csv

# Cross-factor comparisons (using exclude_unknown data)
├── de_analysis_cross_condition_condition2_male_young_vs_male_old_lfc_analysis.csv
├── de_analysis_cross_condition_condition2_male_young_vs_female_young_lfc_analysis.csv
├── de_analysis_cross_condition_condition2_female_young_vs_female_old_lfc_analysis.csv

# Summary and main outputs
├── de_analysis_summary.csv
└── de_analysis_lfc_analysis.csv  # Main output (best comparison)
```

## Benefits

### Normalization Efficiency
- DESeq2 normalization is performed **once** for the entire dataset
- All pairwise comparisons use the same normalized data
- Significant computational time savings compared to running separate analyses

### Comprehensive Analysis
- **Within-factor analysis**: All pairwise comparisons within each factor
- **Cross-factor analysis**: All pairwise comparisons between factor combinations
- **True factorial design**: Captures interaction effects between factors
- Provides summary statistics for easy comparison of effect sizes
- Maintains compatibility with existing downstream analysis

### Flexible Configuration
- Can specify exactly which factors to analyze
- Control complexity with `max_cross_factors`
- Automatic detection of categorical vs continuous factors
- Backwards compatible with existing binary analysis workflows

## Migration from Binary Analysis

### Old Configuration (Binary)
```yaml
deseq2:
    design_factors:
        - "condition"
    # Only analyzes condition with exactly 2 levels
```

### New Configuration (Within-Factor Only)
```yaml
deseq2:
    multi_factor_analysis: true
    design_factors:
        - "condition"
        - "treatment"
        - "genotype"
    factors_to_analyze:
        - "condition"
        - "treatment"
        - "genotype"
    cross_factor_analysis: false  # or omit
```

### New Configuration (With Cross-Factor)
```yaml
deseq2:
    multi_factor_analysis: true
    cross_factor_analysis: true
    max_cross_factors: 2
    design_factors:
        - "condition"
        - "treatment"
        - "genotype"
    factors_to_analyze:
        - "condition"
        - "treatment"
        - "genotype"
```

## Best Practices

1. **Include Control Variables**: Add batch effects and other confounding variables to `design_factors` even if you don't want to analyze them directly.

2. **Specify Continuous Factors**: List quantitative variables in `continuous_factors` to exclude them from pairwise analysis.

3. **Check Sample Balance**: Ensure adequate sample sizes for each level of your factors.

4. **Review Summary**: Use the summary file to identify the most informative comparisons.

## Troubleshooting

### Common Issues

1. **Factor not found in metadata**: 
   - Check spelling in `factors_to_analyze`
   - Verify factor exists in your samples.csv

2. **Insufficient levels for comparison**:
   - Each factor must have at least 2 levels
   - Check for typos or inconsistent naming in samples.csv

3. **LFC shrinkage warnings**:
   - Normal for some comparisons
   - Analysis continues without shrinkage if coefficient not found

### Error Messages

- `Factor 'X' not found in metadata columns`: Check factor name spelling
- `Factor 'X' must have at least 2 levels`: Verify factor has multiple unique values
- `Warning: Could not perform LFC shrinkage`: Non-critical, analysis continues

## Performance Considerations

- **Within-factor analysis** scales as O(k²) where k is the number of levels per factor
- **Cross-factor analysis** scales exponentially: O((k₁×k₂×...×kₙ)²) for n factors
- Large numbers of factors or levels may require increased computational resources
- **Recommendations**:
  - Use `factors_to_analyze` to limit analysis to factors of primary interest
  - Start with `max_cross_factors: 2` and increase cautiously
  - Consider biological significance before adding more factors
  - Monitor memory usage for large datasets with many factor combinations

### Complexity Examples

| Factors | Levels Each | Within-Factor Comparisons | Cross-Factor Comparisons (2-way) |
|---------|-------------|---------------------------|-----------------------------------|
| 2       | 2           | 2                        | 6                                 |
| 2       | 3           | 6                        | 36                                |
| 3       | 2           | 3                        | 21                                |
| 3       | 3           | 9                        | 153                               |
| 4       | 2           | 4                        | 66                                |

**Warning**: 3+ factor combinations grow very quickly. A 3×3×3 design creates 729 possible combinations!
