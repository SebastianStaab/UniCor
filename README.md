## UniCor and UniCorP: A Novel Metric and Hierarchical Feature Selection Algorithm for Microbial Community Analysis

The accepted manuscript elucidating the UniCor metric and the UniCorP hierachical feature selection algorithm is available as open access here: https://doi.org/10.1093/ismeco/ycaf174 

The idea is to utilize the natural hierarchy in high dimensional, hierarchical datasets (like taxonomic hierarchy in microbiome datasets) in order to make them 
appropriate for a bigger variety of methods through a reduction of their feature space without the loss of relevant information. 

The UniCor Metric = |fcc| - ffc identifies UNIquely CORrelated eNtities (UNICORNs) with
- high absolute correlation (feature [cont. target var.] correlation, |fcc|)
- negative or low uniqueness (average feature feature correlation, ffc)

The UniCorP algorithm propagates UNICORNs through mutliple hierarchical levels to create enriched but more focused featuresets in higher hierarchies



## Input Format

UniCorP expects three input files:

1. **Feature Table** (`features.csv`)  
   A sample-by-feature matrix, where rows represent individual samples and columns represent features (e.g., ASVs or gene IDs). This table should contain raw counts or relative abundances depending on the intended transformation. It must not contain any non-numeric metadata columns.

2. **Target Variable** (`target.csv`)  
   A single-column table containing the continuous target variable (e.g., pH, temperature, biomass) for each sample. The index must match the sample IDs in the feature table.

3. **Hierarchy Table** (`hierarchy.csv`)  
   A feature-by-level matrix describing the biological or functional hierarchy of the features. Each column corresponds to a feature (matching the columns in the feature table), and each row represents a hierarchical level, from broad (e.g., Phylum) to fine-grained (e.g., Genus or ASV).  
   Missing values can be left blank or replaced with the feature ID itself. The order of levels should ideally go from highest (left) to lowest (right), although UniCorP can infer and adjust for reversed hierarchies.

### Format Summary

| File | Orientation | Index | Columns | Notes |
|------|-------------|-------|---------|-------|
| `features.csv` | Samples × Features | Sample IDs | Feature IDs | Numeric only |
| `target.csv`   | Samples × 1         | Sample IDs | Target name | One column |
| `hierarchy.csv`| Levels × Features  | Levels (e.g. Phylum, Class, ...) | Feature IDs | Categorical |


### Hierarchy Requirements

- The hierarchy must be a **strict tree**: each child node should have exactly **one parent**, and the number of nodes should decrease toward broader levels (e.g., Phylum > Class > Order).
- UniCorP cannot operate on **DAGs** (directed acyclic graphs) such as those found in **Gene Ontology** or **KEGG**, as these violate the strict parent-child structure.
- **Flat hierarchies** (where each node has only one child per level) provide no exploitable structure for propagation and are therefore not suitable.

### Notes
- Input files must be loaded as `pandas.DataFrame` objects with properly set indices.
- Use `.dropna()` or fill missing values as needed before applying transformations like CLR.
