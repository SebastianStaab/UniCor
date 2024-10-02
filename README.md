UniCor - Hierarchical Feature Selection through Propagation of Uniquely Correlated Entities

The idea is to utilize the natural hierarchy in high dimensional, hierarchical datasets (like taxonomic hierarchy in microbiome datasets) in order to make them 
appropriate for a bigger variety of methods through a reduction of their feature space without the loss of relevant information. 

The UniCor Metric = |fcc| - ffc identifies UNIquely CORrelated eNtities (UNICORNs) with
- high absolute correlation (feature [cont. target var.] correlation, |fcc|)
- negative or low uniqueness (average feature feature correlation, ffc)

The UniCorP algorithm propagates UNICORNs through mutliple hierarchical levels to create enriched but more focused featuresets in higher hierarchies
