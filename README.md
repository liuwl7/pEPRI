# pEPRI

RNA-binding protein-mediated enhancer–promoter RNA interactions link noncoding variants to transcriptional dysregulation

This repository contains the analysis code and original Cytoscape session files used to generate the network visualizations presented in the manuscript.

Software: Cytoscape (v3.8.2) (https://cytoscape.org/)

Abbreviations
	pEPRI map: protein-mediated enhancer–promoter RNA interaction map
	pVTF map: protein-mediated variant-to-function map

0.reference/

	overlap.merge_vs_sub.activePT.txt:
		Promoter-target gene link file (hg19).

1.pEPRI_map_cytoscape_file/

	pEPRI_map.cys: 
		Cytoscape session displaying the pEPRI map for 52 RBPs.

2.pVTF_map_cytoscape_file/

	pVTF_map.cys: 
		Cytoscape session displaying the pVTF map for 52 RBPs.
		Two versions of the PTBP1 pVTF map are provided: 
		(1) PTBP1_pVTF_map_no_affinity_cut-off.cys (Figure S6A): No PTBP1 binding affinity fold-change threshold was applied.
		(2) PTBP1_pVTF_map_affinity_FC_1.5.cys (Figure 6D): Variants shown have predicted PTBP1 binding affinity fold change ≥ 1.5, resulting in a clearer and more interpretable pVTF map.

3.KEGG_cytoscape_file/

	PTBP1_pEPRI_mutation_KEGG_network.cys: 
		Cytoscape network corresponding to the KEGG pathway analysis of genes affected by mutations within PTBP1-mediated pEPRIs.
		Large nodes represent enriched KEGG pathway terms, whereas small nodes represent genes (Entrez Gene IDs).

These files are provided as source data accompanying the manuscript for academic and research purposes.
