# pEPRI

This repository contains the analysis code and Cytoscape session files used to generate the network visualizations presented in the manuscript.


	RNA-binding protein-mediated enhancer–promoter RNA interactions link noncoding variants to transcriptional dysregulation.

## **Analysis Workflow**

<img width="432" height="610" alt="pipeline" src="https://github.com/user-attachments/assets/bb6e1a19-dac0-4450-b581-304cb4d7e0dd" />

## **Dependencies and Software**
	
	Ubuntu 22.04 LTS
	R (v4.3.1)
	
	Trimmomatic (v0.39) (http://www.usadellab.org/cms/?page=trimmomatic)
	cutadapt (v4.4) (https://cutadapt.readthedocs.io/en/stable)
	STAR (v2.5.2b) (https://github.com/alexdobin/STAR)
	bwa (v0.7.17) (https://github.com/lh3/bwa)
	SAMtools (v1.10) (https://github.com/samtools/samtools)
	BEDtools (v2.27.1) (https://github.com/arq5x/bedtools2)
	RICpipe (v1.0) (https://github.com/caochch/RICpipe)
	CRIC-seq (v1.0) (https://github.com/HuNaijing/CRIC-seq)
	Cytoscape (v3.8.2) (https://cytoscape.org/)
	
## **Scripts for pEPRI identification**

	Remove adapters and PCR duplicates, STAR map to genome and remove background (IgG) for CRIC-seq data 
	using the RICpipe (https://github.com/caochch/RICpipe) and CRIC-seq (https://github.com/HuNaijing/CRIC-seq) pipelines.
	After removed background, inter-molecular chimeric reads are mapped to enhancers and promoters using the following scripts:

	

## **Mapping risk variants onto pEPRIs**


## **Cytoscape session files for pEPRI map, pVTF map and KEGG network**

**Abbreviations**

	pEPRI map: protein-mediated enhancer–promoter RNA interaction map
	
	pVTF map: protein-mediated variant-to-function map

**0.reference/**

	overlap.merge_vs_sub.activePT.txt:
		Promoter-target gene link file (hg19).

**1.pEPRI_map_cytoscape_file/**

	pEPRI_map.cys: 
		Cytoscape session displaying the pEPRI map for 52 RBPs.

**2.pVTF_map_cytoscape_file/**

	pVTF_map.cys: 
		Cytoscape session displaying the pVTF map for 52 RBPs.
		Two versions of the PTBP1 pVTF map are provided: 
		(1) PTBP1_pVTF_map_no_affinity_cut-off.cys (Figure S6A): No PTBP1 binding affinity fold-change threshold was applied.
		(2) PTBP1_pVTF_map_affinity_FC_1.5.cys (Figure 6D): Variants shown have predicted PTBP1 binding affinity fold change ≥ 1.5, resulting in a clearer and more interpretable pVTF map.

**3.KEGG_cytoscape_file/**

	PTBP1_pEPRI_mutation_KEGG_network.cys: 
		Cytoscape network corresponding to the KEGG pathway analysis of genes affected by mutations within PTBP1-mediated pEPRIs.
		Large nodes represent enriched KEGG pathway terms, whereas small nodes represent genes (Entrez Gene IDs).

These files are provided as source data accompanying the manuscript for academic and research purposes.

This repository is still being updated.
