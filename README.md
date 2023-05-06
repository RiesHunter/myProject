# MUSCLE is faster and just as accurate as ClustalW in aligning genetically homogenous influenza outbreaks
### Hunter J. Ries
### Department of Pathobiological Sciences, University of Wisconsin–Madison, Madison, Wisconsin, USA

## Abstract
Deep-sequencing viral genomics has provided high-resolution glimpses into respiratory virus evolution and epidemiology. While this has formed the foundation for modern public health infectious disease surveillance, it proves troublesome for resolving phylogenies in low-diversity outbreaks. During seasonal epidemics of influenza virus, acute infections drive transmission but not evolution within communities. Viruses transmitted from one host to another often share near-identical sequences, suggesting the virus is transmitted before it can accumulate mutations or adapt to its host. Thus, many community sequences are genetically similar but not closely related in their host histories. MUSCLE and ClustalW are standard multiple-sequence alignment tools used to align influenza virus sequences; however, the impact of alignment software on maximum likelihood trees needs to be clarified. In this study, I compare the maximum likelihood trees from ClustalW- and MUSCLE-aligned influenza A virus sequences from three separate, relatively genetically homogenous outbreaks. While MUSCLE is traditionally used for large or highly divergent populations and ClustalW for smaller or similar populations, both demonstrated similar tree topology, length, and likelihood. However, MUSCLE outpaced ClustalW and did not do so at the cost of accuracy. This study positions MUSCLE as rapid and accurate in our datasets, especially those with many samples.

Directory tree created with `tree`:
> ├── README.md
> ├── data
> │   ├── HK
> │   ├── cali09
> │   └── perth
> ├── figures
> │   ├── Fig1-All_ggtree_lowqual.png
> │   ├── Fig2-All_cophylo_lowqual.png
> │   ├── Fig3-run_stats.pdf
> │   ├── Fig3-run_stats_lowqual.png
> │   ├── cophylo
> │   └── ggtree
> ├── notebook-log.md
> ├── notes
> │   └── ClassNotes.md
> ├── results
> │   ├── HK
> │   ├── cali09
> │   ├── perth
> │   └── run_data.tsv
> ├── scripts
> │   ├── MLtree_root_score.R
> │   ├── run_clustalw.sh
> │   ├── run_muscle.sh
> │   └── split_concatenate.sh
> └── slides
>     ├── 230502-HJR_Presentation
>     └── 230502-HJR_Presentation.key