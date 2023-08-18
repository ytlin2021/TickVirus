# Characterisation and zoonotic risk of tick viruses in public datasets 
Code and data in the study in _Lin and Pascall. Characterisation and zoonotic risk of tick viruses in public datasets_. The study re-examined currently known tick-borne viruses and identified putative novel viruses associated with ticks in public datasets. 

## File structure
```
└─TickVirus/
   ├─BioinformaticAnalyses/ ................... Genbank accession numbers of the candidate sequences 
   │   ├─SRA_Accession.txt .................... Accession numbers of 19,990 SRA sequences after 
   │   │                                        removing duplicates 
   │   ├─TSA_Accession.txt .................... Accession numbers of 1328 TSA sequences after 
   │   │                                        removing duplicates 
   │   └─WGS_Accession.txt .................... Accession numbers of 1861 WGS sequences after 
   │                                            removing duplicates 
   ├─Phylogenetics/ ........................... Data used in phylogenetic analyses
   │   ├─Alignments/ .......................... Multiple sequence alignments using ClustalW 
   │   │  ├─AlphatetraviridaeMSA.fa ........... MSA of Alphatetraviridae RdRp protein sequences after 
   │   │  │                                     removing poorly aligned sequences and regions
   │   │  ├─ChuviridaeMSA.fa .................. MSA of Chuviridae RdRp protein sequences after 
   │   │  │                                     removing poorly aligned sequences and regions
   │   │  └─OrthomyxoviridaeMSA.fa ............ MSA of Orthomyxoviridae RdRp protein sequences after 
   │   │
   │   ├─ProteinSeqs/ ..........................Raw protein sequences used for alignments 
   │   │  ├─HepeAlphatetraviridae+candidates_RdRp.fa 
   │   │  ├─Chuviridae+TSAcandidate_RdRp.fa 
   │   │  └─Orthomyxoviridae+candidates_PB1.fa
   │   │  
   │   └─Trees/ ............................... Tree files in nexus format
   │      ├─AlphaMCC.tree ..................... The maximum clade credibility tree of Alphatetraviridae
   │      │                                     based on RdRp protein sequences
   │      ├─ChuMCC.tree ....................... The maximum clade credibility tree of Chuviridae
   │      │                                     based on RdRp protein sequences
   │      ├─OrthoMCC.tree ..................... The maximum clade credibility tree of Orthomyxoviridae
   │      │                                     based on RdRp protein sequences
   │      ├─AlphaMCMCMC.xml ................... BEAST2 XML file for the phylogenetic analysis of Alphatetraviridae
   │      ├─ChuMCMCMC.xml ..................... BEAST2 XML file for the phylogenetic analysis of Chuviridae                                
   │      └─OrthMCMCMC.xml .................... BEAST2 XML file for the phylogenetic analysis of Orthomyxoviridae
   │                                            
   └─ZoonoticRiskAnalyses/ .................... Data and scripts used in zoonotic risk analyses
       ├─MakeFigure4.R ........................ Scripts to generate Figure 4 
       ├─TBVCompleteGenome.csv ................ Taxonomy information of tick viruses with complete genomes
       ├─TBVCompleteGenome.predictions.csv .... Zoonotic risk prediction of tick viruses with complete genomes
       └─TickVector.csv ....................... Main tick vectors of tick viruses with complete genomes                                            
```

