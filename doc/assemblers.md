# Assemblers

List of assemblers which might be useful for serratus.

Review article:

https://doi.org/10.1186/s40168-019-0626-5

This table lists 13 virus assemblers with links to code & papers:

https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0626-5/tables/1

Please update this list if you have ideas, corrections, comments. If you don't have commit rights to this repository, add a comment to Issue #71 or send me an email (robert@drive5.com). For each assembler, provide:

1. Name
2. Type (e.g. reference or de-novo)
3. Link to code
4. Link to paper
5. Comments on pros or cons for the serratus project.

Use "??" as a placeholder if not known.

**Kollector**  
Type: Targeted De-novo  
Code: https://github.com/bcgsc/kollector  
Paper: doi: 10.1093/bioinformatics/btx078  
Comments: Orignally tested for genomic DNA assembly, may need refinement to work with transcriptome data.  

**ABySS**  
Type: De-novo Genomic  
Code: https://github.com/bcgsc/abyss  
Paper: doi: 10.1101/gr.214346.116  
Comments: ??  

**RNA-Bloom**  
Type: De-novo Transcriptomic  
Code: https://github.com/bcgsc/RNA-Bloom  
Paper: doi: https://doi.org/10.1101/701607  
Comments: Isoform assembly maybe an unnecessary feature, but our datasets are expected to be transciptomic.

**SPAdes**
Type: De-novo genomic / transcriptomic / metagenomic (different varieties exist - rnaSPAdes, SPAdes meta etc.)
Code: https://github.com/ablab/spades
Paper: doi: 10.1089/cmb.2012.0021
Comments: Well-supported and generally robust assembler. SPAdes meta was highlighted in the review article at the top of the document ("Choice of assembly software has a critical impact on virome characterisation") as performing "consistently well".

**Megahit**
Type: De-novo genomic / metagenomic
Code: https://github.com/voutcn/megahit
Paper: doi: 10.1093/bioinformatics/btv033
Comments: Very memory-efficient.

**IDBA**
Type: De-novo metagenomic
Code: https://github.com/loneknightpy/idba
Paper: doi: 10.1093/bioinformatics/bts174
Comments: Anecdotally (i.e. in my own experience) works well for viral genome assembly. Also positively reviewed in the review paper above.
