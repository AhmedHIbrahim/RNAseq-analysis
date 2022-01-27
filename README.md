# rnaseq-analysis

## RNAseq analysis

- **Background**
    
    > RBM10 (RNA Binding Motif Protein 10) is a Protein Coding gene. Its mutations are associated with various human diseases.[[ref](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2868995/)] The synthesized protein plays a role in regulating alternative splicing with hnRNA. In alternative splicing, different proteins can be translated from the same mRNA. Actually, RBM10 is highly conserved among mammals; the human RBM10 protein shares 96% and 97% aa sequence homology with that of mice and rats, respectively. [[ref](https://www.sciencedirect.com/science/article/pii/S0378111921000573)] In this project, **the goal** is to find the impact on Vascular Smooth Muscle Cells (VSMC) after treating them with Rbm10-associate adenovirus and controls. The VSMC are an important component of blood vessels; These cells are located in the medium part of a blood vessel [[ref](https://www.intechopen.com/books/muscle-cell-and-tissue-current-status-of-research-field/the-role-of-vascular-smooth-muscle-cells-in-the-physiology-and-pathophysiology-of-blood-vessels#:~:text=Vascular%20smooth%20muscle%20cells%20(VSMCs)%20are%20an%20important%20component%20of,lumen%20and%20form%20numerous%20layers.)]. RNA-binding motif protein 10 (RBM10) is reported to induce cell cycle arrest and apoptosis. Several studies have demonstrated that RBM10 exerts negative effects on tumor proliferation.[[ref](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA000712)]
    > 
    
    > The objective of choosing this data is to run different RNA-seq analyses on it and be able to interpret the output of each step in the pipeline and tweak the tools parameters in order to get better results.
    > 
- **Data Retrieval**
    - The original data were obtained from the Gene Expression Omnibus (GEO) repository. The URL of the data is: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188258](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188258)
    - The data consists of 6 different SRRs (3 Target, 3 Control) and it is paired-end. So, the downloaded files are 12 with a compressed size of 19G.
        
        ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image31.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image31.png)
        
    - The data was retrieved using **prefetch** software from The sratoolkit. The downloaded files were in **.sra** format. With the help of **fastq-dump** from the sratoolkit, the **fastq** datasets were obtained from the .sra datasets.
    
- **Initial Quality Check**
    
    ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image7.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image7.png)
    
    > The MultiQC Sequence Quality Histograms show high sequence quality mean scores.
    > 
    
    ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image15.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image15.png)
    
    > The MultiQC Adapter Content indicates that the sequences have adapter contamination that needs to be cleaned before progressing with the next steps of RNA analysis.
    > 
    
    ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image5.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image5.png)
    
    > In the Status Check diagram, as the top blue line (SRR16798080_2) shows an abnormality that crosses 5% in the Adapter Content Diagram, the cross between SRR16798080_2 and Adapter content in the following diagram is highlighted in yellow. Also, the following diagram shows an abnormality in the “Per base Sequence Content” and also in the “Sequence Duplication” in all of the provided fastq data.
    > 
    
    ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image2.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image2.png)
    
- **Adapter Trim**
    
    > In order to trim the adapters and clean the sequences, the TruSeq illumina adapter sequences were obtained from [here](https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/SingleIndexes.htm). The adapters were stored in a **.fa** file to be used later.
    > 
    
    **Read 1 ⇒** AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    
    **Read 2 ⇒** AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    
    > The **flexbar** tool was used to trim the illumina adapters from the original .fastq files.
    > 
- **Quality Check After Adapters Removal**
    
    ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image14.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image14.png)
    
- **Alignment with Ref Genome**
    
    > HISAT2 tool is used to align the trimmed paired-end reads file with the human reference genome (hg38.fasta) that have been retrieved from [here](https://www.ncbi.nlm.nih.gov/genome/guide/human/).
    > 
    
    [HISAT2 Alignment Result](https://www.notion.so/6bf203b481c74fdbae16418f4e06ac7f)
    
- **Expression**
    - **StringTie**
        - Stringtie in 'reference only' mode generates expression estimates from the SAM/BAM files generated by HISAT2, and uses the reference annotation file (.gtf) file that have been retrieved from [here](http://ftp.ensembl.org/pub/release-99-gtf/homo_spiens/Home_sapiens.GRCh38.99.gtf.gz) to guide the assembly process. In order to identify differentially expressed genes between experiments, StringTie's output can be processed by specialized software like *Ballgown*, Cuffdiff, or other programs (DESeq2, *edgeR*, etc.).
            
            > (Ad-Gfb_rep1)gene_abundances.tsv
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image3.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image3.png)
            > 
            
            > (Ad-Gfb_rep1)Transcripts.gtf
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image30.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image30.png)
            > 
        - “Create tidy expression matrix files for the StringTie results. This will be done at both the gene and transcript level and also will take into account the various expression measures produced: coverage, FPKM, and TPM.”
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image12.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image12.png)
            
            > gene_fpkm_all_samples.tsv
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image32.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image32.png)
            > 
        
        **Scripts Ref:**
        
        ./scripts/rnaseq/stringtie-expression.bash
        
        ./scripts/rnaseq/stringtie-expression-matrix.bash
        
    - **HTSEQ-COUNT**
        
        htseq-count produces raw gene counts instead of FPKM/TPM values for differential expression analysis from the created alignments.
        
        /expression/htseq_counts/*
        
        ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image21.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image21.png)
        
- **Differential-Expression**
    - **Ballgown DE Analysis**
        
        Perform Ad-RBM10 vs. Ad-Gfp comparison, using all replicates, for known (reference only mode) transcripts that are created using StringTie.
        
        Firstly run (“/scripts/pre-ballgown-setup.bash”) to create a file that lists our 6 expression files.
        
        - [Ballgown R Tutorial #1](https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_Part1_ballgown.R)
            - **Goals**
                - Run DE analysis without filtering the data.
                - Filter the data by ignoring low-abundance genes, then run DE analysis.
                - Get the significant genes with p-value < 0.05.
            - **Summary**
                1. How many genes are there? **60676**
                2. How many passed filter in Ad-Rbm10 or Ad-Gfp? **8717**
                3. How many differentially expressed genes were found (p-value < 0.05)? **1054**
            
            > de/ballgown/ref_only/**Ad-Rbm10_vs_Ad-Gfb_gene_results_sig.tsv**
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image19.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image19.png)
            > 
            
            > /de/ballgown/ref_only/**Ad-Rbm10_vs_Ad-Gfb_transcript_results_sig.tsv**
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image25.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image25.png)
            > 
    - **edgeR Analysis**
        
        [https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_edgeR.R](https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_edgeR.R)
        
        An alternative approach of finding the differently expressed genes using stringTie and Ballgown. edgeR uses the raw counts that have been created with htseq-counts tool, to do the differential expression of the count-based RNA-seq data.
        
        ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image24.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image24.png)
        
    - **Most significant genes In Ballgown**
        - **Higher abundance Genes in Ad-Rbm10**
            
            ```bash
            grep -v feature Ad-Rbm10_vs_Ad-Gfb_gene_results_sig.tsv | sort -rnk 3 | head -n 20 #Higher abundance in UHR
            ```
            
            [Genes List](https://www.notion.so/99b9017145814499897ac914e057aebe)
            
        - **Higher abundance Genes in Ad-Gfp**
            
            [Genes List](https://www.notion.so/72639f8be639434aa4f354ab0eb5a817)
            
        - **Based on FC and P-value**
            
            In /de/ballgown/ref_only/Ad-Rbm10_vs_Ad-Gfb_gene_results_sig.tsv, get the genes that have a fold-change greater than 10. The table shows the gene sorted manually by hight to low fc value.
            
            ```r
            **geneData <- read.table(file = 'Ad-Rbm10_vs_Ad-Gfb_gene_results_sig.tsv', sep = '\t', header = TRUE)
            topGenes <- geneData[geneData$fc > 10,]
            write.table(topGenes, file="top_ballgowns-gene.txt", sep="\t", row.names=FALSE, quote=FALSE)**
            ```
            
            [Genes List](https://www.notion.so/b906376a00354bb7ade2e4d3bc8e979f)
            
    - **Significant Genes in edgeR**
        
        ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image27.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image27.png)
        
    - **Venn Diagram**
        
        > The following Venn Diagram created using [venny tool](https://bioinfogp.cnb.csic.es/tools/venny/), shows the intersection between the significant genes that have been found using (stringTie-Ballgown) and (htseq-count-edgeR).
        > 
        
        ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image29.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image29.png)
        
        > The common 527 genes are listed here**:** [https://shrib.com/#intersection-between-ballgown-and-edgeR](https://shrib.com/#intersection-between-ballgown-and-edgeR)
        > 
        > 
        > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image4.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image4.png)
        > 
- **Ballgown Visualization**
    - Section 1 ([Tutorial_Part2_ballgown.R](https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_Part2_ballgown.R))
        - **[RBM10 Gene](https://www.genecards.org/cgi-bin/carddisp.pl?gene=RBM10)**
            
            Firstly, We have to know the following “A single gene can produce multiple different RNAs (i.e., transcripts).” Focusing on **RBM10** gene which presents in the high significant genes list. The following screenshot displays a few rows that were a result of the following commands.
            
            ```bash
            cd /expression/stringtie/ref_only/Ad-Rbm10_rep1
            cat transcripts.gtf | grep "RBM10" | less -S
            ```
            
            > It shows that **RBM10** gene (green_box and blue_box) has multiple transcripts (orange_box) and each transcript has multiple exons.
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image9.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image9.png)
            > 
            
            > So, we have to find the corresponding **transcript Ids** in the ballgown object that has been created in the following [tutorial](https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_Part1_ballgown.R).
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image8.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image8.png)
            > 
            
            > So, let's create a box plot of a single transcript to see how the expression differs between our control(Ad-Gfp) and our target Group(Ad-Rbm10) for ENST00000377604 Transcript which is 1 of 6 transcripts that RBM10 gene has.
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image6.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image6.png)
            > 
            
            > Time to see how all of the 6 transcripts differ in expression in all of the control/target replicates. The plot shows that (ENST00000377604) 221704 and (ENST00000329236) 221707 Transcripts of RBM10 genes are showing high expression in Ad-RBM10 replicates.
            > 
            
            ![Untitled](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/Untitled.png)
            
            ![Untitled](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/Untitled%201.png)
            
            ![Untitled](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/Untitled%202.png)
            
            ![Untitled](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/Untitled%203.png)
            
            ![Untitled](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/Untitled%204.png)
            
            ![Untitled](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/Untitled%205.png)
            
    - Section 2 ([Tutorial_Supplementary_R.R](https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_Supplementary_R.R))
        - **Heatmap Plot**
            
            Creating a heatmap to visualize expression differences between the genes that are differentially expressed under a provided p-value threshold. P.val <0.05 => 1054 genes, P.val < 0.005 => 119 genes, P.val < 0.001 => 26 genes.
            
            > (P-value < 0.05)
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image36.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image36.png)
            > 
            
            > (P-value < 0.005)
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image26.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image26.png)
            > 
            
            > (P-value < 0.001)
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image33.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image33.png)
            > 
        - **MDS Plot**
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image28.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image28.png)
            
        - **Other plots**
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image10.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image10.png)
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image1.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image1.png)
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image23.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image23.png)
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image13.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image13.png)
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image35.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image35.png)
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image16.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image16.png)
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image34.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image34.png)
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image18.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image18.png)
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image22.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image22.png)
            
- **Variation Analysis**
    - Introduction
        
        In order to detect any mutations/variants that have occurred, **rnaseqmut** has been used; it is a lightweight C++ program that detects variants (or mutations, including SNPs, indels) from RNA-Seq BAM files[[ref](https://github.com/davidliwei/rnaseqmut)]. If the alignment data is in sam format, **samtools** can be used to convert it to bam format(binary, lightweight version) and also enables the user to index it (creates .bam.bai file).
        
        > Sample Bash Script:
        > 
        > 
        > ```bash
        > samtools sort $RNA_HOME/alignments/hisat2/$i.sam > $i.bam
        > samtools index $i.bam
        > ```
        > 
    - Dataset
        
        > The dataset contains 3 Control Replicates (Ad-Gfb_rep*.bam) and three target replicates (Ad-Rbm10_rep*.bam). Following are the bam and indexed(.bai) files of the alignment dataset after conversion from (.sam) format.
        > 
        
        ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image14%201.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image14%201.png)
        
    - Rnaseqmut ([Full Script](https://github.com/davidliwei/rnaseqmut/blob/master/demo/rundemo.sh))
        - **de-novo mutation calling**
            
            > The target is to find the genetic alteration that is present for the first time in one family member as a result of a variant (or mutation) in a germ cell (egg or sperm) of one of the parents or a variant that arises in the fertilized egg itself during early embryogenesis. Also called de novo variant, new mutation, and new variant.[ref]
            > 
            > 
            > More explanation: [[https://www.youtube.com/watch?v=dzjUDPeIO24](https://www.youtube.com/watch?v=dzjUDPeIO24)]
            > 
            > Sample script (to be run on the 6 sample (3-Controls & 3-Targets):
            > 
            > ```bash
            > file=”Ad-Gfb_rep1.bam”
            > **rnaseqmut.linux.x64** $file > results/$file.1st.txt
            > ```
            > 
            > **Script Explanation:**
            > 
            > “scans the samples individually and gets the mutation list for each sample. In this study, it will output all possible mutations (those even occur in only 1 RNA-Seq read), but in reality, we may just need those with enough read support (controlled by -i/--min_read option).”
            > 
            > **Output files**:
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image7%201.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image7%201.png)
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image3%201.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image3%201.png)
            > 
        - **merging mutations in Step 1 into a candidate mutation list**
            
            > Merges the mutation lists from individual samples. These will be the potentially interesting mutations we would like to investigate.
            > 
            > 
            > ```bash
            > **merge1stfile** results/*.1st.txt > results/ALLMUTLIST.txt
            > ```
            > 
            > **Merge1stfile** is a bash script provided by the rnaseqmut tool, it will run the following command underhood. The explanation of the script can be found [here](https://explainshell.com/explain?cmd=sort+-k+1%2C1+-k+2%2C2n+-k+3%2C3+-k+4%2C4+-u++-m+%22%24%40%22+%7C+perl+-ane+%27print+if+%24F%5B3%5D+ne+%22N%22%27).
            > 
            > > sort -k 1,1 -k 2,2n -k 3,3 -k 4,4 -u -m "$@" | perl -ane 'print if $F[3] ne "N"'
            > > 
            
            Output file (11,076,840 records):
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image9%201.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image9%201.png)
            
        - **mutation calling from the merged lists**
            
            > Step 3, scan the samples again, using the provided mutation list in Step 2.
            > 
            > 
            > The following script should be run for the (Control & Target) replicates one by one.
            > 
            > ```bash
            > file="Ad-Gfb_rep1.bam"
            > rnaseqmut.linux.x64 -l results/ALLMUTLIST.txt $file > results/$file.2nd.txt
            > ```
            > 
            > What is -l param?
            > 
            > **-l** <mutation_list>, --mutation_list <mutation_list>
            > 
            > > The text file of a given, sorted list of mutations. Each line in a file records one mutation, with chromosome, location, reference and alternative sequence (separated by a tab). The output will only include mutations within a given mutation list.
            > > 
            
            Output files:
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image1%201.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image1%201.png)
            
        - **merge the second pass of mutations into a big table**
            
            Step 4, merge the mutations in step 3 into a big table. In this table, each row represents a mutation, and columns record the number of supporting reads (and the number of reference reads that do not support this mutation) in all samples. Based on this table, we can use our own criteria to search for interesting mutations.
            
            ```bash
            LABELS="Ad-Gfb_rep1,Ad-Gfb_rep2,Ad-Gfb_rep3,Ad-Rbm10_rep1,Ad-Rbm10_rep2,Ad-Rbm10_rep3"
            merge2ndvcf.py -l $LABELS results/*.2nd.txt > results/ALLMUT.txt
            ```
            
            > Comparison between number of records in merged files in step2 and step4
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image10%201.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image10%201.png)
            > 
            
            > ALLMUT.txt File Overview
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image4%201.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image4%201.png)
            > 
        - **filter mutations based on user-defined parameters**
            
            ```bash
            LABELS="Ad-Gfb_rep1,Ad-Gfb_rep2,Ad-Gfb_rep3,Ad-Rbm10_rep1,Ad-Rbm10_rep2,Ad-Rbm10_rep3"
            CONTROLGROUP="0,1,2"
            python script/filtermut.py -d **5** -f 0.0 -b 0 -c $CONTROLGROUP -l $LABELS < results/ALLMUT.txt > results/ALLMUT_FILTERED.vcf
            ```
            
            > Step 5 illustrates an example of filtering interesting mutations using a python script "filtermut.py". In this study, we define "control" samples as two normal samples, and would like to look for mutations that satisfy the ALL of the following criteria:
            > 
            > - with at least 0% frequency (controlled by -f/--min-recfrac option)
            > - do not have any supporting reads in control samples (-b/--max-alt)
            > - and 5 supporting reads (-d/--min-recread);
            
            > The final output is a VCF file recording mutations and read coverages in all 4 samples. You can also output tab-delimited file instead of VCF file for further downstream analysis (use -z/--no-vcf option in Step 5). For interpreting results, see the next section.
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image12%201.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image12%201.png)
            > 
            
            > ALLMUT_FILTERED.vcf File {[Interpretation](https://github.com/davidliwei/rnaseqmut#interpreting-results)}
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image2%201.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image2%201.png)
            > 
    - Annovar ([Full Flow](http://eng1.mu.edu.tr/~tugba/SeqAnalysis/variation.pipeline))
        
        > The goal is to annotate the mutation that are listed under ALLMUT_FILTERED.vcf file using Annovar.
        > 
        - Download a few annotation databases
            
            ```bash
            CMD="annotate_variation.pl -buildver hg38 -downdb"
            
            perl $CMD -webfrom annovar refGene humandb/
            perl $CMD cytoBand humandb/
            perl $CMD genomicSuperDups humandb/
            perl $CMD -webfrom annovar esp6500siv2_all humandb/
            perl $CMD -webfrom annovar 1000g2015aug humandb/
            perl $CMD -webfrom annovar exac03 humandb/
            perl $CMD -webfrom annovar avsnp150 humandb/
            perl $CMD -webfrom annovar dbnsfp30a humandb/
            perl $CMD -webfrom annovar clinvar_20200316 humandb/
            perl $CMD -webfrom annovar cosmic70 humandb/
            ```
            
        - Add a filter column.
            
            ```bash
            IN_FILE="ALLMUT_FILTERED.vcf"
            OUT_FILE="ALLMUT_FILTERED.filtercoladded.vcf"
            awk '{FS="\t";print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t.\t"$7}' $IN_FILE | tail -13547 > $OUT_FILE
            ```
            
            > ALLMUT_FILTERED.filtercoladded.vcf
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image13%201.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image13%201.png)
            > 
        - Run annovar on the variants we found to check if they exist in any of the databases.
            
            ```bash
            RNASEQMUT_PATH="~/workspace/seq-analysis-course/rnaseqmut/demo/results"
            DATABASES="refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eur,exac03,avsnp150,dbnsfp30a,cosmic70,clinvar_20200316"
            OPERATIONS="g,r,r,f,f,f,f,f,f,f,f"
            
            perl table_annovar.pl $RNASEQMUT_PATH/ALLMUT_FILTERED.filtercoladded.vcf humandb/ \
            - buildver hg38 -out myanno -remove -protocol $DATABASES \
            - operation $OPERATIONS -nastring . -vcfinput
            ```
            
            [Annovar Operations](https://www.notion.so/4c03b6043cac4c91a41570d36e352d8e)
            
        - Output List
            
            > [help] Somatic Mutation Annotators Output File Interpretation
            > 
            > 
            > [https://brb.nci.nih.gov/seqtools/colexpanno.html](https://brb.nci.nih.gov/seqtools/colexpanno.html)
            > 
            > [https://brb.nci.nih.gov/seqtools/colexpanno.html#dbnsfp](https://brb.nci.nih.gov/seqtools/colexpanno.html#dbnsfp)
            > 
            
            > **Results myanno.hg38_multianno.vcf file:**
            > 
            > 
            > [https://docs.google.com/spreadsheets/d/1kZ1IY0Y_uFAiCycVeo9LXXHd4g9uNKAl/edit?usp=sharing&ouid=103975271912997279441&rtpof=true&sd=true](https://docs.google.com/spreadsheets/d/1kZ1IY0Y_uFAiCycVeo9LXXHd4g9uNKAl/edit?usp=sharing&ouid=103975271912997279441&rtpof=true&sd=true)
            > 
            
            > Look-up dbSNP database with **avsnp150** id. An example is [rs992574219](https://www.ncbi.nlm.nih.gov/snp/?term=rs992574219)
            > 
        - Filtering Annovar Output
            1. Select only records Func.refGene = “exonic
            2. Select only records with sift score = 0
            3. Select only records with Polyphen2_HDIV_score = 1
            4. Select only records with Polyphen2_HVAR_score = 1
            5. Remove records with 1000g2015aug_all and 1000g2015aug_eur values if exist
            
            > After the filtering, only **141** unique genes remain with 326 variants.
            > 
        - Venn Diagram
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image11.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image11.png)
            
            The only common gene found between the two datasets is **RBM10**
            
        - COSMIC related findings
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image5%201.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image5%201.png)
            
            > [https://cancer.sanger.ac.uk/cosmic/mutation/overview?id=106076077](https://cancer.sanger.ac.uk/cosmic/mutation/overview?id=106076077)
            > 
            > 
            > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image8%201.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image8%201.png)
            > 
    
- **Enrichment Analysis**
    - **Introduction**
        
        > **Gene set enrichment analysis (GSEA)** (also **functional enrichment analysis**)
        > 
        > 
        > is a method to identify classes of [genes](https://en.wikipedia.org/wiki/Genes) or [proteins](https://en.wikipedia.org/wiki/Proteins) that are over-represented in a large set of genes or proteins, and may have an association with disease [phenotypes](https://en.wikipedia.org/wiki/Phenotype).[[ref](https://en.wikipedia.org/wiki/Gene_set_enrichment_analysis#:~:text=Gene%20set%20enrichment%20analysis%20(GSEA,an%20association%20with%20disease%20phenotypes.)]
        > 
        
        > **Over Representation Analysis (ORA)**
        > 
        > 
        > is a widely used approach to determine whether known biological functions or processes are over-represented (= enriched) in an experimentally-derived gene list. [[ref](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html)]
        > 
        
        > **The Gene Ontology (GO)**
        > 
        > 
        > provides a system for hierarchically classifying genes or gene products into terms organized in a graph structure (or an ontology). The terms are groups into three categories: **molecular function** (describing the molecular activity of a gene), **biological process** (describing the larger cellular or physiological role carried out by the gene, coordinated with other genes), and **cellular component** (describing the location in the cell where the gene product executes its function). [[ref](https://en.wikipedia.org/wiki/Gene_Ontology_Term_Enrichment)]
        > 
        
        > **KEGG: Kyoto Encyclopedia of Genes and Genomes**
        > 
        > 
        > KEGG is a database resource for understanding high-level functions and utilities of the biological system, such as the cell, the organism, and the ecosystem, from molecular-level information, especially large-scale molecular datasets generated by genome sequencing and other high-throughput experimental technologies.[[ref](https://www.genome.jp/kegg/)]
        > 
        
        > **C7: immunologic signature gene sets**
        > 
        > 
        > It contains Gene sets that represent cell states and perturbations within the immune system. [[ref](http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)]
        > 
        
        > **C6: oncogenic signature gene sets**
        > 
        > 
        > It represents signatures of cellular pathways that are often dis-regulated in cancer. The majority of signatures were generated directly from microarray data from NCBI GEO or from internal unpublished profiling experiments involving the perturbation of known cancer genes. [[ref](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C6)]
        > 
    - **Dataset**
        - **DE Significant Genes with EntrezIds**
            
            ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image8%202.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image8%202.png)
            
        - **Variant Analysis Filtered genes**
            
            [https://docs.google.com/spreadsheets/d/1fOZXgFKTsf6GB_m4rm_ovSw57JiwkRCA/edit?usp=sharing&ouid=103975271912997279441&rtpof=true&sd=true](https://docs.google.com/spreadsheets/d/1fOZXgFKTsf6GB_m4rm_ovSw57JiwkRCA/edit?usp=sharing&ouid=103975271912997279441&rtpof=true&sd=true)
            
    - **Over Representation Analysis (ORA)**
        
        Script Ref: over-representation-analysis.R
        
        - **Based on Ballgowns Significant Genes**
            - **Go (Gene Ontology)**
                
                ```r
                ## GO over-representation analysis
                ego <- enrichGO(gene = gene.df$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
                
                ego_df = data.frame(ego)
                ```
                
                **# NO Results Founds**
                
            - **KEGG (Kyoto Encyclopedia of Genes and Genomes)**
                
                ```r
                keggEnrich<-enrichKEGG(gene= gene.df$ENTREZID,
                organism= "hsa",
                pAdjustMethod="BH",
                pvalueCutoff = 0.05)
                keggDf<-data.frame(keggEnrich)
                ```
                
                **# NO Results**
                
            - **C7:** **Immunologic signature gene sets**
                
                clusterProfiler package in R provides an enricher function for hypergeometric test and GSEA function for gene set enrichment analysis that is designed to accept user-defined annotation. [[ref](http://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html)].
                
                ```r
                c7 <- read.gmt("c7.all.v7.1.entrez.gmt")
                
                ImmunSigEnrichC7 <- enricher(gene.df$ENTREZID,
                TERM2GENE=c7,
                pvalueCutoff = 0.05)
                ImmunSigEnrichC7 <- setReadable(ImmunSigEnrichC7,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID")
                
                write.csv(data.frame(ImmunSigEnrichC7), "ballgown-based-ora-immune.csv")
                ```
                
                > Result
                > 
                > 
                > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image5%202.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image5%202.png)
                > 
                
                > 1)
                > 
                > 
                > [GSE11864_UNTREATED_VS_CSF1_IFNG_IN_MAC_UP](https://www.notion.so/7eeab3ff2fcb4e70902d31838e2df71a)
                > 
                > > **Probable Relation:**
                > > 
                > > 
                > > Since in several vascular diseases abnormal vascular smooth-muscle cell (VSMC) proliferation is often associated with the presence of macrophages, we examined whether macrophage-colony-stimulating factor (M-CSF) might play a role in the control of VSMC growth. ...... It was further demonstrated that M-CSF could act in synergy with thrombin, platelet-derived growth factor or basic fibroblast growth factor in promoting VSMC DNA synthesis. These results support the hypothesis that M-CSF affects the growth of cultured **rat** VSMCs through paracrine/autocrine mechanisms. **Its effects at both the macrophage and the VSMC level confer to M-CSF a central role in the development of vascular lesions that occurs during atherosclerotic progression.** [[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1218536/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1218536/)]
                > > 
                
                > 2)
                > 
                > 
                > [GSE13484_12H_VS_3H_YF17D_VACCINE_STIM_PBMC_DN](https://www.notion.so/08d7c5f0af464b6295c672914db45163)
                > 
                > > **Probable Relation:**
                > > 
                > > 
                > > The data of this study suggest that expression alterations in PBMCs of microRNAs, which are involved in the regulation of VSMCs phenotype, possibly reflect the disease status in hypertensive patients. [https://www.nature.com/articles/jhh2013117]
                > > 
            - **C6: oncogenic signature gene sets**
                
                ```r
                c6 <- read.gmt("c6.all.v7.1.entrez.gmt")
                
                oncogSigEnricC6 <- enricher(gene.df$ENTREZID,
                TERM2GENE=c6,
                pvalueCutoff = 0.05)
                oncogSigEnrichC6 <- setReadable(oncogSigEnricC6,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID")
                write.csv(data.frame(oncogSigEnrichC6), "ballgown-based-ora-oncog.csv")
                ```
                
                **# NO Results**
                
        - **Based on Annovars Genes**
            - **Go (Gene Ontology)**
                
                ```r
                ## GO over-representation analysis
                
                ego <- enrichGO(gene = gene.df$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.01,
                readable = TRUE)
                
                ego_df = data.frame(ego)
                write.csv(ego_df, "annovar-based-ora-go-p-0.01.csv")
                ```
                
                ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image3%202.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image3%202.png)
                
                > Full 55 records:
                > 
                > 
                > [https://docs.google.com/spreadsheets/d/1tpJrKTPFwJgxXC9g7dCW0Fjln3veXg8ILqeKXNKZ4eE/edit?usp=sharing](https://docs.google.com/spreadsheets/d/1tpJrKTPFwJgxXC9g7dCW0Fjln3veXg8ILqeKXNKZ4eE/edit?usp=sharing)
                > 
                > > To explore the found results:
                > > 
                > > 
                > > https://www.ebi.ac.uk/QuickGO/term/GO:******
                > > 
            - **KEGG (Kyoto Encyclopedia of Genes and Genomes)**
                
                ```r
                keggEnrich<-enrichKEGG(gene= gene.df$ENTREZID,
                organism= "hsa",
                pAdjustMethod="BH",
                pvalueCutoff = 0.01)
                
                keggDf<-data.frame(keggEnrich)
                ```
                
                ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image6%201.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image6%201.png)
                
                > Full List: [https://drive.google.com/file/d/1uLV7l7fBTbIDvZLVOy7o18-0PHamtktn/view?usp=sharing](https://drive.google.com/file/d/1uLV7l7fBTbIDvZLVOy7o18-0PHamtktn/view?usp=sharing)
                > 
                > 
                > > To explore the results e.g.:
                > > 
                > > 
                > > https://www.kegg.jp/kegg-bin/show_pathway?hsa05165/5728/3106/9636/998/3688/8325/7132/3678/1282/1277/2335/1293/4599/5894/1499/3912/8323/1288
                > > 
                > > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image2%202.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image2%202.png)
                > > 
                
            - **C7:** **Immunologic signature gene sets**
                
                ```r
                c7 <- read.gmt("c7.all.v7.1.entrez.gmt")
                
                ImmunSigEnrichC7 <- enricher(gene.df$ENTREZID,
                TERM2GENE=c7,
                pvalueCutoff = 0.05)
                ImmunSigEnrichC7 <- setReadable(ImmunSigEnrichC7,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID")
                
                write.csv(data.frame(ImmunSigEnrichC7), "annovar-based-ora-immune-p-0.05.csv")
                ```
                
                **# NO Results**
                
            - **C6: oncogenic signature gene sets**
                
                ```r
                c6 <- read.gmt("c6.all.v7.1.entrez.gmt")
                
                oncogSigEnricC6 <- enricher(gene.df$ENTREZID,
                TERM2GENE=c6,
                pvalueCutoff = 0.05)
                oncogSigEnrichC6 <- setReadable(oncogSigEnricC6,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID")
                
                write.csv(data.frame(oncogSigEnrichC6), "annovar-based-ora-oncog-p-0.05.csv")
                ```
                
                ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image1%202.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image1%202.png)
                
                > [https://www.gsea-msigdb.org/gsea/msigdb/cards/ESC_V6.5_UP_EARLY.V1_DN](https://www.gsea-msigdb.org/gsea/msigdb/cards/ESC_V6.5_UP_EARLY.V1_DN)
                > 
                > 
                > ![rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image7%202.png](rnaseq-analysis%208db9d84b7dc849f193141b9945468838/image7%202.png)
                >