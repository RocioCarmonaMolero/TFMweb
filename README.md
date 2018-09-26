# Master's Final Thesis: 'Identification and Characterization of the *Cicer arietinum* aldehyde dehydrogenase superfamily'


In this work we will see the next points:

* [Identification of the *C. arietinum* aldehyde dehydrogenase suprefamily](https://github.com/RocioCarmonaMolero/TFMweb#identification-of-the-cicer-arietinum-aldehyde-dehydrogenase-superfamily)
  * [Search for candidate proteins](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/README.md#search-for-candidate-proteins)
  * [Checking ALDH conserved domains](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/README.md#checking-aldh-conserved-domains)
  * [BLASTP of 'aldhcicer' against *C. arietinum* genome](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/README.md#blastp-of-aldhcicer-against-carietinum-genome)
* [Characterization of aldehyde dehydrogenases of *C. arietinum*](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/README.md#characterization-of-aldehyde-dehydrogenases-of-c-arietinum)
* [Classification of *C. arietinum* aldehyde dehydrogenases](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/README.md#classification-of-c-arietinum-aldehyde-dehydrogenases)
* [Phylogeny](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/README.md#phylogeny)
* [Duplication Analysis](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/README.md#duplication-Analisis)
* [Expression in silico](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/README.md#expression-in-silico)
* [Code availability](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/README.md#code-availability)


## Identification of the *Cicer arietinum* aldehyde dehydrogenase superfamily

### Search for candidate proteins

First, a keyword-based search was carried out in the 'Protein' database of the National Biotechnology Information Center ([NCBI](https://www.ncbi.nlm.nih.gov/)) specifying *Arabidopsis* as organism and 'Aldehyde dehydrogenase' as title. 136 proteins were obtained and their sequences were downloaded into a text file. A BLASTP (Altschul et al. 1990) was run against *Cicer arietinum* proteome in the Reference Protein database (refseq_protein). The thresholds established to select candidate aldolases from *C. arietinum* were: Query cover ≥ 25%, E-value ≥ e-25, Identity ≥ 25%. The chickpea sequences that fulfilled these parameters were downloaded into a csv file.

As each BLASTP was done indepedently, the existence of the same XP *Cicer* aldolase was likely in several different csv files. For this, an R code script was written to clean those files so that a vector called 'candidates' would be created. In this 'candidates' vector the accessions of the *Cicer* proteins of the 136 csv were included without being repeated. It was possible executing the function [get_IDprot()](https://github.com/RocioCarmonaMolero/ScriptProteinas/blob/master/get_IDprot.R).
The results of applying the function get_IDprot () are added to the vector 'candidates'; avoiding the accession repetition using the following commands:

> dat <-  read.csv( 1_2cicer_D9V2GJ7401R-Alignment-HitTable.csv , header = FALSE, stringsAsFactors = F)
 
> dat <- dat[,2]
 
> IDprot <- get_IDprot(dat, 2)
 
> candidatos <- c(candidatos, IDprot[!IDprot %in% candidatos])


First of all, the second column of the csv, downloaded directly from the BLASTP result, is taken; then, the function is applied and the resulting vector is saved in 'IDprot'. The accessions of 'IDprot' not found in 'candidates' are included in it. After executing it for the 136 csv, we have all the proteins of *C.arietinum* in which later we will check the aldolase domains.
The same process was repeated twice, being *Medicago truncatula* and *Glycine max* (soybean) the query in each of them. There were 36 aldehyde dehydrogenases found in [Phytozome v12.1](https://phytozome.jgi.doe.gov/pz/portal.html) and 21 in Refseq_Protein for *Medicago* and 55 those found for soybean in the NCBI.
The final 'candidate' vector contains 43 proteins.

Considering that the search was carried out against the model plant *Arabidopsis thaliana*, and against two legumes: *Glycine max* and *Medicago truncatula*, which is a model in turn of this family; we can conclude to have covered an adequate spectrum of species that, by sequence homology, will allow us to identify all the aldehyde dehydrogenases of *Cicer arietinum*.


### Check ALDH conserved domains

To verify which of our candidate proteins are aldolases, the presence of the ALDH-superfamily conserved domains was checked: 'PF00171.21', 'PF07368.10' and 'PF05893.13' in Pfam; 'PS00687' and 'PS00070' in ScanProsite; the accession 'cl11961' in the database of Conserved Domains of the NCBI; and the accession '53720' in the Superfamily database.
After this screen 36 proteins are selected, which will constitute a new vector called 'aldhcicer'.

### BLASTP 'aldhcicer' against *Ca* genome 

Once the aldehyde dehydrogenases of chickpea are identified, a BLASTP is performed against its own genome. In this way, it will be possible to include any protein that does not have high homology with the other species, since in this case the query proteins are from the same organism.
One more protein was identified with the Pfam domain 'PF00171.21'. Thus, we conclude that the chickpea ALDH superfamily has 37 members.


## Characterization of CaALDH
 
For the characterization of the aldolases, attention was focused on obtaining the following properties:
* Number of amino acids (naa)
* Accession of your mRNA (accRNA)
* Length in bp of its RNA (RNAlen)
* Name of the gene locus (LOC)
* Chromosome in which it is (chr)
* Number of exons (exon)
* Position of the start of the protein (startp)
* Position of the end of the protein (stop)
* Molecular weight of the protein (mol_wg)

In order to collect this information, the function [caracterization()](https://github.com/RocioCarmonaMolero/ScriptProteinas/blob/master/InformationProteins.R) is applied; naming the function [weight_prot()](https://github.com/RocioCarmonaMolero/ScriptProteinas/blob/master/get_mol_weight.R).

Open access bioinformatic tools were also used to predict the subcellular location and the active center. This knowledges will also be helpful in the next step: the classification of proteins.
For cell localization, [DeepLoc 1.0](http://www.cbs.dtu.dk/services/DeepLoc/), which differentiates between 10 eukaryotic protein locations: Cytoplasm, Extracellular, Mitochondria, Cell Membrane, Chloroplast, Golgi apparatus, Lysosome/Vacuole, Nucleus, Endoplasmic reticulum and Peroxisome; [SLP-Local](http://sunflower.kuicr.kyoto-u.ac.jp/~smatsuda/slplocal.html) and [TargetP 1.1](http://www.cbs.dtu.dk/services/TargetP/) were used. Moreover, DeepLoc 1.0 indicates whether it is a soluble or a membrane protein. In the membrane protein results, the transmembrane domain was verified by [Smart](http://smart.embl-heidelberg.de/)(Simple Modular Architecture Research Tool). In order to consolidate the chloroplast and mitochondria location results, the presence of chloroplast transit peptides (cTP) was confirmed by [ChloroP 1.1](http://www.cbs.dtu.dk/services/ChloroP/) and mitochondrial-targeting sequences (mTP) by [Mitoprot](https://ihg.gsf.de/ihg/mitoprot.html). 
For the active center, the protein sequences were analyzed with [PROSITE](https://prosite.expasy.org/), looking up the identification of the active center PS00687 (glutamic acid active site, E) and PS00070 (cysteine active site, C). In case another active site was found, it was also annotated.

[PROPSEARCH](http://abcis.cbs.cnrs.fr/propsearch/) was used to determine the molecular function. Each protein was analyzed individually and the hits with highest probability were written down. For the verification of these, it was checked whether the functional residues were conserved between query and hit. First through a database Blast search and after through an alignment between them with [Clustal Omega](http://www.clustal.org/omega/) program.

![Scheme of the bioinformatic tools used for the characterization step](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/function%20scheme.jpg)
**Figure 1**.Scheme of the bioinformatic tools used for the characterization and its aim


## Classification of CaALDHs

For the classification, the criteria stablished by the ALDH Gene Nomenclature Committee (AGNC) in 1999 was applied. Two proteins belong to the same gene family if they have more than 40% identity; to the same subfamily if they have more than 60% identity. In the annotation, the root 'ALDH' is followed by a family descriptor number, a capital letter to describe the subfamily, a number specifying the individual gene within the subfamily and a lowercase letter if necessary to designate variants (Zhu et al., 2014).
Frequent methods for classifying this family in plants are based on homology with other plant species that are already described. For this, BLASTP of all 37 chickpea aldolases was made against the Protein Refseq database of *Medicago truncatula* and *Glycine max*, both being legumes and *Medicago* the model plant of them. All the results were downloaded and filtered to eliminate those whose identity was <40% and whose query length was less than the hit length: [Blast_Sieve.R](https://raw.githubusercontent.com/RocioCarmonaMolero/ScriptProteinas/master/Blast_Sieve.R)


## Phylogeny

A multiple alignment of the sequences of the aldehyde dehydrogenases of *Medicago truncatula* (MtALDH) and *Cicer arietinum* (CaALDH) was performed in the form of one gene per locus with the MUSCLE program using the default parameters (Edgar, 2004). To deduce the evolutionary history, the Neighbor-Joining method was used (Saitou & Nei, 1987). The consensus bootstrap tree was inferred from 1000 repetitions (Felsenstein, 1985). Evolutionary distances were calculated using the Poisson correction method (Zuckerkandl & Pauling, 1965), being in the units of the number of amino acid substitutions per site. Evolutionary analyzes were performed in MEGA6 (Tamura et al., 2013). We consider sister pairs those proteins grouped on the basis of high bootstrap values (> 65%) [(Die et al., 2018)](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4695-9).


## Duplication Analysis

For the duplication analysis of the CaALDH, the Circoletto program (Darzentas, 2010) was used. Both the FASTA query and the FASTA database are the protein sequences of CaALDH; with ultra-strict E-value values (10 -180), using absolute score / colored bands and colors: green for identities ≤ 95%, orange ≤ 99% and red> 99%.


## Expression in silico

The coding sequences of the CaALDH genes were used as a query against the chickpea NCBI EST database. The search parameters were established as follows: megablast, identity> 90%, length> 180 bp and E-value <10-10 [(Die et al., 2018)](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4695-9).


## Code availability

The codes of this TFM have been written in the R programming language, with the R software (Team R.C., 2017) and the free access RStudio interface (Team R, 2016, http://www.rstudio.com/). With these scripts we have collected data from the NCBI and performed the analysis of them. They are available in the repository: https://github.com/RocioCarmonaMolero/ScriptProteinas. The code is distributed under the MIT open source license.

Every step is deeply explained in the [TFM](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/TFM.md) file.



# Results


## Identification and characterization of aldehyde dehydrogenases of C. arietinum

The identification process of the CaALDH is shown in Figure 1. Considering the species of the search, the model plant *Arabidopsis thaliana* and the legumes *Glycine max* and *Medicago truncatula*, we can conclude to have covered an adequate spectrum of species that, by homology of sequences, will allow us to identify all the aldehyde dehydrogenases of *Cicer arietinum*.

The availability of the complete chickpea genome together with database searches allowed the identification of 37 chickpea ALDH proteins encoded by 29 *CaALDH* genes (Table 1). These genes encode proteins with a range of 134 (CaALDH3H1j) to 755 (CaALDH18B3d) aa of length. The number of exons of the *CaALDH* genes varies from 3 (CaALDH3H1j) to 21 (CaALDH18B3a, CaALDH18B3d and CaALDH18B3e). The corresponding molecular weights range between 15.07 and 81.90 kDa; and the predicted isoelectric points between 4.34 and 9.49 (Table 2). The wide range of pI suggests that chickpea ALDH proteins work in very different subcellular environments.

We observed that families 5, 6, 11, 12 and 22 are defined by a single gene in chickpea, similar to *Arabidopsis, O. sativa, S. italica, S. bicolor, E. parvulum* and *E. salsugineum* (Table 3); suggesting that these families constitute ALDH house-keeping genes, involved in the plants central metabolism and the preservation of non-toxic aldehydes levels. Families 2, 3 and 18 are comparatively abundant in *C. arietinum*. Family 3 (10 genes) and family 18 (6 genes) have a higher number of members in chickpea than any plant species described so far. ALDH3 proteins constitute one of the most widespread and diverse groups of ALDH genes in plant species (Brocker et al., 2013).
There are no *CaALDH* genes in families 19, 21, 23 and 24. Families 21 and 23 have been found to contain only primitive terrestrial plant genes (Chen et al., 2002), and ALDH24 appears only in the single-celled alga *C. reinhardtii* (Wood & Duff, 2009); suggesting that these three families may have played an important role in the evolution of lower plants and subsequently were lost in higher plants. Family 19 is found only in tomato (*S. lycopersicum*), so it evolved specifically in this lineage (Jimenez-Lopez et al., 2016).

![](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/Esquema_CaALDH1.jpg)
**Figure 2**. CaALDH identification scheme. Own functions appear in color. Proteins from *Arabidopsis, G. max* and *Medicago* were used as a query in BLASTP searches against of the chickpea reference genome. The result of this was cleaned with R scripts to avoid repetitions. Chickpea proteins without the ALDH domains were deleted and another BLASTP was run against the *Cicer* genome to ensure the inclusion of every CaALDH and the detection of the unpredicted ones. This resulted in 37 aldehyde dehydrogenases encoded by the chickpea genome.

**Tabla 1**. Chickpea ALDH proteins classification.
![](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/TABLA1.jpg)


*ALDH18* genes encode Δ1 -pyrroline-5-carboxylate synthetase (P5CS), defined as ALDH-like proteins (Sophos & Vasiliou, 2003). They are involved in proline biosynthesis (Igarashi et al., 1997), whose accumulation has adaptive roles in tolerance to biotic/abiotic stresses (Verbruggen & Hermans, 2008). Li et al. (2013) conclude that family 18 is the group that differs the most between species. The gene structure of the members of this family is different from that of the other families and presents additional domains AA-kinase; lacking the preserved ALDH active sites as shown by the results (Table 1). In addition, it is believed that the balance between proline biosynthesis and degradation is essential for its osmoprotective functions. Δ1 -pyrroline-5-carboxylate dehydrogenase (P5CDH) acts on the degradation, which is an ALDH12 protein. This degradation occurs in the mitochondria (Verbruggen & Hermans, 2008), which agrees with the results obtained in the prediction of subcellular localization (Table 1).

**Table 2**. CaALDH Characterization.
![](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/Tabla%202.jpg)

**Table 3**. Members of ALDH families identified in plants, humans and fungi.
![](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/All%20species%20ALDH%20families.jpg)



Most genes in family 3 seem to be regulated by the abscisic acid (ABA) stress response pathway, they get expressed in response to abiotic stresses (Brocker et al., 2013, Kirch et al., 2004). Family 3 in plants have diverged to six subfamilies: 3E, 3F, 3H, 3I, 3J and 3K. Two of these subfamilies are found in chickpea (3F and 3H). The 3H1 subfamily is very expanded with 10 members encoded by 7 genes, while 3F1 has 3 members encoded by 3 genes. Missihoun et al. (2012) have postulated the abundance of ALDH3s proteins as result of a complex expression pattern of their genes regulated by gene-splicing or alternative promoters. Based on the differences in expression and response in *A. thaliana* (Kirch et al., 2004), it was suggested that the isoforms of ALDH3 have evolved as a consequence of functional specialization in specific tissues and subcellular organelles; which fits with the diversity of subcellular localization predicted for the chickpea proteins of this family (Table 1). The algae described so far, *C. reinhardtii* unicellular and *V. carteri*  colonial, lack family 3; suggesting that the expansion and diversification of this gene family occurred along with the evolutionary movement of aquatic plants (Brocker et al., 2013).

ScanProsite analysis showed that the characteristic domains PS00687 and PF00070 were not found on the 37 ALDH sequences: 12 of the 37 contain both domains; 5 of them contain only domain PS00687 and 2 of them contain only PS00070. Some of the proteins that do not contain these domains are found in subfamilies 3F1, 3H1, and family 18. Therefore, it was necessary to make other alternative searches to identify these ALDH. In family 18, the domains PS00902 (Glutamate 5 kinase signature) and PS01223 (γ-glutamyl phosphate reductase signature) appear.

All ALDH gene families identified in higher plants such as *Arabidopsis* are present in *C. arietinum*, with the exception of ALDH19 present only in *S. lycopersicum* (Table 3). The chickpea, with 29 CaALDH genes, is the third one along with tomato with the highest number of ALDH genes. Above tham is apple tree with 39 genes and cotton with 30 genes. *C. arietinum* seems to have additional ALDH proteins for stress response. Families 3 and 18 are particularly abundant, which may be significant for carrying out the detoxification of aldehyde molecules generated under different stresses and maintaining the homeostasis of reducing equivalents.

*CaALDH* genes were mapped to the chickpea genome in order to understand the chromosomal distribution. Based on the available genome assembly of *C. arietinum*, 23 of the 29 genes are distributed in seven of the eight chromosomes. We find genes that encode ALDH from different families on the same chromosome. We could not map genes whose proteins have short sequence length (<272 aa). Chromosomes 6 and 7 contain 52% of the mapped genes. No *CaALDH* gene is mapped on chromosome 2 (Figure 3). The six *CaALDH* genes that can not be mapped in the reference genome (Table 2) could be mitochondrial genes; since the mitochondria is not found in the distribution of the reference genome in the NCBI.

![](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/Todos%20chr.jpg)
**Figure 3**. Genomic distribution of *CaALDH* genes on the chickpea chromosomes. Chromosome number and sizes (Mb) are indicated above and below each bar, respectively. Only chromosomes with *CaALDH* genes are represented.
