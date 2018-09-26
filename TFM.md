## Identification of the *Cicer arietinum* aldehyde dehydrogenase superfamily


### Search for candidate proteins

First, a keyword-based search was carried out in the 'Protein' database of the National Biotechnology Information Center ([NCBI](https://www.ncbi.nlm.nih.gov/)) specifying *Arabidopsis* as organism and 'Aldehyde dehydrogenase' as title. 136 proteins were obtained and their sequences were downloaded into a text file. A blastp (Altschul et al. 1990) was run against *Cicer arietinum* proteome in the Reference Protein database (refseq_protein). The thresholds established to select candidate aldolases from *C. arietinum* were: Query cover ≥ 25%, E-value ≥ e-25, Identity ≥ 25%. The chickpea sequences that fulfilled these parameters were downloaded into a csv file.

As each blastp was done indepedently, the existence of the same XP *Cicer* aldolase was likely in several different csv files. For this, an R code script was written to clean those files so that a vector called 'candidates' would be created. In this 'candidates' vector the accessions of the *Cicer* proteins of the 136 csv were included without being repeated. It was possible executing the function [get_IDprot()](https://github.com/RocioCarmonaMolero/ScriptProteinas/blob/master/get_IDprot.R).
The results of applying the function get_IDprot () are added to the vector 'candidates'; avoiding the accession repetition using the following commands:

> dat <-  read.csv( 1_2cicer_D9V2GJ7401R-Alignment-HitTable.csv , header = FALSE, stringsAsFactors = F)
 
> dat <- dat[,2]
 
> IDprot <- get_IDprot(dat, 2)
 
> candidatos <- c(candidatos, IDprot[!IDprot %in% candidatos])


First of all, the second column of the csv, downloaded directly from the blastp result, is taken; then, the function is applied and the resulting vector is saved in 'IDprot'. The accessions of 'IDprot' not found in 'candidates' are included in it. After executing it for the 136 csv, we have all the proteins of *C.arietinum* in which later we will check the aldolase domains.
The same process was repeated twice, being *Medicago truncatula* and *Glycine max* (soybean) the query in each of them. There were 36 aldehyde dehydrogenases found in [Phytozome v12.1](https://phytozome.jgi.doe.gov/pz/portal.html) and 21 in Refseq_Protein for *Medicago* and 55 those found for soybean in the NCBI.
The final 'candidate' vector contains 43 proteins.

Considering that the search was carried out against the model plant *Arabidopsis thaliana*, and against two legumes: *Glycine max* and *Medicago truncatula*, which is a model in turn of this family; we can conclude to have covered an adequate spectrum of species that, by sequence homology, will allow us to identify all the aldehyde dehydrogenases of *Cicer arietinum*.


### Checking ALDH conserved domains

To verify which of our candidate proteins are aldolases, the presence of the ALDH-superfamily conserved domains was checked: 'PF00171.21', 'PF07368.10' and 'PF05893.13' in Pfam; 'PS00687' and 'PS00070' in ScanProsite; the accession 'cl11961' in the database of Conserved Domains of the NCBI; and the accession '53720' in the Superfamily database.
After this screen 36 proteins are selected, which will constitute a new vector called 'aldhcicer'.

### Blastp of 'aldhcicer' against the genome of *C.arietinum*

Once the aldehyde dehydrogenases of chickpea are identified, a blastp is performed against its own genome. In this way, it will be possible to include any protein that does not have high homology with the other species, since in this case the query proteins are from the same organism.
One more protein was identified with the Pfam domain 'PF00171.21'. Thus, we conclude that the chickpea ALDH superfamily has 37 members.


### Characterization of aldehyde dehydrogenases of C. arietinum
 
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


## Classification of *C. arietinum* aldehyde dehydrogenases

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
