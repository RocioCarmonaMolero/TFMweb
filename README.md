# 'Identification and Characterization of the *Cicer arietinum* aldehyde dehydrogenase superfamily'



## 1. Identification of the *Cicer arietinum* aldehyde dehydrogenase superfamily

### 1.1. Search for candidate proteins

First, a keyword-based search was carried out in the 'Protein' database of the National Biotechnology Information Center ([NCBI](https://www.ncbi.nlm.nih.gov/)) specifying *Arabidopsis* as organism and 'Aldehyde dehydrogenase' as title. 136 proteins were obtained and their sequences were downloaded into a text file. A BLASTP (Altschul et al. 1990) was run against *Cicer arietinum* proteome in the Reference Protein database (refseq_protein). The thresholds established to select candidate *C. arietinum* ALDH were: Query cover ≥ 25%, E-value ≥ e-25, Identity ≥ 25%. The chickpea sequences that fulfilled these parameters were downloaded into a csv file.

As each BLASTP was done indepedently, the existence of the same XP was likely in different csv files. In order to clean those files and avoid the accession repetition enclosing the unique XPs on a vector called 'candidates', the R written function [get_IDprot()](https://github.com/RocioCarmonaMolero/ScriptProteinas/blob/master/get_IDprot.R) followed by the next commands was executed:

> dat <-  read.csv(1_2cicer_D9V2GJ7401R-Alignment-HitTable.csv , header = FALSE, stringsAsFactors = F) 
 
> dat <- dat[,2]     
#the second column of the csv, downloaded directly from the BLASTP result, is taken
 
> IDprot <- get_IDprot(dat, 2)     
#the function is applied and the result is saved in the 'IDprot' vector
 
> candidates <- c(candidates, IDprot[!IDprot %in% candidates])     
#the accessions of 'IDprot' not found in 'candidates' are included in it


This step was carried out with every 136 csv.

The same process was repeated twice, being *Medicago truncatula* and *Glycine max* (soybean) the query. There were 36 aldehyde dehydrogenases found in [Phytozome v12.1](https://phytozome.jgi.doe.gov/pz/portal.html) and 21 in Refseq_Protein for *Medicago* and 55 found for soybean in the NCBI.
The final 'candidates' vector contains 43 proteins.

Considering that the search was carried out against the model plant *Arabidopsis thaliana*, and against the legumes *Glycine max* and *Medicago truncatula*, which is the model of this family; we can conclude to have covered an adequate spectrum of species that, by sequence homology, will allow us to identify all the aldehyde dehydrogenases of *Cicer arietinum*.


### 1.2. Check ALDH conserved domains

To verify whether our candidate proteins are ALDH, the presence of the ALDH-superfamily conserved domains was manually checked: 'PF00171.21', 'PF07368.10' and 'PF05893.13' in Pfam; 'PS00687' and 'PS00070' in ScanProsite; the accession 'cl11961' in the database of Conserved Domains of the NCBI; and the accession '53720' in the Superfamily database.
After this screen 36 proteins are selected, which will constitute a new vector called 'aldhcicer'.

### 1.3. BLASTP 'aldhcicer' against chickpea genome 

Once the aldehyde dehydrogenases of chickpea are identified, a BLASTP is performed against its own genome. In this way, it will be possible to include any protein that does not have high homology with the other species and detect the unpredicted ones, since the query sequences are from the same organism.
One more protein was identified with the Pfam domain 'PF00171.21'. Thus, we conclude that the chickpea ALDH superfamily (CaALDHs) has 37 members.


## 2. Characterization of CaALDH
 
For the characterization of CaALDH, attention was focused on obtaining the following properties:
* Amino acids number (naa)
* mRNA accession (accRNA)
* RNA length in bp (RNAlen)
* Gene locus identifier (LOC)
* Chromosome location (chr)
* Exon number (exon)
* Start position of the protein (startp)
* Stop position (stop)
* Molecular weight (mol_wg)

In order to collect this information, the function [caracterization()](https://github.com/RocioCarmonaMolero/ScriptProteinas/blob/master/InformationProteins.R) is applied; naming the function [weight_prot()](https://github.com/RocioCarmonaMolero/ScriptProteinas/blob/master/get_mol_weight.R).

Open access bioinformatic tools were also used to predict the subcellular location and the active center. This knowledges will also be helpful for the classification of proteins.
For cell localization, [DeepLoc 1.0](http://www.cbs.dtu.dk/services/DeepLoc/), which differentiates between 10 eukaryotic protein locations: cytoplasm, extracellular, mitochondria, cell membrane, chloroplast, golgi apparatus, lysosome/vacuole, nucleus, endoplasmic reticulum and peroxisome; [SLP-Local](http://sunflower.kuicr.kyoto-u.ac.jp/~smatsuda/slplocal.html) and [TargetP 1.1](http://www.cbs.dtu.dk/services/TargetP/) were used. Moreover, DeepLoc indicates whether it is a soluble or membrane protein. In the membrane protein results, the transmembrane domain was verified by [Smart](http://smart.embl-heidelberg.de/)(Simple Modular Architecture Research Tool). In order to consolidate the chloroplast and mitochondria location results, the presence of chloroplast transit peptides (cTP) was confirmed by [ChloroP 1.1](http://www.cbs.dtu.dk/services/ChloroP/) and mitochondrial-targeting sequences (mTP) by [Mitoprot](https://ihg.gsf.de/ihg/mitoprot.html). 
For the active center, the protein sequences were analyzed with [PROSITE](https://prosite.expasy.org/), looking up the identification of the active center PS00687 (glutamic acid active site, E) and PS00070 (cysteine active site, C). In case any other active site was found, it was also noted down.

[PROPSEARCH](http://abcis.cbs.cnrs.fr/propsearch/) was used to determine the molecular function. Each protein was individually analyzed and the hits with the highest probability were written down. For the verification of these, it was checked whether the functional residues were conserved between query and hit. First through a database BLAST search and after through an alignment between them with the [Clustal Omega](http://www.clustal.org/omega/) program.

![Scheme of the bioinformatic tools used for the characterization step](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/function%20scheme.jpg)
**Figure 1**.Scheme of the bioinformatic tools used for the characterization and its aim


## 3. Classification of CaALDHs

For the classification, the criteria stablished by the ALDH Gene Nomenclature Committee (AGNC) in 1999 was applied. Two proteins belong to the same gene family if they have more than 40% identity; to the same subfamily if they have more than 60% identity. In the annotation, the root 'ALDH' is followed by a family descriptor number, a capital letter to describe the subfamily, a number specifying the individual gene within the subfamily and a lowercase letter if necessary to designate variants (Zhu et al., 2014).
Frequent methods for classifying this family in plants are based on homology with other plant species that are already described. For this, BLASTP of all 37 chickpea aldolases was made against the Protein Refseq database of *Medicago truncatula* and *Glycine max*, both being legumes and *Medicago* the model plant of them. All the results were downloaded and filtered to eliminate those whose identity was <40% and whose query length was less than the hit length: [Blast_Sieve.R](https://raw.githubusercontent.com/RocioCarmonaMolero/ScriptProteinas/master/Blast_Sieve.R)


## 4. Phylogeny

A multiple alignment of the sequences of the aldehyde dehydrogenases of *Medicago truncatula* (MtALDH) and *Cicer arietinum* (CaALDH) was performed in the form of one gene per locus with the MUSCLE program using the default parameters (Edgar, 2004). To deduce the evolutionary history, the Neighbor-Joining method was used (Saitou & Nei, 1987). The consensus bootstrap tree was inferred from 1000 repetitions (Felsenstein, 1985). Evolutionary distances were calculated using the Poisson correction method (Zuckerkandl & Pauling, 1965), being in the units of the number of amino acid substitutions per site. Evolutionary analyzes were performed in MEGA6 (Tamura et al., 2013). We consider sister pairs those proteins grouped on the basis of high bootstrap values (> 65%) [(Die et al., 2018)](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4695-9).


## 5. Duplication Analysis

For the duplication analysis of the CaALDH, the Circoletto program (Darzentas, 2010) was used. Both the FASTA query and the FASTA database are the protein sequences of CaALDH; with ultra-strict E-value values (10 -180), using absolute score / colored bands and colors: green for identities ≤ 95%, orange ≤ 99% and red> 99%.


## 6. Expression in silico

The coding sequences of the CaALDH genes were used as a query against the chickpea NCBI EST database. The search parameters were established as follows: megablast, identity> 90%, length> 180 bp and E-value <10-10 [(Die et al., 2018)](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4695-9).


## 7. Code availability

The codes of this TFM have been written in the R programming language, with the R software (Team R.C., 2017) and the free access RStudio interface (Team R, 2016, http://www.rstudio.com/). With these scripts we have collected data from the NCBI and performed the analysis of them. They are available in the repository: https://github.com/RocioCarmonaMolero/ScriptProteinas. The code is distributed under the MIT open source license.

Every step is deeply explained in the [TFM](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/TFM.md) file.


Find out the [RESULTS](https://github.com/RocioCarmonaMolero/TFMweb/blob/master/Results%20TFM.md) clicking here.
