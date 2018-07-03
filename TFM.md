## Identification of the *Cicer arietinum* aldehyde dehydrogenase

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


### Checking the conserved domains ALDH

To verify which of our candidate proteins are aldolases, the presence of the ALDH-superfamily conserved domains was checked: 'PF00171.21', 'PF07368.10' and 'PF05893.13' in Pfam; 'PS00687' and 'PS00070' in ScanProsite; the accession 'cl11961' in the database of Conserved Domains of the NCBI; and the accession '53720' in the Superfamily database.
After this screen 36 proteins are selected, which will constitute a new vector called 'aldhcicer'.

### Blastp of 'aldhcicer' against the genome of *C.arietinum*

Once the aldehyde dehydrogenases of chickpea are identified, a blastp is performed against its own genome. In this way, it will be possible to include any protein that does not have high homology with the other species, since in this case the query proteins are from the same organism.
One more protein was identified with the Pfam domain 'PF00171.21'. Thus, we conclude that the chickpea ALDH superfamily has 37 members.
