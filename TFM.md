## Identification of the *Cicer arietinum* aldehyde dehydrogenase

### Search for candidate proteins

First, a keyword-based search was carried out in the 'Protein' database of the National Biotechnology Information Center ([NCBI](https://www.ncbi.nlm.nih.gov/)) specifying *Arabidopsis* as organism and 'Aldehyde dehydrogenase' as title. 136 proteins were obtained and their sequences were downloaded into a text file. A blastp (Altschul et al. 1990) was run against *Cicer arietinum* proteome in the Reference Protein database (refseq_protein). The thresholds established to select candidate aldolases from *C. arietinum* were: Query cover ≥ 25%, E-value ≥ e-25, Identity ≥ 25%. The chickpea sequences that fulfilled these parameters were downloaded into a csv file.

As each blastp was done indepedently, the existence of the same XP *Cicer* aldolase was likely in several different csv files. For this, an R code script was written to clean those files so that a vector called 'candidates' would be created. In this 'candidates' vector the accessions of the *Cicer* proteins of the 136 csv were included without being repeated. It was possible executing the function [get_IDprot()](https://github.com/RocioCarmonaMolero/ScriptProteinas/blob/master/get_IDprot.R).

