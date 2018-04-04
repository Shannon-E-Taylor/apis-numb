Important preprocessing of data: 
- SCRMshaw requires the scaffold names in the .gff file and .fa file to be the same- ie both 'Group1.1'. Remove any junk parts of scaffold names. 
- we used the Amel_4.5 and Amel_2.0 assemblies from Beebase. The .gff files differ in their assignment of genes and exons- Amel_4.5 calls them 'genes' and 'exons' while Amel_2.0 has "mRNA" and "CDS". The R script extract_genes_and_exons.R may need to be modified accordingly. (See also Amel_2.0 branch). 
- SCRMshaw fails quietly. Check all the output folders. 
- need to change line105 of source SCRMshaw code
