ParaHaplo5 has a variety of programs for conducting GWAS and SKAT.
GroupingSNPsByCDS is an example of it.

The purpose of GroupingSNPsByCDS.py is to classify SNPs into each gene. This program reads a GFF file defining gene regions and a VCF file containing genetic polymorphism information. The program only targets SNPs within CDS regions.

To run GroupingSNPsByCDS.py, installation of the bitarray package and the succinct package is required
pip install bitarray
pip install succinct

To run GroupingSNPsByCDS.py, specify the name of the GFF file defining gene regions as the first argument, and the name of the VCF file containing genetic polymorphism information as the second argument.

python GroupingSNPsByCDS.py genes.gff variants.vcf


