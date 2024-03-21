import sys
import random
from bitarray import bitarray
from cyvcf2 import VCF
import gffutils
from succinct import rle_bit_array

# load gff files
gff = sys.argv[1]

# load vcf
variantFile = sys.argv[2]

dg = gffutils.create_db(gff, dbfn="test.db", force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

product = dict()
# read chromosome length
chrLength = dict()
for gene in dg.features_of_type('region'):
    # 1-based
    # print(gene.seqid, gene.start, gene.end, sep="\t" ) 
    chrLength[ str(gene.seqid) ] = int (gene.end) - int (gene.start ) + 1 

# create bit_array
bitArrayDict = dict() 
for i in chrLength:
    # size =chrLength[i] + 1
    bitArrayDict[i] = bitarray( chrLength[i] + 1 ) 


#for i in range(10):
#    j = random.randint(0,k-1)
#    bit_array[j] = True
#a = rle_bit_array.RunLengthEncodedBitArray(bit_array)

# target = CDS, exon, gene, ...
target = "CDS"
for gene in dg.features_of_type( target):
    # write boundaries of target region to bit_array
    chromosome, start, end = gene.seqid, int(gene.start), int(gene.end) 
    #print( chromosome, start, end )
    bitArrayDict[i][ start   ] = True
    bitArrayDict[i][ end + 1 ] = True
    #print(gene.seqid, gene.start, gene.end, gene.attributes["protein_id"],  sep="\t" ) 
    #print(gene.seqid, gene.start, gene.end, gene.featuretype, gene.attributes, sep="\t" ) 

# convert bit_array to succinct data structure
succinctArrayDict = dict() 
for i in bitArrayDict:
    succinctArrayDict[i] = rle_bit_array.RunLengthEncodedBitArray( bitArrayDict[i] )

#for i in succinctArrayDict:
#    print( i, len(succinctArrayDict[i] ) ) 

# create the CDS dictionary using the rank function
CDSdict = dict()
for gene in dg.features_of_type( target):
    # write boundaries of target region to bit_array
    chromosome, start, end = gene.seqid, int(gene.start), int(gene.end) 
    i = succinctArrayDict[ chromosome ].rank( start ) 
    CDSdict[i] =  gene.attributes["protein_id"][0] 

#for i in CDSdict:
#    print( i, CDSdict[i] )


inVCF = VCF( variantFile )
for line in inVCF:
    chrom, pos, ID = line.CHROM, line.POS, line.ID
    if ID==None or ID=="." or ID == "None":
        ID = chrom + ":" + str(pos)
    CDSid = succinctArrayDict[chrom].rank( pos )
    if CDSid in CDSdict: # if variant is on CDS
        print( ID, CDSdict[CDSid] )
