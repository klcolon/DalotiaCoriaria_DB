#!/usr/bin/python
import pandas as pd
import click
import numpy as np
from Bio.Seq import Seq

class pw2_fasta:
    def __init__(self, fasta, gff, bed=None):
        self.fasta = fasta
        self.bed = bed
        self.gff = gff

    def generate_fasta_df(self):
        """
        Read fasta file and convert to df.
        """
        _list = []
        with open(self.fasta, "r") as f:
            temp_gene = []
            temp_seq = []
            line = "init"
            while line != "":
                line = f.readline()
                if ">" in line and temp_gene == []:
                    temp_gene.append(line[1:].replace("\n",""))
                elif ">" in line and temp_gene != []:
                    temp_seq = ["".join(temp_seq)]
                    _list.append(temp_gene + temp_seq)
                    temp_gene = []
                    temp_seq = []
                    temp_gene.append(line[1:].replace("\n",""))
                elif line == "":
                    temp_seq = ["".join(temp_seq)]
                    _list.append(temp_gene + temp_seq)
                else:
                    temp_seq.append(line.replace("\n",""))
            f.close()
            
        df = pd.DataFrame(_list)
        df.columns = ["Gene","Sequence"]
        
        return df
    
    def generate_bed_df(self):
        """
        Create df from bed file
        """
        if type(self.bed) == type(None):
            print('No BED file provided..')
            return None
        df = pd.read_csv(self.bed, sep = "\t", header=None)
        df.columns = ["Gene", "Start", "End"]
        
        return df
    
    def generate_gff_df(self):
        """
        Create df from gff file
        """
        df = pd.read_csv(self.gff, sep = "\t", header=None)
        df.columns = ["Gene", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame","Attribute"]
        
        return df
    
    def create_fasta(filename, sequences, id_, product, biotype='mRNA'):
        """
        function to create fasta file
        """
        #fill fasta file
        with open(filename, "w+") as f:
            for i,seq in enumerate(sequences):
                f.write(f">{id_[i]} product:{product[i]}" \
                f" transcript_biotype:{biotype} \n")
                f.write(seq + "\n")
        f.close()
        # reference text file of ids and names
        with open('reference.csv','w+') as f:
            for i,x in enumerate(id_):
                f.write(x +',' + product[i]+'\n')
        f.close()

        
    def generate_probe_fasta(self, filename, biotype):
        
        """
        Use df generated from fasta to cut out sequences based on bed and gff file.
        """
        #create df for fasta, bed and gff file
        df_fasta = pw2_fasta.generate_fasta_df(self)
        df_bed = pw2_fasta.generate_bed_df(self)
        df_gff = pw2_fasta.generate_gff_df(self)
        
        #add extra info to df bed file regarding strand orientation and actual gene name
        strand_info = []
        name = []
        if biotype == "introns":
            print("Generating intron fasta sequences")
            assert type(self.bed) != type(None), 'Missing BED file'
        elif biotype == "exons":
            print("Generating exon fasta sequences")
            assert type(self.bed) != type(None), 'Missing BED file'
        else:
            df_gff = df_gff[df_gff.Feature == 'mRNA'].reset_index(drop=True)
            print("Generating mRNA fasta sequences using GFF only")
        
        if biotype != 'mRNA':
            for gene, end in df_bed[["Gene", "End"]].values:
                if biotype == "introns":
                    #start of gene exon is the end of the gene intron
                    strand = df_gff[(df_gff.Gene == gene) & (df_gff.Start == end)].Strand.iloc[0]
                    attribute = df_gff[(df_gff.Gene == gene) & (df_gff.Start == end)].Attribute.iloc[0]
                else:
                    #any the exon end should match with gff end
                    strand = df_gff[(df_gff.Gene == gene) & (df_gff.End == end)].Strand.iloc[0]
                    attribute = df_gff[(df_gff.Gene == gene) & (df_gff.End == end)].Attribute.iloc[0]  
                #id
                id_ = attribute.split(";")[0].replace("ID=", "")   
                #product
                if "product" in attribute:
                    product = attribute.split(";")[-1].replace("product=", "")
                      
                if 'product' in locals():
                    name.append([id_, product])
                else:
                    name.append([id_, 'NA'])
                strand_info.append(strand)
        else:
            for i in range(len(df_gff)):
                attribute = df_gff.Attribute.iloc[i]
                #id
                id_ = attribute.split(";")[0].replace("ID=", "")   
                #product
                if "product" in attribute:
                    product = attribute.split(";")[-1].replace("product=", "")    
                if 'product' in locals():
                    name.append([id_, product])
                else:
                    name.append([id_, 'NA'])
            
        if biotype != 'mRNA':     
            df_bed.loc["Strand"] = strand_info
            df_bed.loc["id"] = np.array(name)[:, 0]
            df_bed.loc["product"] = np.array(name)[:, 1]
        else:
            df_bed = df_gff[['Gene','Start', 'End','Strand']].copy()
            df_bed["id"] = np.array(name)[:, 0]
            df_bed["product"] = np.array(name)[:, 1]
        
        #isolate sequences based on bed file
        gene_sequence_final = []
        for i in range(len(df_bed)):
            #grab all info regarding gene product
            gene_name = df_bed.iloc[i].Gene
            gene_start = df_bed.iloc[i].Start
            gene_end = df_bed.iloc[i].End
            gene_orientation = df_bed.iloc[i].Strand   
            gene_seq = df_fasta[df_fasta["Gene"] == gene_name].Sequence.iloc[0]
            #cut fragment while accounting for python indexing
            gene_fragment = gene_seq[gene_start-1:gene_end-1]
            #if the anti sense strand is read, then take reverse complement 
            if gene_orientation == "-":
                gene_fragment = str(Seq(gene_fragment).reverse_complement())
            gene_sequence_final.append(gene_fragment)
            
        pw2_fasta.create_fasta(filename, gene_sequence_final, 
                               df_bed["id"].values.tolist(), df_bed["product"].values.tolist(), biotype)

@click.command()
@click.option('--fasta','-f', "fasta", help='FASTA File of organism genome',required=False, type=str)
@click.option('--gff','-gf', "gff", help='gff file of organism', type=str, required=False)
@click.option('--biotype','-bio', "biotype",help='what type of molecule is this? (ex.introns, exon, mRNA)', required=True, type=str)
@click.option('--sequences','-s', "sequences",help='Probe sequence csv file', required=False, type=str, default="")
@click.option('--genenames','-g', "genenames",help='Gene names csv file for probe sequences', required=False, type=str)
@click.option('--filename','-fi', "filename",help='output file name', required=True, type=str)
@click.option('--bed','-b', "bed",help='BED file containing genomic coordinates', required=False, type=str)

#overall wrapper function
def main(fasta, gff, biotype, sequences, genenames, filename, bed=None):
    if sequences == "":
        print("Creating FASTA file...")
        pw2f = pw2_fasta(fasta = fasta, gff=gff, bed=bed)
        pw2f.generate_probe_fasta(filename, biotype)
    else:
        print("Creating FASTA file...")
        sequences = pd.read_csv(sequences, header=None)[0].values.tolist()
        genenames = pd.read_csv(genenames, header=None)[0].values.tolist()
        pw2_fasta.create_fasta(filename, sequences, genenames, genenames, biotype)
    
if __name__ == '__main__':
    main()
