#!/usr/local/package/python/3.6.5/bin/python3

# Update gene name for protein ID in the OrthoFinder-produced matrix

from argparse import ArgumentParser
from sys import exit
import os
import subprocess

def update_gene_name_to_orthofinder_table(orthofinder_table,genome_dir,species_genome_indices,output_dir):

    # Create species-genome dict
    species_genome_dict={}
    species_genome_lines=[line.strip('\n') for line in open(args.species_genome_indices)]
    for i in range(2,len(species_genome_lines)):
        species_genome_line_cmp=species_genome_lines[i].split('\t')
        species_name='_'.join(species_genome_line_cmp[7].split(' '))
        species_genome=species_genome_line_cmp[0]+'_'+'_'.join(species_genome_line_cmp[15].split(' '))
        species_genome_dict[species_name]=species_genome

    genomic_lines=[line.strip('\n') for line in open(orthofinder_table)]
    species_num=(len(genomic_lines[0].split('\t'))-1)
    species_name_list=genomic_lines[0].split('\t')[1:]
    new_genomic_lines=[]
    new_genomic_lines.append(genomic_lines[0]+'\n')

    # Create gene dictionary for each species
    for j in range(0,species_num):
        ### Gather information of all proteins and their gene origins from mapping files
        if species_name_list[j]=='Paguma_larvata':
            continue
        prot_gene_map_file=genome_dir+'/'+species_genome_dict[species_name_list[j]]+'/'+species_genome_dict[species_name_list[j]]+'_prot_gene_map.tsv'
        prot_gene_map_lines=[line.strip('\n') for line in open(prot_gene_map_file)]
        prot_gene_map_dict={}
        for k in range(1,len(prot_gene_map_lines)):
            prot_gene_map_line_cmp=prot_gene_map_lines[k].split('\t')
            if prot_gene_map_line_cmp[1] not in prot_gene_map_dict:
                prot_gene_map_dict[prot_gene_map_line_cmp[1]]=prot_gene_map_line_cmp[3]
            else:
                print('There are more than 1 gene name for the protein accession no.',prot_gene_map_line_cmp[1])
                exit(1)

        # Update gene name into the charater-state matrix
        for n in range(1,len(genomic_lines)):
            genomic_line_component=genomic_lines[n].split('\t')
            gene_list=genomic_line_component[j+1].split(', ')
            for o in range(0,len(gene_list)):
                if gene_list[o].split('.')[0] in prot_gene_map_dict:
                    gene_list[o]=prot_gene_map_dict[gene_list[o].split('.')[0]].upper()
            gene_list=list(set(gene_list))
            genomic_line_component[j+1]=','.join(gene_list)
            genomic_lines[n]='\t'.join(genomic_line_component)
            
    # Write result file - gene matrix
    new_genomic_lines=[]
    for p in range(0,len(genomic_lines)):
        new_line_component=[]
        for r in range(0,species_num):
            new_line_component.append(genomic_lines[p].split('\t')[r+1])
        new_line=genomic_lines[p].split('\t')[0]+'\t'+'\t'.join(new_line_component)+'\n'
        new_genomic_lines.append(new_line)

    output_file=open(output_dir,'w')
    output_file.writelines(new_genomic_lines)
    output_file.close()

if __name__ == '__main__':
    parser = ArgumentParser(description='Update gene name for protein ID in the OrthoFinder-produced matrix')
    parser.add_argument('orthofinder_table',help='Input gene_matrix_with_only_changed_states.tsv')
    parser.add_argument('genome_dir',help='Input a directory containing all genome GPFF files annotated from NCBI or Ensembl')
    parser.add_argument('species_genome_indices',help='Input the 009_2_species_genome_indices.txt file')
    parser.add_argument('output_dir',help='Add output name and directory')
    args = parser.parse_args()
    if not args.orthofinder_table:
        parser.print_help()
        exit(1)
    if not args.genome_dir:
        parser.print_help()
        exit(1)
    if not args.species_genome_indices:
        parser.print_help()
        exit(1)
    if not args.output_dir:
        parser.print_help()
        exit(1)
    update_gene_name_to_orthofinder_table(args.orthofinder_table,args.genome_dir,args.species_genome_indices,args.output_dir)
