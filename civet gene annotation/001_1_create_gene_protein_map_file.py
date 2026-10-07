#!/usr/local/package/python/3.12.0/bin/python3

# Create protein-gene mapping file

from argparse import ArgumentParser
from sys import exit
import os
import subprocess
import math
import multiprocessing as mp
import queue
import time

# Create protein-gene mapping file
def queue_tasks(tasks_to_accomplish,tasks_that_are_done):
    while True:
        try:
            task=tasks_to_accomplish.get_nowait()
        except queue.Empty:
            break
        else:
            q_target_gpff_dir=task[0]
            q_target_map_output=task[1]
            q_target_prot_faa_dir=task[2]
            q_species_name=task[3]
                        
            if os.path.exists(q_target_map_output)==False:
                            
                ##### Uncompress gzip files if the files are in zip format
                if os.path.exists(q_target_gpff_dir+'.gz')==True:
                    cmd_gunzip='gunzip '+q_target_gpff_dir+'.gz'
                    print(cmd_gunzip)
                    subprocess.run(cmd_gunzip, shell=True)
                if os.path.exists(q_target_prot_faa_dir+'.gz')==True:
                    cmd_gunzip='gunzip '+q_target_prot_faa_dir+'.gz'
                    print(cmd_gunzip)
                    subprocess.run(cmd_gunzip, shell=True)
                                    
                ### Gather protein sequences of all proteins
                print('\nCreate protein-gene map for the species',q_species_name)
                prot_fasta_lines=[line.strip('\n') for line in open(q_target_prot_faa_dir)]
                fasta_dict={}
                for a in range(0,len(prot_fasta_lines)):
                    if prot_fasta_lines[a][0]=='>':
                        prot_fasta_id=prot_fasta_lines[a].strip('>').split(' ')[0]
                        prot_loop_index,prot_add_index=0,1
                        prot_seq=''
                        while prot_loop_index==0:
                            if (a+prot_add_index)!=len(prot_fasta_lines) and prot_fasta_lines[a+prot_add_index][0]!='>':
                                prot_seq+=prot_fasta_lines[a+prot_add_index]
                                prot_add_index+=1
                            else:
                                fasta_dict[prot_fasta_id]=prot_seq
                                prot_loop_index=1
                
                ### Gather information of proteins from the GPFF file
                target_gpff_lines=[line.strip('\n') for line in open(q_target_gpff_dir)]
                mapping_dict={}
                gene_dict={}
                unknown_id=1
                for k in range(0,len(target_gpff_lines)):
                    target_gpff_line_cmp=target_gpff_lines[k].split(' ')
                    if target_gpff_line_cmp[0]=='LOCUS':
                        locus_id=target_gpff_line_cmp[7]
                        prot_info_list=[locus_id]
                        other_accession_list=[]
                        loop_index,add_index=0,1
                        gene_name='UNK:'+'g'+str(unknown_id)+':'+q_species_name
                        coded_by='unknown'
                        db_xref='unknown'
                        gene_id='unknown'
                        synonyms='none'
                        while loop_index==0:
                            prot_line_cmp=target_gpff_lines[k+add_index].split(' ')
                            if prot_line_cmp[0]=='DEFINITION':
                                prot_def=target_gpff_lines[k+add_index][12:]
                                def_loop_index,def_add_index=0,1
                                while def_loop_index==0:
                                    if target_gpff_lines[k+add_index+def_add_index][0:9]=='ACCESSION':
                                        def_loop_index=1
                                    else:
                                        prot_def+=(' '+target_gpff_lines[k+add_index+def_add_index][12:])
                                        def_add_index+=1
                            if prot_line_cmp[0]=='ACCESSION':
                                if len(prot_line_cmp)>5:
                                    other_accession_list=prot_line_cmp[4:]
                                accession_number=prot_line_cmp[3]
                                prot_info_list.append(accession_number)
                            if prot_line_cmp[0]=='VERSION':
                                prot_id_ver=prot_line_cmp[5]
                                prot_info_list.append(prot_id_ver)
                            if '     CDS     ' in target_gpff_lines[k+add_index]:
                                feature_loop_index,feature_add_index=0,1
                                while feature_loop_index==0:
                                    feature_line=target_gpff_lines[k+add_index+feature_add_index]
                                    if '/gene=' in feature_line:
                                        gene_name=feature_line.split('/gene=')[1].strip('"').upper()
                                    if '/coded_by=' in feature_line:
                                        if "join(" not in feature_line:
                                            coded_by=feature_line.split('/coded_by=')[1].strip('"')
                                        elif feature_line.count('"')==2:
                                            coded_by=feature_line.split('/coded_by=')[1].strip('"')
                                        else:
                                            coded_by=feature_line.split('/coded_by=')[1].strip('"')
                                            coded_by_loop_index,coded_by_add_index=0,1
                                            while coded_by_loop_index==0:
                                                if '"' not in target_gpff_lines[k+add_index+feature_add_index+coded_by_add_index]:
                                                    coded_by+=target_gpff_lines[k+add_index+feature_add_index+coded_by_add_index]
                                                    coded_by_add_index+=1
                                                else:
                                                    coded_by+=target_gpff_lines[k+add_index+feature_add_index+coded_by_add_index].strip('"')
                                                    coded_by=coded_by.replace(' ','')
                                                    coded_by_loop_index=1
                                    ##### /gene_synonym can split into more than one line
                                    if '/gene_synonym=' in feature_line:
                                        if feature_line[-1]!='"':
                                            to_add_synonym=feature_line
                                            synonym_loop_index,synonym_add_index=0,1
                                            while synonym_loop_index==0:
                                                if to_add_synonym[-1]==';':
                                                    to_add_synonym+=' '
                                                to_add_synonym+=' '.join(target_gpff_lines[k+add_index+feature_add_index+synonym_add_index].split(' ')[21:])
                                                if target_gpff_lines[k+add_index+feature_add_index+synonym_add_index][-1]!='"':
                                                    synonym_add_index+=1
                                                else:
                                                    synonyms=('|'.join(to_add_synonym.split('/gene_synonym=')[1].strip('"').split('; '))).upper()
                                                    synonym_loop_index=1
                                        else:
                                            synonyms=('|'.join(feature_line.split('/gene_synonym=')[1].strip('"').split('; '))).upper()
                                    ##### There can be more than one /db_xref line.
                                    if '/db_xref=' in feature_line:
                                        if db_xref=='unknown':
                                            db_xref=feature_line.split('/db_xref=')[1].strip('"')
                                        else:
                                            db_xref+=(','+feature_line.split('/db_xref=')[1].strip('"'))
                                        if 'GeneID' in feature_line:
                                            gene_id=feature_line.split('/db_xref=')[1].strip('"').split("GeneID:")[1]
                                    feature_add_index+=1
                                    if target_gpff_lines[k+add_index+feature_add_index][0:6]=='ORIGIN':
                                        feature_loop_index=1
                                loop_index=1
                            add_index+=1
                            if target_gpff_lines[k+add_index]=='//':
                                loop_index=1
                        prot_info_list=prot_info_list+[gene_name,str(len(fasta_dict[prot_id_ver])),coded_by,prot_def]
                        if len(other_accession_list)>0:
                            prot_info_list=prot_info_list+['|'.join(other_accession_list)]
                        else:
                            prot_info_list=prot_info_list+['none']
                        prot_info_list=prot_info_list+[synonyms,db_xref,gene_id]
                        if 'UNK:' in gene_name:
                            unknown_id+=1
                        ### Ignore redundancy
                        mapping_dict[locus_id]=prot_info_list
                        
    #                    ### Check for redundancy of the protein sequences whether they are redundant proteins that were translated from the same genes or not
    #                    ### In the case that gene name is unknown, redundancy of all sequences can still be checked! :)
    #                    if gene_name not in gene_dict:
    #                        gene_dict[gene_name]=set([fasta_dict[prot_id_ver]])
    #                        mapping_dict[locus_id]=prot_info_list
    #                    else:
    #                        if fasta_dict[prot_id_ver] not in gene_dict[gene_name]:
    #                            gene_dict[gene_name].add(fasta_dict[prot_id_ver])
    #                            mapping_dict[locus_id]=prot_info_list
                
                mapping_lines=['locus_name\taccession_no\tversion\tgene_name\tprot_length\tcoded_by\tdefinition\tother_versions\tsynonyms\tdb_xref\tGeneID\n']
                mapping_dict_keys=sorted(list(mapping_dict.keys()))
                for l in range(0,len(mapping_dict_keys)):
                    mapping_line='\t'.join(mapping_dict[mapping_dict_keys[l]])+'\n'
                    mapping_lines.append(mapping_line)
                output_file=open(q_target_map_output,'w')
                output_file.writelines(mapping_lines)
                output_file.close()

                ##### Compress the unzipped files
                if os.path.exists(q_target_gpff_dir)==True:
                    cmd_gzip='gzip '+q_target_gpff_dir
                    print(cmd_gzip)
                    subprocess.run(cmd_gzip, shell=True)
                if os.path.exists(q_target_prot_faa_dir)==True:
                    cmd_gzip='gzip '+q_target_prot_faa_dir
                    print(cmd_gzip)
                    subprocess.run(cmd_gzip, shell=True)

            tasks_that_are_done.put('A protein-gene map of the species '+q_species_name+' is already created by '+mp.current_process().name)
            time.sleep(.5)

if __name__ == '__main__':
    parser = ArgumentParser(description='Update gene name for protein ID in the OrthoFinder-produced matrix')
    parser.add_argument('genome_dir',help='Input a directory containing all genome GPFF files annotated from NCBI or Ensembl')
    parser.add_argument('species_genome_indices',help='Input the 009_2_species_genome_indices.txt file')
    args = parser.parse_args()
    if not args.genome_dir:
        parser.print_help()
        exit(1)
    if not args.species_genome_indices:
        parser.print_help()
        exit(1)

    number_of_processes=30
    tasks_to_accomplish=mp.Queue()
    tasks_that_are_done=mp.Queue()
    processes=[]

    # Gather all directories of GPFF file from species-genome index file
    species_genome_lines=[line.strip('\n') for line in open(args.species_genome_indices)]
    for i in range(1,len(species_genome_lines)):
        species_genome_line_cmp=species_genome_lines[i].split('\t')
        species_name='_'.join(species_genome_line_cmp[4].split(' '))
        species_genome=species_genome_line_cmp[0]+'_'+'_'.join(species_genome_line_cmp[1].split(' '))
        target_gpff_dir=args.genome_dir+'/'+species_genome+'/'+species_genome+'_protein.gpff'
        target_map_output=args.genome_dir+'/'+species_genome+'/'+species_genome+'_prot_gene_map.tsv'
        target_prot_faa_dir=args.genome_dir+'/'+species_genome+'/'+species_genome+'_protein.faa'
        tasks_to_accomplish.put((target_gpff_dir,target_map_output,target_prot_faa_dir,species_name))

    # Creating processes
    for w in range(number_of_processes):
        p=mp.Process(target=queue_tasks, args=(tasks_to_accomplish,tasks_that_are_done))
        processes.append(p)
        p.start()
        time.sleep(1)

    # Ensure all of the processes have finished
    for p in processes:
        p.join()

    # Print the output
    while not tasks_that_are_done.empty():
        print(tasks_that_are_done.get())
    print('\nORFs of all RTs were already extracted from the proteome set')
