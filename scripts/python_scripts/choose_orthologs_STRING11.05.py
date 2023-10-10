import os, sys
import argparse
from collections import defaultdict

def preload_data(groups_folder,species_file,tree_file):

    
    # # this three files can be downloaded directly from string and eggnog files directly 

    
    
    print("Loading taxids for all STRING species...", file=sys.stderr)
    
    all_sps = set()
    
    
    for line in open(species_file):
        if line[0] == "#":
            continue
        l = line.strip().split("\t")
        sp = l[0]
        domain = l[4]
        if domain == "Bacteria":
            all_sps.add(line.split()[0])


    ##  
    ## Load groups
    ##  
    
    print("Loading groups...", file=sys.stderr)
    
    level_prot_og = {} # contain information for each taxonomic level in eggNOG, mapping from prot to og
                       # prot is key 
    level_og_prot = {} # contain information for each taxonomic level in eggNOG, mapping from prot to og
                       # og is key 
    
    for group_file in os.listdir(groups_folder):
    
        if "info" in group_file:
            continue
        if os.path.isdir(groups_folder+"/"+group_file):
            continue
    
        level = group_file.split(".")[0]
    
        level_prot_og[level] = {}
        level_og_prot[level] = {}
        
        for line in open(groups_folder+"/"+group_file):
            if not line.startswith("#"):
                l = line.strip().split("\t")
    
                og, sp, prot = l[1], l[2], l[3]
                
                if sp not in all_sps:
                    continue
    
                if prot not in level_prot_og[level]:
                    level_prot_og[level][prot] = set()
    
                level_prot_og[level][prot].add(og)
        
                if og not in level_og_prot[level]:
                    level_og_prot[level][og] = set()
                level_og_prot[level][og].add(prot)
    
    
    
    
    ##  
    ## Load full tree
    ##  
    
    print("Loading NCBI tree...", file=sys.stderr)
    
    tree_up_full = {}
    for line in open(tree_file):
        if line.startswith("#"): continue

        l = line.split("\t")
        child, parent = l[0], l[1]
        tree_up_full[child] = parent
        
        
    ##  
    ## Prune levels to sps and eggNOG levels and species 
    ##  
    
    print("Pruning NCBI tree...", file=sys.stderr)
    
    levels = set(level_prot_og.keys())
    
    tree_up = {}
    for sp in all_sps:
    
        taxid = sp
        prev_taxid = taxid
    
        while taxid != "1":
    
            taxid = tree_up_full[taxid]
            if taxid in levels:
                tree_up[prev_taxid] = taxid
                prev_taxid = taxid

        tree_up[prev_taxid] = "1"


    return level_prot_og, level_og_prot, all_sps, levels, tree_up





def choose_intital_protein(proteins):
    return proteins[0]

def is_child(child, parent):
    if child == parent:
        return True
    else:
        while child != "1":
            child = tree_up[child]
            
            if child == parent:
                return True
        return False
    

def resolve(input_protein, max_level, orthologs):
    

    input_species = input_protein.split(".")[0]
    orthologs[input_species] = [input_protein, "init-"+max_level]

    #print('START', 'input_protein:', input_protein, 'max_level', max_level, 'is:', is_child(input_species, max_level))
    
    current_level = input_species
    
    while current_level != max_level:
        
        
        current_level = tree_up[current_level]
        
#         if not is_child(input_species, current_level):
#             break
            
        if input_protein in level_prot_og[current_level]: 
            
            ## check for multiple ogs for protein per level

            current_species_proteins = {}
            
            for og in level_prot_og[current_level][input_protein]:    #here considering possibility 
                
                # that input ecoli protein may esist in multiple OGs???
                # two situation, in one level, one ecoli in two OG, is it possible?  to make it worse, two OG has proteins from same speices ??/ yes it is possible ; this problem is not sovle in this script ; here protins from both OG are choose, 
                # not correct, as they are from different domains. ,problem could be solved by find overlaping genes in 
                # OG in our wanted taxonimics level and then choose OG with more overlaped proteins 
                # or in one OG, one ecoli in one OG has two paralog from same species 
                for ortholog in level_og_prot[current_level][og]:
                    
                    sp_ortholog = ortholog.split(".")[0] # sp_ortholog is a speice name 
                    
                    
                    if sp_ortholog not in orthologs:
                        if sp_ortholog not in current_species_proteins:
                            current_species_proteins[sp_ortholog] = []
                        current_species_proteins[sp_ortholog].append(ortholog)
                
            
            not_resolved_paralogs = []
            for sp_ortholog in current_species_proteins:
                                
                if len(current_species_proteins[sp_ortholog]) > 1:
                    
                    not_resolved_paralogs.append(current_species_proteins[sp_ortholog][0])
                    
                else:
                    if sp_ortholog not in orthologs: # herre ensure one species only appear once ?? 
                        
                        orthologs[sp_ortholog] = [current_species_proteins[sp_ortholog][0], current_level]
                        
            if not_resolved_paralogs:
                # here i have a question, how to make sure every proein in not_resolved_paralogs will be included in the end ???
                initial_protein = choose_intital_protein(not_resolved_paralogs)
                
                if not is_child(initial_protein.split(".")[0], current_level):
                    print('ERROR:', initial_protein, current_level, max_level)
                    print("not_resolved_paralogs len : ", len(not_resolved_paralogs))
                
                else:
                    resolve(initial_protein, current_level, orthologs) #???? what does this mean, where is assigned to ??
                # only used to change orthologs ? what is name of this python mechanism 
                # whye is this command necessary, who not just add intial protien ahere ot orthologs 
            

if __name__ == '__main__':
# run code in this way : python choose_orthologs_STRING11.05.py -o "fsd" -i "asdf" -m "asdfa"fsd asdf asdfa
    
    parser = argparse.ArgumentParser(description='getOGsFromEggNOGForAllProteins')
    parser.add_argument('-o','--currentSpe_currentMaxLevel_orthologs', type=str, help='output folder')
    parser.add_argument('-i','--currentSpe_fastaData', type=str, help='fasta file of current speceis')
    parser.add_argument('-m','--max_orthology_level', type=str, help='max eggnog level belong which we find ogs of protein in current species')
    parser.add_argument('-g','--groups_folder', type=str, help='eggnog orthoglog group information')
    parser.add_argument('-s','--species_file', type=str, help='STRING species file')
    parser.add_argument('-t','--tree_file', type=str, help='STRING species tree file')
    args = parser.parse_args()
    currentSpe_currentMaxLevel_orthologs=args.currentSpe_currentMaxLevel_orthologs
    currentSpe_fastaData=args.currentSpe_fastaData
    max_orthology_level=args.max_orthology_level
    groups_folder=args.groups_folder
    species_file=args.species_file
    tree_file=args.tree_file

    

    print(currentSpe_currentMaxLevel_orthologs,currentSpe_fastaData,max_orthology_level)

    ##
    ## Load species, tree and groups
    ##

    #edit here with the input 
    level_prot_og, level_og_prot, all_sps, levels, tree_up = preload_data(groups_folder,species_file,tree_file)


    ##
    ## Load inital proteins
    ##

    input_proteins = []
    for line in open(currentSpe_fastaData):
        if line.startswith(">"):
            prot = line.strip()[1:].split()[0]
            input_proteins.append(prot)


    ##
    ## Compute orthologs
    ##


    for input_protein in input_proteins:

        fh = open(currentSpe_currentMaxLevel_orthologs+input_protein+".orthologs", "w")
        orthologs = {}


        input_protein_species = input_protein.split(".")[0]

        if not is_child(input_protein_species, max_orthology_level):
            #here happend when NCBI level is not an eggNOG level 
            sys.exit(f"Species {input_protein_species} is not child of chosen maximum orthology level... {max_orthology_level}")

        print("input_protein is : " , input_protein)
        resolve(input_protein, max_orthology_level, orthologs)

        for sp_ortholog in orthologs:
            print(sp_ortholog, orthologs[sp_ortholog][0], orthologs[sp_ortholog][1], sep="\t", file=fh)
        fh.close()

                    





