import numpy as np 
from Bio.PDB import PDBParser
import xpdb
# import py3Dmol


# also check  in  local computer http://localhost:8888/lab/tree/Documents/PhD_Tao/Classes/PyMOL/plot_3d_PDBStructure.py

def read_allres_pdbformat_data(file_name):
    """
    read  all reasidues from pdb format file 
    """

    # Create sloppy parser
    sloppy_parser = PDBParser(structure_builder=xpdb.SloppyStructureBuilder())
    pro_stuc = sloppy_parser.get_structure(id=None, file=file_name)

    # Structure 1
    sloppyio= xpdb.SloppyPDBIO()
    sloppyio.set_structure(pro_stuc)
    #sloppyio_1.save("sloppyio_1.pdb")

    # Get protein residue structures
    # https://alphafold.ebi.ac.uk/faq
    # Note that the PDB and mmCIF files contain coordinates for all regions, regardless of their pLDDT score. 
    # It is up to the user to interpret the model judiciously, in accordance with the guidance above.
    pro_allRes = [x for x in sloppyio.structure.get_residues()]
    
    return pro_allRes


# Calculate residue distance
def calculate_residue_dist(residue_one, residue_two):
    diff_vector = residue_one["CA"].coord - residue_two["CA"].coord
    sq_dist = np.sqrt(np.sum(diff_vector * diff_vector))
    return sq_dist

def calculate_dist_matrix(chain_one, chain_two):
    answer = np.zeros((len(chain_one), len(chain_two)), np.float)
    for row, residue_one in enumerate(chain_one):
        for col, residue_two in enumerate(chain_two):
            euclidean_dist = calculate_residue_dist(residue_one, residue_two)
            answer[row, col] = euclidean_dist 
    return answer

def generate_proximity_matrix( seq_1, seq_2, angstroms=10):

    # Select the residues from maps that are less than 'n' ansgtoms apart
    contact_map = calculate_dist_matrix(seq_1, seq_2)
    adjacency_matrix = np.zeros(np.shape(contact_map))
    adjacency_matrix[contact_map < angstroms] = 1


    return contact_map,adjacency_matrix



def get_pdb_stuc_chainRes_dict(pdb_file,selected_chainIDs=[]):

    sloppy_parser = PDBParser(structure_builder=xpdb.SloppyStructureBuilder())
    pro_stuc = sloppy_parser.get_structure(id=None, file=pdb_file)

    pro_stuc_chainRes=dict()
    for chain in pro_stuc.get_chains():
        print(chain.get_id())
        if chain.get_id() in selected_chainIDs:
            pro_stuc_chainRes[chain.get_id()]=[x for x in chain.get_residues()]
            print(chain.get_id(),len(pro_stuc_chainRes[chain.get_id()]))

    
    return(pro_stuc_chainRes)


def map_topRankingBetValue_list_to_realPDB_CACoords(topRankingBetValue_list,pro_stuc_chainRes,chainList):
    start_coords_list=list()
    end_coords_list=list()
    chain1,chain2=chainList
    for i in range(0,len(topRankingBetValue_list),3):
        node1=topRankingBetValue_list[i+1]
        node2=topRankingBetValue_list[i+2]
        try:
            node1_corrd=pro_stuc_chainRes[chain1][node1]["CA"].coord.tolist()
            node2_corrd=pro_stuc_chainRes[chain2][node2]["CA"].coord.tolist()
            #print(node1,node2,node1_corrd, node2_corrd)
        except:
            print("some mapping wrong",chain1,node1,chain2,node2)
            continue
        start_coords_list.append({"x":node1_corrd[0],"y":node1_corrd[1],"z":node1_corrd[2]})
        end_coords_list.append({"x":node2_corrd[0],"y":node2_corrd[1],"z":node2_corrd[2]})
    return(start_coords_list,end_coords_list)


def show_pdb(pdb_file,selected_chains,
             # start,end,
             hiden_chainList=None,
             colors=None,
             show_sidechains=False, show_mainchains=False,
            view_width=600,
            view_height=600):
    import py3Dmol
    # adapted from http://localhost:8206/lab/workspaces/auto-K/tree/code/MNF/src/tao_utilities/py3Dmol_functions.py show_pdb()
    
    #view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)
    # view = py3Dmol.view(js="https://cdn.jsdelivr.net/npm/3dmol@1.8.0/build/3Dmol-min.min.js") # https://github.com/3dmol/3Dmol.js/issues/635 and https://github.com/TaoDFang/MNF/issues/73
    #py3Dmol.view(query='mmtf:1ycr',js='https://cdn.jsdelivr.net/npm/3dmol@1.8.0/build/3Dmol-min.min.js'), download to /code/MNF/src/tao_utilities/3dmol%401.8.0%3Abuild%3A3Dmol-min.min.js?_xsrf=2%7Cda02d2d0%7Ca680e32e07e55cf9e801615f13db74f6%7C1672739124
    # view = py3Dmol.view(width=view_width, height=view_height,)
    view = py3Dmol.view(width=view_width, height=view_height,
                       js="https://cdn.jsdelivr.net/npm/3dmol@1.8.0/build/3Dmol-min.min.js")
    
    
    view.removeAllModels()
    view.addModel(open(pdb_file,'r').read(),'pdb')
    
    #view.getModel().hide()
    
    #https://pymolwiki.org/index.php/Color_Values
    if colors is None:
        colors=["purple	","blue","red","green","yellow","gray","black",][0:len(selected_chains)]
    
    for n,chain,color  in zip(range(len(selected_chains)),selected_chains,colors):
        print(f"chain {chain} use colorscheme {color}")
        view.setStyle({'chain':chain},{'cartoon': {'color':color}})
        #view.setStyle({'chain':chain},{'cartoon': {'opacity':0}})
        
    if hiden_chainList is not None:
        for n,chain in zip(range(len(hiden_chainList)),hiden_chainList):
            view.setStyle({'chain':chain},{'cartoon': {'opacity':0}})

    if show_sidechains:
        BB = ['C','O','N']
        view.addStyle({'and':[{'resn':["GLY","PRO"],'invert':True},{'atom':BB,'invert':True}]},
                        {'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
        view.addStyle({'and':[{'resn':"GLY"},{'atom':'CA'}]},
                        {'sphere':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
        view.addStyle({'and':[{'resn':"PRO"},{'atom':['C','O'],'invert':True}]},
                        {'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})  
    if show_mainchains:
        BB = ['C','O','N','CA']
        view.addStyle({'atom':BB},{'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
        
    # show residue label by javascript functions : 
    # https://colab.research.google.com/drive/1T2zR59TXyWRcNxRgOAiqVPJWhep83NV_?usp=sharing#scrollTo=cwGY0eGdGx5z
    # It isn't possible to convert Python functions to Javascript functions, but Javascript code can be provided in string form to click/hover callbacks.
    view.setHoverable({},True,'''function(atom,viewer,event,container) {
                   if(!atom.label) {
                    atom.label = viewer.addLabel(atom.resn+":"+atom.atom,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                   }}''',
               '''function(atom,viewer) { 
                   if(atom.label) {
                    viewer.removeLabel(atom.label);
                    delete atom.label;
                   }
                }''')

    # add lable by coordinate  https://github.com/3dmol/3Dmol.js/issues/391
    # viewer.addLabel("Hydrogen Donor", {position: {x:-5.4217, y:-0.4795, z:0.6395}, backgroundColor: 0x808080, backgroundOpacity:0.8});
    # or just add in other software like powerpoint or affinity design

        
    view.zoomTo()
    # view.render()
    # view.zoom(0.8, 2000)
    return view