import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import py3Dmol
from Bio.PDB import PDBParser
import sys

import xpdb


def show_pdb(pdb_file,chains,
             # start,end,
             show_sidechains=False, show_mainchains=False, color="lDDT",
            view_width=600,
            view_height=600,
            chainNames=[],
            chaincolors=[]):
    view = py3Dmol.view(width=view_width, height=view_height,
                       js="https://cdn.jsdelivr.net/npm/3dmol@1.8.0/build/3Dmol-min.min.js")
    
    view.removeAllModels()
    view.addModel(open(pdb_file,'r').read(),'pdb')

    if color == "lDDT":
        view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
    elif color == "rainbow":
        view.setStyle({'cartoon': {'color':'spectrum'}})
    elif color == "chain":
        for n,chain,color in zip(range(chains),list("ABCDEFGH"),
                     ["lime","cyan","magenta","yellow","salmon","white","blue","orange"]):
            view.setStyle({'chain':chain},{'cartoon': {'color':color}})
    elif color == "chain_withNames":
        for chain,color in zip(chainNames,chaincolors):
            view.setStyle({'chain':chain},{'cartoon': {'color':color}})  
    # check http://3dmol.csb.pitt.edu/doc/types.html#ColorschemeSpec for cllorscheme
    #http://3dmol.csb.pitt.edu/doc/$3Dmol.Gradient.ROYGB.html
    elif color == "chain_lDDT":
        for n,chain,gradient in zip(range(chains),list("ABCDEFGH"),
                     ["roygb","roygb","sinebow","sinebow","sinebow","sinebow","sinebow","sinebow"]):
            #remember here roygb is blue, 
            # sinebow is purple 
            #https://3dmol.csb.pitt.edu/doc/types.html#ColorschemeSpec
            # https://3dmol.csb.pitt.edu/doc/$3Dmol.Gradient.ROYGB.html
            #https://3dmol.csb.pitt.edu/doc/$3Dmol.Gradient.Sinebow.html
            print(f"chain {chain} use colorscheme {gradient}")
            view.setStyle({'chain':chain},{'cartoon':{'colorscheme': {'prop':'b','gradient': gradient,'min':50,'max':90}}},{"label":chain})
            
            
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
        

    # print(start,end)
    # print(  {"start":{"x":start["x"],"y":start["y"],"z":start["z"]},
    #          "end":{"x":start["x"],"y":start["y"],"z":start["z"]}}
    # )
    # view.addCylinder({"color":'red',
    #                   "radius":0.15,
    #                   "dashed":False,  # why this line "True" not working ?
    #                   "fromcamp":2,
    #                   "toCAP":2,
    #                   "start":{"x":start["x"],"y":start["y"],"z":start["z"]},
    #               "end":{"x":end["x"],"y":end["y"],"z":end["z"]}})
    
    view.zoomTo()
    # view.render()
    # view.zoom(0.8, 2000)
    return view


def py3Dmol_addOneLine(view,
                     start,
                     end,
                    color="red",
                       radius=0.15,
                     # start={"x":0,"y":0,"z":0},
                     # end={"x":0,"y":0,"z":0},
                    ):
    # add edges in between : 
    #view.addLine({color:'blue',start:{x:0,y:0,z:0},end:{x:0,y:5,z:0}});
    #print(start,end)
    view.addCylinder({"color":color,
                      "radius":radius,    #0.15,
                      "dashed":False, # True
                      "fromcamp":1,  #2,1
                       "toCAP":1,  #2,1
                      "start":{"x":start["x"],"y":start["y"],"z":start["z"]},
                  "end":{"x":end["x"],"y":end["y"],"z":end["z"]}})
    #view.show()
    
    return(view)

def py3Dmol_addLines(view,
                     start_list,
                     end_list,
                     color="red",
                     radius=0.15,
                    ):

    print(f"lines color is {color}")
    for i in range(len(start_list)):
        start=start_list[i]
        end=end_list[i]
        view=py3Dmol_addOneLine(view,start,end,color,radius)
    
    
    return(view)

def show_alphfold_plot(figure_fileName):

    img = mpimg.imread(figure_fileName)
    imgplot = plt.imshow(img)
    plt.show()

def get_pro_stuc_chainRes_dict(len1,len2,pdb_file):

    sloppy_parser = PDBParser(structure_builder=xpdb.SloppyStructureBuilder())
    pro_stuc = sloppy_parser.get_structure(id=None, file=pdb_file)

    pro_stuc_chainRes=dict()
    for chain in pro_stuc.get_chains():
        pro_stuc_chainRes[chain.get_id()]=[x for x in chain.get_residues()]
        print(chain.get_id(),len(pro_stuc_chainRes[chain.get_id()]))

    assert len1==len(pro_stuc_chainRes["B"])
    assert len2==len(pro_stuc_chainRes["C"])
    
    return(pro_stuc_chainRes)


def tripleMSA_get_pro_stuc_chainRes_dict(len1,len2,len3,pdb_file):

    sloppy_parser = PDBParser(structure_builder=xpdb.SloppyStructureBuilder())
    pro_stuc = sloppy_parser.get_structure(id=None, file=pdb_file)

    pro_stuc_chainRes=dict()
    for chain in pro_stuc.get_chains():
        pro_stuc_chainRes[chain.get_id()]=[x for x in chain.get_residues()]
        print(chain.get_id(),len(pro_stuc_chainRes[chain.get_id()]))

    assert len1==len(pro_stuc_chainRes["B"])
    assert len2==len(pro_stuc_chainRes["C"])
    assert len3==len(pro_stuc_chainRes["D"])
    
    return(pro_stuc_chainRes)

# this one actuall only works for alphafold pdb files 
def map_topRankingBetValue_list_to_PDB_CACoords(topRankingBetValue_list,pro_stuc_chainRes,
                                               chain1_ID="B",
                                               chain2_ID="C"):
    start_coords_list=list()
    end_coords_list=list()
    for i in range(0,len(topRankingBetValue_list),3):
        node1=topRankingBetValue_list[i+1]
        node2=topRankingBetValue_list[i+2]
        node1_corrd=pro_stuc_chainRes[chain1_ID][node1]["CA"].coord.tolist()
        node2_corrd=pro_stuc_chainRes[chain2_ID][node2]["CA"].coord.tolist()
        #print(node1,node2,node1_corrd, node2_corrd)
        start_coords_list.append({"x":node1_corrd[0],"y":node1_corrd[1],"z":node1_corrd[2]})
        end_coords_list.append({"x":node2_corrd[0],"y":node2_corrd[1],"z":node2_corrd[2]})
    return(start_coords_list,end_coords_list)
        
