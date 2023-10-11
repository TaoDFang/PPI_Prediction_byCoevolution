from collections import defaultdict
from Bio import AlignIO

#limited frquence of protein in protein pair list 
def get_PPIwithLimitedProFreByOr(pp_pairs,limitedProFre=30):

    limitedProFre_PP_pairs=list()
    
    limitation_countDict=defaultdict(int)
    #for p1,p2,l1,l2 in pp_pairs:
    for l in pp_pairs:
        p1,p2=l[0],l[1]
        #if (limitation_countDict[p1]<limitedProFre) and (limitation_countDict[p2]<limitedProFre):
        if (limitation_countDict[p1]<limitedProFre) or (limitation_countDict[p2]<limitedProFre):
            #limitedProFre_PP_pairs.append((p1,p2,l1,l2))
            limitedProFre_PP_pairs.append(l)
        limitation_countDict[p1] +=1
        limitation_countDict[p2] +=1

    return(limitedProFre_PP_pairs)



#heck paired MSA of homologus , if same proteins  from same species in both side
def getSameProteinRatio(record):
    pp1, pp2,final_msa_folder=record
    msa = AlignIO.read(final_msa_folder+pp1+"and"+pp2+".fasta", "fasta")

    sameProtein_count=0
    for record in msa:
        record_id=str(record.id)
        record_ids=record_id.split("and")
        if record_ids[0]==record_ids[1]:
            sameProtein_count +=1
    sameProtein_ratio=sameProtein_count/len(msa)
    #print(sameProtein_ratio)
    #sameProtein_ratio_dict[(pp1,pp2)]=sameProtein_ratio
    return([pp1,pp2,sameProtein_ratio])
