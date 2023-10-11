from collections import defaultdict

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