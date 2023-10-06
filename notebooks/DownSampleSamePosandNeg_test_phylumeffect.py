import subprocess



##################################################################
newSTRING_rootFolder="/mnt/mnemo6/tao/PPI_Coevolution/STRING_data_11.5/"
#phylum2spe_list=[('1236', '511145')] 
# phylum2spe_list=[('1224', '511145')] 
phylum2spe_list=[('1236', '511145'), ('1224', '511145'),("2", "511145"),] 


DownSample_strategy = "num"#"percent"
#DownSample_sizes=[3,10,50,100,200,300]
# DownSample_sizes=[300]
DownSample_sizes=[5]
for DownSample_size in DownSample_sizes:
    for current_EggNOG_maxLevel,currentSpe_TaxID in phylum2spe_list: 
        print(current_EggNOG_maxLevel,currentSpe_TaxID)
        cmd = [
            "python", 
            "script_DownSampleSamePosandNeg_test_phylumeffect.py",
            '-l',current_EggNOG_maxLevel,
            '-i', currentSpe_TaxID,
            '-s',DownSample_strategy,
            '-m',str(DownSample_size), # because subprocess only accpet string 
        ]


        #skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        log_name = newSTRING_rootFolder+currentSpe_TaxID+"DownSample_"+DownSample_strategy+str(DownSample_size)+"_DownsampleSamePosandNeg_EggNOGmaxLevel"+current_EggNOG_maxLevel+".log"
        print(log_name)
        with open(log_name, 'w') as logfp:
            skw = dict(stdout=logfp, stderr=logfp)
            p = subprocess.run(cmd, universal_newlines=True, **skw)

        print("p.returncode",p.returncode)
        

# ###################################################################
# newSTRING_rootFolder="/mnt/mnemo6/tao/PPI_Coevolution/STRING_data_11.5/"
# phylum2spe_list=[("2", "511145"),] 


# DownSample_strategy = "num"#"percent"
# #DownSample_sizes=[10,50,100,200,300]
# DownSample_sizes=[3,500]
# for DownSample_size in DownSample_sizes:
#     for current_EggNOG_maxLevel,currentSpe_TaxID in phylum2spe_list: 
#         print(current_EggNOG_maxLevel,currentSpe_TaxID)
#         cmd = [
#             "python", 
#             "script_DownSampleSamePosandNeg_test_phylumeffect.py",
#             '-l',current_EggNOG_maxLevel,
#             '-i', currentSpe_TaxID,
#             '-s',DownSample_strategy,
#             '-m',str(DownSample_size), # because subprocess only accpet string 
#         ]


#         #skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)

#         log_name = newSTRING_rootFolder+currentSpe_TaxID+"DownSample_"+DownSample_strategy+str(DownSample_size)+"_DownsampleSamePosandNeg_EggNOGmaxLevel"+current_EggNOG_maxLevel+".log"
#         print(log_name)
#         with open(log_name, 'w') as logfp:
#             skw = dict(stdout=logfp, stderr=logfp)
#             p = subprocess.run(cmd, universal_newlines=True, **skw)

#         print("p.returncode",p.returncode)

        
        
        
        
        
###################################################################
# newSTRING_rootFolder="/mnt/mnemo6/tao/PPI_Coevolution/STRING_data_11.5/"


# phylum2spe_list=[("1224", "511145"),("2", "511145"),("1236", "511145")] 


# #RUN ONCE

# for current_EggNOG_maxLevel,currentSpe_TaxID in phylum2spe_list: 
#     print(current_EggNOG_maxLevel,currentSpe_TaxID)
#     cmd = [
#         "python", 
#         "script_DownSampleSamePosandNeg_test_phylumeffect.py",
#         '-l',current_EggNOG_maxLevel,
#         '-i', currentSpe_TaxID,
#     ]


#     #skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)

#     log_name = newSTRING_rootFolder+currentSpe_TaxID+"_DownsampleSamePosandNeg_EggNOGmaxLevel"+current_EggNOG_maxLevel+".log"
#     print(log_name)
#     with open(log_name, 'w') as logfp:
#         skw = dict(stdout=logfp, stderr=logfp)
#         p = subprocess.run(cmd, universal_newlines=True, **skw)
        
#     print("p.returncode",p.returncode)


