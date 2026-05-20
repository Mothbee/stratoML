import sys
import node
import tree_reader,read_fasta,tree_utils,stratlike,mfc
import numpy as np
import qmat
from scipy.optimize import minimize
import time


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("usage: "+ sys.argv[0]+ " <newick> <trait fasta file> <stratigraphic data> <stratigraphic model> <morphologic model> <geo or nogeo>")
        sys.exit()

    geo = True if sys.argv[6] == 'geo' else False

    traits,ss = read_fasta.read_fasta(sys.argv[2]) # all traits together
    retraits  = read_fasta.recode_poly_traits(traits,ss)

    for line in open(sys.argv[1],"r"):
        nwk = line.strip().split()[-1]
        tree = tree_reader.read_tree_string(nwk)
       # for n in tree.iternodes():
       #     print(n.label)

        tree_utils.map_strat_to_tree(tree,sys.argv[3])  ###
        #stratlike.calibrate_brlens_strat(tree,0.3)
        #print(tree.get_newick_repr(True))
        tree_utils.map_tree_disc_traits(tree,retraits,ss,geo) # traits mapped HERE
        tree_utils.fix_obs_lv(tree,geo)  ### add geo
        ###tree_utils.fix_obs_lvgeo(tree)  ### add geo

        #tree_utils.sort_children_by_age(tree)
        #tree_utils.init_budd_marginals(tree,len(ss))
        qmats = qmat.Qmat(0.01,0.05) # create the qmatrix with the provided m and l
        if geo == True:
            qmatsgeo = qmat.Qmat(0.01,0.05) # diff default m and l? rn exactly same as qmats
        #print(qmatsgeo, qmats)
        ### qmat agnostic to chars, just abt n states and rates


        for n in tree.iternodes():
            if geo == True:
                n.update_pmat(qmats,max(ss[0:-1]),"mid") # get p matrix for each node
                n.update_pmatgeo(qmatsgeo,ss[-1],"mid") # for geo
            else:
                n.update_pmat(qmats,max(ss),"mid") # get p matrix for each node
            #n.update_pmat(qmats,max(ss),"mid") # get p matrix for each node
            #n.update_pmatgeo(qmats,max(ss),"mid")
        ### transition probabilities depend on transition rates, n char states, and taxon duration


        #pmat = qmats.calc_single_p_mat(1.0, 2)
        #for i in pmat:
        #    print(list(i))

        #print(tree)
        #treell = -mfc.evaluate_m_l2(np.array([0.01,0.05]),tree,qmats,ss)
        #print("TREELL 1", treell) 
        #treell = -mfc.evaluate_m_l2(np.array([0.01,0.001]),tree,qmats,ss)
        #print("TREELL 2", treell)
        t1 = time.time()
        if geo == True:
            ###aic,traitll,bdsll,geoll = tree_utils.calc_tree_ll2geo(tree,qmats,qmatsgeo,ss,"hr97") ### added geoll here
            tree_utils.calc_tree_ll2geo(tree,qmats,qmatsgeo,ss,"hr97")
        else:  
            aic,traitll,bdsll = tree_utils.calc_tree_ll2(tree,qmats,ss,"hr97") 
        #   aic,traitll,bdsll = tree_utils.calc_tree_ll2perchar(tree,qmats,ss,"hr97") ### spit out char by char
            
        ### ADDED TO SPIT OUT SUBLIKES
            treell,sublikes = mfc.mfc2_treellperchar(tree,qmats,ss) ### corrected MFC treell, char by char likelihoods
            print("AIC: ", aic, "TRAITLL: ", traitll, "STRATLL: ", bdsll)
            sublikes.insert(0, sys.argv[1]) ### name of tree and sublikes
           # print(','.join(str(x) for x in sublikes)) ### add >> <flnm.txt> at end of command to append to csv

            

        t2 = time.time()
        #print(aic,traitll,bdsll)
        ###print(aic,nwk)
        ###print("TIME OPTIMIZING",t2-t1) 

   
    
        '''if geo == True: ### COMMENT THESE OUT TO DO 1 ONLY
            tree_utils.tree_search3geo(tree,ss,qmats,qmatsgeo,"hr97",False) ### either 1 function with all or add a separate one for geo 
        else:
            tree_utils.tree_search3(tree,ss,qmats,"hr97",False)'''
    

        


"""
        res_st = minimize(stratlike.poisson_neg_ll,x0=np.array([1.0]),args=(tree),method="Nelder-Mead")
        bdsll = -res_st.fun
        print("stratlike:",bdsll)
        nparam = 1.0 + 2.0
        res_tr = minimize(mfc.evaluate_m_l2,x0=np.array([0.01,0.01]),args=(tree,qmats,ss),method="L-BFGS-B",bounds=((0.00001,0.5),(0.00001,0.5)))
        print("mfc rates",res_tr.x)
        traitll = -res_tr.fun
        print(traitll, len([n for n in tree.iternodes()]))
        tree_ll = traitll + bdsll

        nparam += float(len([n for n in tree.iternodes()]) - 1)
        aic = (2. * nparam) - (2. * tree_ll) 
        print("AIC",aic)
        #aic = tree_utils.single_tree_aic(tree,ss,sys.argv[4],sys.argv[5])
        #print(aic,tree.get_newick_repr()+";")
"""

"""
for n in tree.iternodes():
    for desc in n.timeslice_lv:
        for char in desc:
            print(list(char))
"""

#print(dir(tree)) ###