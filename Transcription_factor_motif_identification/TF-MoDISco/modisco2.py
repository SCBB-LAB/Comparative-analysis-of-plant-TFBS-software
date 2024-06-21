import deeplift
import keras
from collections import OrderedDict
from deeplift.conversion import kerasapi_conversion as kc
from deeplift.layers import NonlinearMxtsMode  
from deeplift.util import get_shuffle_seq_ref_function
from deeplift.dinuc_shuffle import dinuc_shuffle
from deeplift.util import get_hypothetical_contribs_func_onehot
import modisco
import modisco.util
from modisco.visualization import viz_sequence
import modisco.affinitymat.core
import modisco.cluster.phenograph.core
import modisco.cluster.phenograph.cluster
import modisco.cluster.core
import modisco.aggregator
from collections import Counter
import os
import sys
import time
import h5py
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

'''
sample_path: the path of samples (fasta)
species: one of three options (Arabidopsis_models/Oryza_sativa_models/Zea_mays_models)
tf: the name of TF

'''

def one_hot_encode_along_channel_axis(sequence):
    to_return = np.zeros((len(sequence),4), dtype=np.int8)
    seq_to_one_hot_fill_in_array(zeros_array=to_return,
                                 sequence=sequence, one_hot_axis=1)
    return to_return

def seq_to_one_hot_fill_in_array(zeros_array, sequence, one_hot_axis):
    assert one_hot_axis==0 or one_hot_axis==1
    if (one_hot_axis==0):
        assert zeros_array.shape[1] == len(sequence)
    elif (one_hot_axis==1): 
        assert zeros_array.shape[0] == len(sequence)
    #will mutate zeros_array
    for (i,char) in enumerate(sequence):
        if (char=="A" or char=="a"):
            char_idx = 0
        elif (char=="C" or char=="c"):
            char_idx = 1
        elif (char=="G" or char=="g"):
            char_idx = 2
        elif (char=="T" or char=="t"):
            char_idx = 3
        elif (char=="N" or char=="n"):
            continue #leave that pos as all 0's
        else:
            raise RuntimeError("Unsupported character: "+str(char))
        if (one_hot_axis==0):
            zeros_array[char_idx,i] = 1
        elif (one_hot_axis==1):
            zeros_array[i,char_idx] = 1

def modisco_test(sample_path,species,tf):
    with open(sample_path) as f:
        data = []
        while True:
            line = f.readline()
            if not line:
                break
            if not line[0] == '>':
                data.append(line.strip()) 
    f.close()
    #densenet_path='/home/hlcheng/corn/desenet_without_parameter/64_64/bHLH145/checkmodel.hdf5'
    model_path=species+'/'+tf+'_checkmodel.hdf5'
    print(model_path)
    deeplift_model =\
        kc.convert_model_from_saved_files(
            model_path,
            nonlinear_mxts_mode=deeplift.layers.NonlinearMxtsMode.DeepLIFT_GenomicsDefault)
    l=list(deeplift_model.get_name_to_layer().keys())
    
def modisco_test(sample_path,species,tf):
    with open(sample_path) as f:
        data = []
        while True:
            line = f.readline()
            if not line:
                break
            if not line[0] == '>':
                data.append(line.strip()) 
    f.close()
    #densenet_path='/home/hlcheng/corn/desenet_without_parameter/64_64/bHLH145/checkmodel.hdf5'
    model_path=species+'/'+tf+'_checkmodel.hdf5'
   # print(model_path)
    deeplift_model =\
        kc.convert_model_from_saved_files(
            model_path,
            nonlinear_mxts_mode=deeplift.layers.NonlinearMxtsMode.DeepLIFT_GenomicsDefault)
    l=list(deeplift_model.get_name_to_layer().keys())

    deeplift_contribs_func = deeplift_model.get_target_contribs_func(
                                find_scores_layer_name=l[0],
                                pre_activation_target_layer_name=l[-2])

    rescale_conv_revealcancel_fc_many_refs_func = get_shuffle_seq_ref_function(
        #score_computation_function is the original function to compute scores
        score_computation_function=deeplift_contribs_func,
        #shuffle_func is the function that shuffles the sequence
        #technically, given the background of this simulation, randomly_shuffle_seq
        #makes more sense. However, on real data, a dinuc shuffle is advisable due to
        #the strong bias against CG dinucleotides
        shuffle_func=dinuc_shuffle,
        one_hot_func=lambda x: np.array([one_hot_encode_along_channel_axis(seq) for seq in x]))

    num_refs_per_seq=10 #number of references to generate per sequence
    multipliers_func = deeplift_model.get_target_multipliers_func(find_scores_layer_name=l[0],
                            pre_activation_target_layer_name=l[-2])
    hypothetical_contribs_func = get_hypothetical_contribs_func_onehot(multipliers_func)

    #Once again, we rely on multiple shuffled references
    hypothetical_contribs_many_refs_func = get_shuffle_seq_ref_function(
        score_computation_function=hypothetical_contribs_func,
        shuffle_func=dinuc_shuffle,
        one_hot_func=lambda x: np.array([one_hot_encode_along_channel_axis(seq)
                                        for seq in x]))
    onehotdata=np.array([one_hot_encode_along_channel_axis(seq) for seq in data])
    task_to_contrib_scores = {}
    task_to_hyp_contrib_scores = {}
    all_tasks = [0]
    for task_idx in all_tasks:
        print("On task",task_idx)
        task_to_contrib_scores[task_idx] =\
            np.sum(rescale_conv_revealcancel_fc_many_refs_func(
                task_idx=task_idx,
                input_data_sequences=data,
                num_refs_per_seq=num_refs_per_seq,
                batch_size=200,
                progress_update=1000,
            ),axis=2)[:,:,None]*onehotdata
        task_to_hyp_contrib_scores[task_idx] =\
            hypothetical_contribs_many_refs_func(
                task_idx=task_idx,
                input_data_sequences=data,
                num_refs_per_seq=num_refs_per_seq,
                batch_size=200,
                progress_update=1000,
            )  
    for task_idx in all_tasks:
        max_diff = np.max(np.abs(task_to_hyp_contrib_scores[task_idx]*onehotdata - task_to_contrib_scores[task_idx]))
        print("task",task_idx,"max diff:",max_diff)
        assert max_diff < 1e-2 #assert the difference is within numerical precision

    if not os.path.exists('./%s_modisco'%(tf)):
        os.mkdir('./%s_modisco'%(tf))
    result_dir='./%s_modisco'%(tf)
    score_path=result_dir+'/scores.h5'
    f = h5py.File(score_path,'w')
    g = f.create_group("contrib_scores")
    for task_idx in all_tasks:
        g.create_dataset("task"+str(task_idx), data=task_to_contrib_scores[task_idx])
    g = f.create_group("hyp_contrib_scores")
    for task_idx in all_tasks:
        g.create_dataset("task"+str(task_idx), data=task_to_hyp_contrib_scores[task_idx])
    f.close()
    task_to_scores = OrderedDict()
    task_to_hyp_scores = OrderedDict()
    f = h5py.File(score_path,"r")
    tasks = f["contrib_scores"].keys() # 'task0'
    n =len(data)  #since this is just a test run, for speed I am limiting to 100 sequences
    for task in tasks:
        task_to_scores[task] = [np.array(x) for x in f['contrib_scores'][task][:n]]
        task_to_hyp_scores[task] = [np.array(x) for x in f['hyp_contrib_scores'][task][:n]]
        print(len(task_to_scores[task]))
    onehot_data = [one_hot_encode_along_channel_axis(seq) for seq in data ][:n]


    #TF-MoDISco workflow
    null_per_pos_scores = modisco.coordproducers.LaplaceNullDist(num_to_samp=5000)
    tfmodisco_results = modisco.tfmodisco_workflow.workflow.TfModiscoWorkflow(
                        #Slight modifications from the default settings
                        sliding_window_size=15,
                        flank_size=5,
                        target_seqlet_fdr=0.15,
                        seqlets_to_patterns_factory=
                        modisco.tfmodisco_workflow.seqlets_to_patterns.TfModiscoSeqletsToPatternsFactory(
                            #Note: as of version 0.5.6.0, it's possible to use the results of a motif discovery
                            # software like MEME to improve the TF-MoDISco clustering. To use the meme-based
                            # initialization, you would specify the initclusterer_factory as shown in the
                            # commented-out code below:
                            #initclusterer_factory=modisco.clusterinit.memeinit.MemeInitClustererFactory(    
                            #    meme_command="meme", base_outdir="meme_out",            
                            #    max_num_seqlets_to_use=10000, nmotifs=10, n_jobs=1),
                            trim_to_window_size=15,
                            initial_flank_to_add=5,
                            kmer_len=5, num_gaps=1,
                            num_mismatches=0,
                            final_min_cluster_size=60)
                    )(
                    task_names=["task0"],
                    contrib_scores=task_to_scores,
                    hypothetical_contribs= task_to_hyp_scores,
                    one_hot=onehot_data,
                    null_per_pos_scores = null_per_pos_scores)
    results_path=result_dir+'/results.hdf5'
    grp = h5py.File(results_path, 'w')
    tfmodisco_results.save_hdf5(grp)
    grp.close()
    # get pwm 
    
    hdf5_results = h5py.File(results_path,'r')
    print("Metaclusters heatmap")
    
    activity_patterns = np.array(hdf5_results['metaclustering_results']['attribute_vectors'])[
                        np.array(
            [x[0] for x in sorted(
                    enumerate(hdf5_results['metaclustering_results']['metacluster_indices']),
                key=lambda x: x[1])])]
    sns.heatmap(activity_patterns, center=0)
   # plt.show()

    metacluster_names = [
        x.decode("utf-8") for x in 
        list(hdf5_results["metaclustering_results"]
            ["all_metacluster_names"][:])]

    all_patterns = []
    background = np.mean(onehot_data, axis=(0,1))

    all_pwms = []
    for metacluster_name in metacluster_names:
        print(metacluster_name)
        metacluster_grp = (hdf5_results["metacluster_idx_to_submetacluster_results"]
                                    [metacluster_name])
        print("activity pattern:",metacluster_grp["activity_pattern"][:])
        all_pattern_names = [x.decode("utf-8") for x in 
                            list(metacluster_grp["seqlets_to_patterns_result"]
                                                ["patterns"]["all_pattern_names"][:])]
        if (len(all_pattern_names)==0):
            print("No motifs found for this activity pattern")
        for pattern_name in all_pattern_names:
            print(metacluster_name, pattern_name)
            all_patterns.append((metacluster_name, pattern_name))
            pattern = metacluster_grp["seqlets_to_patterns_result"]["patterns"][pattern_name]
            print("total seqlets:",len(pattern["seqlets_and_alnmts"]["seqlets"]))
            print("Task 0 hypothetical scores:")
            #viz_sequence.plot_weights(pattern["task0_hypothetical_contribs"]["fwd"])
            print("Task 0 actual importance scores:")
            #viz_sequence.plot_weights(pattern["task0_contrib_scores"]["fwd"])
            print("onehot, fwd and rev:")
            #viz_sequence.plot_weights(pattern["sequence"]["fwd"])
            #viz_sequence.plot_weights(pattern["sequence"]["rev"])
            #viz_sequence.plot_weights(viz_sequence.ic_scale(np.array(pattern["sequence"]["fwd"]),
                                                            #background=background))
            #viz_sequence.plot_weights(viz_sequence.ic_scale(np.array(pattern["sequence"]["rev"]),
                                                            #background=background)) 
            all_pwms.append(np.array(pattern["sequence"]["fwd"]))
            #Plot the subclustering too, if available
            
    hdf5_results.close()


    def NUMPY2STRING(input_array):
        # convert numpy to string for 2 dimension numpy array.
        output_str = ""
        
        
        for i in range(input_array.shape[0]):
            for j in range(input_array.shape[1]):                
                output_str = output_str + str(input_array[i, j]) + "\t"
            output_str += "\n"
#            print(output_str)
         
        return output_str
    
    pwmmeme_result=result_dir+'/pwm.meme'
    with open( pwmmeme_result,'w') as f:
        f.write("MEME version 5.0.4\n\n")
        f.write("ALPHABET= ACGT\n\n")
        f.write("strands: + -\n\n")
        f.write("Background letter frequencies\n")
        f.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")
        for i, pwm in enumerate(all_pwms):
            f3 = open(result_dir+"/pwm_"+str(i+1)+".txt", 'w')
            one = "MEME version 5.0.4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
            f.write("MOTIF" + "\t" + "filter" + str(i) + "\n")
            f.write('letter-probability matrix: alength= 4 w=%s'%pwm[3:22,:].shape[0]+'\n')            
            outline=NUMPY2STRING(pwm[3:22,:])
            f.write(outline+'\n')
            
            f3.writelines(one+"MOTIF" + "\t" + "filter" + str(i) + "\n")
            f3.writelines('letter-probability matrix: alength= 4 w=%s'%pwm[3:22,:].shape[0]+'\n')            
            f3.writelines(outline+'\n')            
            f3.close()
#            os.system("python3.8 m5.py "+ "pwm_" + str(i) + " "+ sys.argv[1])
            os.system("python3.8 x6.py "+ result_dir+"/pwm_"+str(i+1) + " "+ sys.argv[1] + " temp") 
    f.close()

#os.system("python3.8 m5.py pwm_4 " + sys.argv[1])


def main():
    sample_path=sys.argv[1]
    species=su=sys.argv[2]
    tf=sys.argv[3]
    modisco_test(sample_path,species,tf)
if __name__=='__main__':
    main()
