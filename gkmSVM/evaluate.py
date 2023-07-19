from sklearn.metrics import average_precision_score, roc_auc_score
import argparse
from imblearn.metrics import specificity_score
from sklearn.metrics import confusion_matrix, accuracy_score, matthews_corrcoef, f1_score  
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pos", type=str, help="gkmpredict output on positive test examples", required=True)
args = parser.parse_args()

def get_scores(fname):
    f = open(fname)
    d = [float(x.strip().split('\t')[1]) for x in f]
    f.close()
    return d

pos_scores = get_scores(args.pos+"_pos_gkmpredict.txt")
neg_scores = get_scores(args.pos+"_neg_gkmpredict.txt")
labels = [1]*len(pos_scores) + [0]*len(neg_scores)

auprc = average_precision_score(labels, pos_scores+neg_scores)
auroc = roc_auc_score(labels, pos_scores+neg_scores)

pred_labels = [int(p>=0) for p in pos_scores+neg_scores]

tn, fp, fn, tp = confusion_matrix(labels, pred_labels).ravel()
specificity = tn / (tn+fp)
specificity=specificity_score(labels, pred_labels)
sensitivity = tp / (tp+fn)
cm=confusion_matrix(labels, pred_labels)
MCC=matthews_corrcoef(labels, pred_labels)
acc=accuracy_score(labels, pred_labels)

f1=f1_score(labels, pred_labels)

f3 = open(args.pos+"_result.txt",'w')
f3.writelines("TP"+"\t"+"TN"+"\t"+"FN"+"\t"+"FP"+"\t"+"TPR"+"\t"+"TNR"+"\t"+"ACC"+"\t"+"F1"+"\t"+"MCC"+"\t"+"auc"+"\n")
f3.writelines(str(tp)+"\t"+str(tn)+"\t"+str(fn)+"\t"+str(fp)+"\t"+str(sensitivity)+"\t"+str(specificity)+"\t"+str(acc)+"\t"+str(f1)+"\t"+str(MCC)+"\t"+str(auroc)+"\n")

