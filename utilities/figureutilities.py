#oddt
###
import numpy as np
from numpy import linspace
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import metrics
from scipy.stats import wilcoxon


def wilcoxon_signed(data1, data2):

    stat, p = wilcoxon(data1, data2)
    print('Statistics=%.3f, p=%.5f' % (stat, p))
    # interpret
    alpha = 0.05
    if p > alpha:
        print('Same distribution (fail to reject H0)')
    else:
        print('Different distribution (reject H0)')


def plot_benchmark_graph_ppi(data, column1, column2):
    sns.regplot(x=column1, y=column2, data=data, fit_reg=False,marker='o', scatter_kws={"color":"purple",'s':80,"alpha":0.7});
    plt.plot([0,2000], [0,2000], color='g')
    plt.xlim((0,2000)) 
    plt.ylim((0,2000))
    plt.xlabel(column1,fontsize=18,fontweight='bold')
    plt.ylabel(column2,fontsize=18,fontweight='bold')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.title(column2+' vs '+column2)
    plt.show()

def plot_benchmark_graph_dekois(data, column1, column2):
    sns.regplot(x=column1, y=column2, data=data, fit_reg=False,marker='o', scatter_kws={"color":"purple",'s':80,"alpha":0.7});
    plt.plot([0,30], [0,30], color='g')
    plt.xlim((0,30))
    plt.ylim((0,30))
    plt.xlabel(column1,fontsize=18,fontweight='bold')
    plt.ylabel(column2,fontsize=18,fontweight='bold')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.title(column2+' vs '+column2)
    plt.show()

def predict_rank_and_enrichment(model_all,model_rf, test_x, predictors,predictors_rf, scores):
    
    random_hit_rate = test_x['Label'].sum()/float(len(test_x['Label']))*100

    
    dtest_pred = model_all.predict_proba(test_x[predictors])[:,1]
    test_x["vScreenML"] = dtest_pred
    
    dtest_pred_rf = model_rf.predict_proba(test_x[predictors_rf])[:,1]
    test_x["RFSCORE_V1_remade"] = dtest_pred_rf

    scores.extend(["vScreenML","Interface_Energy","RFSCORE_V1_remade"])
    

    
    enrichment_dict = {}
    rank_dict = {}
    for score in scores:
        #print (y_test.dtypes, test_x[score].dtypes)
                    
        if score == 'vina_affinity' or score == 'Interface_Energy':
                sort_final = test_x.sort_values(by=score, ascending=True)
                rank=list(sort_final["Label"]).index(1)+1
                #print (sort_final.head())
                one_percent = int (len(sort_final['Name'].tolist()) * 0.01)
                one_percent_df = sort_final.head(one_percent)
                true_positive=float(len(one_percent_df[(one_percent_df.Label == 1)]))
                hit_rate = true_positive/one_percent*100    
                enrichment = hit_rate/random_hit_rate
                enrichment_dict[score] = enrichment
                rank_dict[score] = rank
        else:
                sort_final = test_x.sort_values(by=score, ascending=False)
                rank=list(sort_final["Label"]).index(1)+1
                one_percent = int (len(sort_final['Name'].tolist()) * 0.01)
                one_percent_df = sort_final.head(one_percent)
                true_positive=float(len(one_percent_df[(one_percent_df.Label == 1)]))
                hit_rate = true_positive/one_percent*100    
                enrichment = hit_rate/random_hit_rate
                enrichment_dict[score] = enrichment
                rank_dict[score] = rank
                
    return enrichment_dict, rank_dict
