#oddt
###
import os
import sys
import pickle
import pandas as pd
import numpy as np
from numpy import linspace
import xgboost as xgb
import seaborn as sns
from xgboost.sklearn import XGBClassifier
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
import matplotlib.style as style
from sklearn.metrics import confusion_matrix
from numpy.random import seed
from numpy.random import randn
from scipy.stats import wilcoxon
import warnings


def predict_rank_active_other_ppi(test_x,score,name):

    no_of_points =float(len(test_x))
    if score != 'vina_affinity':
        sort_final = test_x.sort_values(by=score, ascending=False)
    else:
        sort_final = test_x.sort_values(by=score, ascending=True)
    rank=list(sort_final["Label"]).index(1) +1
    #print "Rank for", name , "=", rank
    print ('Rank for the active native ligand {} is {}'.format(name, rank+1))
    #print sort_final.head()
    return rank

#rfscorevs
def predict_rank_active_rf_ppi(test_x,name):

    no_of_points =float(len(test_x))

    sort_final = test_x.sort_values(by='Score', ascending=False)
    rank=list(sort_final["Label"]).index(1)+1
    #print "Rank for", name , "=", rank
    print ('Rank for the active native ligand in {} is {}'.format(name, rank+1))
    #print sort_final.head()
    return rank

#vscreenml
def predict_rank_vscreenml_ppi(model,test_x,predictors,name):
    actual_positive=float(len(test_x[test_x["Label"] == 1]))

    no_of_points =float(len(test_x))
    
    random_hit_rate = float(actual_positive/no_of_points) *100
    #print actual_positive
    dtest_predict = model.predict(test_x[predictors])
    dtest_prob = model.predict_proba(test_x[predictors])[:,1]
    test_x_pred = test_x.copy()
    test_x_pred['prediction'] = dtest_predict
    test_x_pred['probability'] = dtest_prob
    
    test_x_final1 = pd.DataFrame()
    test_x_final1['Name'] = test_x_pred['Name']
    test_x_final1['Score'] = test_x_pred['probability']
    test_x_final1['Prediction'] = test_x_pred['prediction']
    test_x_final1['Label'] = test_x_pred["Label"]
    sort_final = test_x_final1.sort_values(by='Score', ascending=False)
    
    rank=list(sort_final["Label"]).index(1)+1
    #print "Rank for", name , "=", rank
    print ('Rank for the active native ligand in {} is {}'.format(name, rank+1))
    #print sort_final.head()
    return rank

#rosetta
def predict_rank_active_rosetta_ppi(test_x,score,name,order=False):

    no_of_points =float(len(test_x))

    sort_final = test_x.sort_values(by=score, ascending=order)
    rank=list(sort_final["Label"]).index(1) +1
    print ('Rank for the active native ligand in {} is {}'.format(name, rank+1))
    #print sort_final.head()
    return rank

def predict_from_rf_otherscore_dekoise(test_x,score,name):
    actual_positive=float(len(test_x[test_x["Label"] == 1]))
    
    no_of_points =float(len(test_x))
    
    random_hit_rate = float(actual_positive/no_of_points) *100

    sort_final = test_x.sort_values(by=score, ascending=False)
    one_percent = int (len(sort_final['Name'].tolist()) * 0.01)
    #print one_percent
    #Print model report:
    #print "\nPrediction Report"
    #print "Accuracy(Test) : %.4g" % metrics.accuracy_score(test_x['label'].values, dtest_predict)
    #print "AUC Score (Test): %f" % metrics.roc_auc_score(test_x['label'], dtest_prob)
    #print len(sort_final['PDB'].tolist()),"complexes are active"
    one_percent_df = sort_final.head(one_percent)
    #print one_percent_df
    #true_positive=float(len(one_percent_df[(one_percent_df[score] > 6.0) & (one_percent_df.Label == 1)]))
    true_positive=float(len(one_percent_df[(one_percent_df.Label == 1)]))
    #print true_positive
    hit_rate = true_positive/one_percent*100
    
    enrichment = hit_rate/random_hit_rate
    
    print ("Enrichment for", name , "=", enrichment)
    return enrichment

def predict_1percent_enrichment_dekois(model, test_x, predictors, scores):
    
    random_hit_rate = test_x['Label'].sum()/float(len(test_x['Label']))*100

    
    dtest_pred = model.predict_proba(test_x[predictors])[:,1]
    test_x["vScreenML"] = dtest_pred

    scores.extend(["vScreenML","Interface_Energy"])
    

    
    enrichment_dict = {}
    for score in scores:
        #print (y_test.dtypes, test_x[score].dtypes)
                    
        if score == 'vina_affinity' or score == 'Interface_Energy':
                sort_final = test_x.sort_values(by=score, ascending=True)
                #print (sort_final.head())
                one_percent = int (len(sort_final['Name'].tolist()) * 0.01)
                one_percent_df = sort_final.head(one_percent)
                true_positive=float(len(one_percent_df[(one_percent_df.Label == 1)]))
                hit_rate = true_positive/one_percent*100    
                enrichment = hit_rate/random_hit_rate
                enrichment_dict[score] = enrichment
        else:
                sort_final = test_x.sort_values(by=score, ascending=False)
                one_percent = int (len(sort_final['Name'].tolist()) * 0.01)
                one_percent_df = sort_final.head(one_percent)
                true_positive=float(len(one_percent_df[(one_percent_df.Label == 1)]))
                hit_rate = true_positive/one_percent*100    
                enrichment = hit_rate/random_hit_rate
                enrichment_dict[score] = enrichment
                
    return enrichment_dict

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
    plt.xlim((-20,2000)) 
    plt.ylim((-20,2000))
    plt.xlabel(column1,fontsize=18,fontweight='bold')
    plt.ylabel(column2,fontsize=18,fontweight='bold')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.title(column2+' vs '+column2)
    plt.show()
