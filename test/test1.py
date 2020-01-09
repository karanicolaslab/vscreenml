#!/usr/bin/python
import os,re, csv, subprocess, sys

import os
import pickle
import pandas as pd
import numpy as np
import xgboost as xgb
from xgboost.sklearn import XGBClassifier

VSCREENML = os.getenv('VSCREENML')

data = np.array([-32.6344,7.80936,2.76329,-0.938994,0,0,-27.8004,842.605,0,0.723475,6,0.333019,0.840562,0.159438,3534,871,1000,57,478,117,138,8,344,80,103,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.7,0.0,0.0,0.0,0.0]).reshape((1, -1))

predictors = ['fa_atr','fa_rep','fa_sol','fa_elec','hbond_bb_sc','hbond_sc','Interface_Energy','Total_BSA','Interface_HB','Total_packstats','Interface_unsat','Total_pose_exposed_SASA','interface_hydrophobic_sasa','interface_polar_sasa','6.6','7.6','8.6','16.6','6.7','7.7','8.7','16.7','6.8','7.8','8.8','16.8','6.9','7.9','8.9','16.9','6.15','7.15','8.15','16.15','6.16','7.16','8.16','16.16','6.17','7.17','8.17','16.17','6.35','7.35','8.35','16.35','6.53','7.53','8.53','16.53','SIDE_flex_ALPHA','SIDE_flex_BETA','SIDE_flex_OTHER','BACK_flex_ALPHA','BACK_flex_BETA','BACK_flex_OTHER','TotalElec','TotalHbond','Hphobe_contact','Pi_Pi','T-stack','Cation-pi','Salt_Bridge','diff_entropy','pienergy','fsp3','polarsurfacearea','vdwsa']


classes = ['decoy', 'active']
def predict(model, tag, input_arr):
    print (input_arr)
    print (input_arr.shape)
    prediction = model.predict(input_arr)
    prediction_prob = model.predict_proba(input_arr)[:,1][0]
    print('{}\t{}\t{}'.format('Name','vscreenml_score','Label'))
    print('{}\t{}\t{}'.format(tag, prediction_prob, classes[int(prediction)]))


with open(VSCREENML+'/model/XGB_CLASSIFIER_alldata.pickle.dat', 'rb') as f:
    clf = pickle.load(f)


df = pd.DataFrame(data=data, columns=predictors)
predict(clf, 'filetag', df)

