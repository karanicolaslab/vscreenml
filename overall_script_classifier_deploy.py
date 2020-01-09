#!/usr/bin/python
import os,re, csv, subprocess, sys

import os
import pickle
import pandas as pd
import numpy as np
import xgboost as xgb
from xgboost.sklearn import XGBClassifier


#get environment variable
BABEL = os.getenv('BABEL')
CHEMAXON = os.getenv('CHEMAXON')
ROSETTA_DIR = os.getenv('ROSETTA_DIR')
ROSETTA_BUILD = os.getenv('ROSETTA_BUILD')
MGLTOOLS = os.getenv('MGLTOOLS')
BINANA = os.getenv('BINANA')
OMEGA = os.getenv('OMEGA')
SZYBKI = os.getenv('SZYBKI')
CXCALC = os.getenv('CXCALC')
VSCREENML = os.getenv('VSCREENML')

if len (sys.argv) != 3 :
    print "Usage: python overall_classifier_deploy.py folderectory filetag_base_name "
    sys.exit (1)

#pass folderectory and filetag
folder =sys.argv[1]
filetag =sys.argv[2]

#Extract ligands and protein structures
p=subprocess.Popen('grep -E \"ATOM|ZN|MG|MN\" '+folder+'/mini_complex_'+filetag+'.pdb > ' +folder+'/mini_complex_'+filetag+'_unb.pdb && grep HETATM '+folder+'/mini_complex_'+filetag+'.pdb|grep -v \'ZN\'|grep -v \'MG\'|grep -v \'MN\' > '+folder+'/mini_complex_'+filetag+'_lig.pdb',shell=True)

p.wait()

#reformatting the ligands
p=subprocess.Popen(BABEL+'/bin/obabel '+folder+'/mini_complex_'+filetag+'_lig.pdb -O '+folder+'/mini_complex_'+filetag+'_lig1.sdf>/dev/null',shell=True)
p.wait()

p=subprocess.Popen(CHEMAXON+'/structure-representation-toolkit/bin/structurecheck -c \"valenceerror\" -m fix -f sdf -o '+folder+'/mini_complex_'+filetag+'_lig.sdf '+folder+'/mini_complex_'+filetag+'_lig1.sdf',shell=True)
p.wait()

#generating params filetag
p=subprocess.Popen(ROSETTA_DIR+'/main/source/scripts/python/public/molfile_to_params.py '+folder+'/mini_complex_'+filetag+'_lig1.sdf -p '+folder+'/'+filetag, shell=True)
p.wait()


#minimization of the complex
try:
    p_rosetta = subprocess.Popen(ROSETTA_DIR+'/main/source/bin/minimize_ppi.'+ROSETTA_BUILD+' -s '+folder+'/mini_complex_'+filetag+'.pdb -extra_res_fa '+folder+'/'+filetag+'.params -database '+ROSETTA_DIR+'/main/database -ignore_zero_occupancy false -score_only |grep \'Interface_Scores:mini\' |awk \'{print $4,$5,$6,$7,$8}\' ', shell=True, stdout= subprocess.PIPE)
    rosetta = p_rosetta.communicate()[0].split()
	
except:
    rosetta= ['FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL']

Interface_Energy = rosetta[0]
Total_BSA = rosetta[1]
Interface_HB = rosetta[2]
Total_packstats = rosetta[3]
Interface_unsat =rosetta[4]


#calculating sasa   
try:
    p_sasa = subprocess.Popen(ROSETTA_DIR+'/main/source/bin/ragul_get_ligand_sasa.'+ROSETTA_BUILD+' -database '+ROSETTA_DIR+'/main/database -input_bound_pdb '+folder+'/mini_complex_'+filetag+'.pdb -input_unbound_pdb '+folder+'/mini_complex_'+filetag+'_unb.pdb -input_ligand_pdb '+folder+'/mini_complex_'+filetag+'_lig.pdb -extra_res_fa '+folder+'/'+filetag+'.params |grep \'Scores:\' |awk \'{print $2,$3,$4}\'', shell=True, stdout= subprocess.PIPE)
    sasa = p_sasa.communicate()[0].split()

except:
    sasa = ['FAIL', 'FAIL', 'FAIL']

Total_pose_exposed_SASA = sasa[0]
interface_hydrophobic_sasa = sasa[1]
interface_polar_sasa = sasa[2]


#calculating componenet scores    
try:
    p_compo = subprocess.Popen(ROSETTA_DIR+'/main/source/bin/yusuf_interface_ddg.'+ROSETTA_BUILD+' -database '+ROSETTA_DIR+'/main/database -s '+folder+'/mini_complex_'+filetag+'.pdb -extra_res_fa '+folder+'/'+filetag+'.params -ignore_zero_occupancy false   |grep \'Component_score:\' |awk \'{print $4,$5,$6,$10,$11,$12}\'', shell=True, stdout= subprocess.PIPE)
    compo = p_compo.communicate()[0].split()

except:
    compo = ['FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL']

fa_atr = compo[0] 
fa_rep = compo[1]
fa_sol = compo[2]
fa_elec = compo[3]
hbond_bb_sc = compo[4]
hbond_sc= compo[5]


#generating .pdbqt filetags for BINANA features
p=subprocess.Popen(MGLTOOLS+'/bin/pythonsh '+MGLTOOLS+'/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r '+folder+'/mini_complex_'+filetag+'_unb.pdb -A hydrogens -o '+folder+'/mini_complex_'+filetag+'_unb.pdbqt',shell=True)
p.wait()

p=subprocess.Popen(MGLTOOLS+'/bin/pythonsh '+MGLTOOLS+'/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l '+folder+'/mini_complex_'+filetag+'_lig.pdb -A hydrogens -o '+folder+'/mini_complex_'+filetag+'_lig.pdbqt',shell=True)
p.wait()

#calculating BINANA features
try:
    p=subprocess.Popen('python '+BINANA+'/binana_1_2_0.py -receptor '+folder+'/mini_complex_'+filetag+'_unb.pdbqt -ligand '+folder+'/mini_complex_'+filetag+'_lig.pdbqt > '+folder+'/mini_complex_'+filetag+'_binana.out',shell=True)
    p.wait()
    p_binana = subprocess.Popen('python '+BINANA+'/parse_binana.py '+folder+'/mini_complex_'+filetag+'_binana.out', shell=True, stdout=subprocess.PIPE)
    binana = p_binana.communicate()[0].split()
except:
    binana = ['FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL']
SIDE_flex_ALPHA = binana[0]
SIDE_flex_BETA= binana[1]
SIDE_flex_OTHER=binana[2]
BACK_flex_ALPHA= binana[3]
BACK_flex_BETA= binana[4]
BACK_flex_OTHER= binana[5]
TotalElec= binana[6]
TotalHbond= binana[7]
Hphobe_contact= binana[8]
Pi_Pi= binana[9]
Tstack= binana[10]
Cationpi= binana[11]
Salt_Bridge= binana[12]

#generating OMEGA conformers for entropy calculations
p=subprocess.Popen(OMEGA+' -in '+folder+'/mini_complex_'+filetag+'_lig1.sdf -out '+folder+'/'+filetag+'.oeb.gz -strictstereo false -strictatomtyping false -warts true -maxconfs 1000 -rms 0.1',shell=True)
p.wait()

try:
    p = subprocess.Popen(SZYBKI+' -entropy AN  -prefix '+ folder+'/bou_'+filetag+' -complex '+folder+'/mini_complex_'+filetag+'.pdb -strict false 2>/dev/null', shell=True, stdout=subprocess.PIPE)
    p.wait()
    p = subprocess.Popen(SZYBKI+' -entropy AN  -prefix '+ folder+'/unb_'+filetag+' -sheffield -in '+folder+'/'+filetag+'.oeb.gz -strict false 2>/dev/null', shell=True, stdout=subprocess.PIPE)
    p.wait()

    p_bound = subprocess.Popen('grep \'At T=298.15 K,\' '+folder+'/bou_'+filetag+'.log|awk \'{print $8}\'|tail -1', shell=True, stdout=subprocess.PIPE)
    bound = p_bound.communicate()[0].split()[0]

    p_unbound = subprocess.Popen('grep \'At T=298.15 K,\' '+folder+'/unb_'+filetag+'.log|awk \'{print $8}\'|tail -1', shell=True, stdout=subprocess.PIPE)
    unbound = p_unbound.communicate()[0].split()[0]
except:
    bound = 'FAIL'
    unbound = 'FAIL'

entropy = float(bound) - float(unbound)

    
#calculating ligand descriptors
try:
    p_chemax = subprocess.Popen(CXCALC+' -N i logp pienergy acceptorcount donorcount fsp3 polarsurfacearea vdwsa '+folder+'/mini_complex_'+filetag+'_lig.sdf |tail -1' , shell=True, stdout= subprocess.PIPE)
    chemax = p_chemax.communicate()[0].split()
    #chemax = ['0.0','0.0','0.0','0.0','0.0','0.0','0.0' ]
except:
    chemax = ['FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL']

logp= chemax[0]
pienergy= chemax[1]
acceptorcount= chemax[2]
donorcount= chemax[3]
fsp3= chemax[4]
polarsurfacearea= chemax[5]
vdwsa= chemax[6]

#calculating rfscore
try:
    p_rfscore = subprocess.Popen('python  '+VSCREENML+'/pdb_ligand_distance_rf_feature.py '+folder+'/mini_complex_'+filetag+'.pdb |tail -1 ', shell=True, stdout= subprocess.PIPE)
    rfscore = p_rfscore.communicate()[0].split()

except:
    rfscore = ['FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL']

C_C = rfscore[0]
N_C = rfscore[1]
O_C = rfscore[2]
S_C = rfscore[3]
C_N = rfscore[4]
N_N = rfscore[5]
O_N = rfscore[6]
S_N = rfscore[7]
C_O = rfscore[8]
N_O = rfscore[9]
O_O = rfscore[10]
S_O= rfscore[11]
C_F = rfscore[12]
N_F = rfscore[13]
O_F = rfscore[14]
S_F = rfscore[15]
C_P =rfscore[16]
N_P = rfscore[17]
O_P = rfscore[18]
S_P = rfscore[19]
C_S = rfscore[20]
N_S = rfscore[21]
O_S = rfscore[22]
S_S = rfscore[23]
C_Cl = rfscore[24]
N_Cl = rfscore[25]
O_Cl = rfscore[26]
S_Cl = rfscore[27] 
C_Br = rfscore[28]
N_Br = rfscore[29]
O_Br = rfscore[30]
S_Br = rfscore[31]
C_I = rfscore[32] 
N_I = rfscore[33]
O_I = rfscore[34]
S_I = rfscore[35]

#>>> print 'We are the {} who say "{}!"'.format('knights', 'Ni')
#We are the knights who say "Ni!"

print ('{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}'.format(filetag, fa_atr, fa_rep, fa_sol, fa_elec, hbond_bb_sc, hbond_sc, Interface_Energy, Total_BSA, Interface_HB, Total_packstats, Interface_unsat, Total_pose_exposed_SASA, interface_hydrophobic_sasa, interface_polar_sasa, C_C, N_C, O_C, S_C, C_N, N_N, O_N, S_N, C_O, N_O, O_O, S_O, C_F, N_F, O_F, S_F, C_P, N_P, O_P, S_P, C_S, N_S, O_S, S_S, C_Cl, N_Cl, O_Cl, S_Cl, C_Br, N_Br, O_Br, S_Br, C_I, N_I, O_I, S_I, SIDE_flex_ALPHA, SIDE_flex_BETA, SIDE_flex_OTHER, BACK_flex_ALPHA, BACK_flex_BETA, BACK_flex_OTHER, TotalElec, TotalHbond, Hphobe_contact, Pi_Pi, Tstack, Cationpi, Salt_Bridge, entropy, pienergy, fsp3, polarsurfacearea, vdwsa))

data = np.array([fa_atr, fa_rep, fa_sol, fa_elec, hbond_bb_sc, hbond_sc, Interface_Energy, Total_BSA, Interface_HB, Total_packstats, Interface_unsat, Total_pose_exposed_SASA, interface_hydrophobic_sasa, interface_polar_sasa, C_C, N_C, O_C, S_C, C_N, N_N, O_N, S_N, C_O, N_O, O_O, S_O, C_F, N_F, O_F, S_F, C_P, N_P, O_P, S_P, C_S, N_S, O_S, S_S, C_Cl, N_Cl, O_Cl, S_Cl, C_Br, N_Br, O_Br, S_Br, C_I, N_I, O_I, S_I, SIDE_flex_ALPHA, SIDE_flex_BETA, SIDE_flex_OTHER, BACK_flex_ALPHA, BACK_flex_BETA, BACK_flex_OTHER, TotalElec, TotalHbond, Hphobe_contact, Pi_Pi, Tstack, Cationpi, Salt_Bridge, entropy, pienergy, fsp3, polarsurfacearea, vdwsa]).reshape((1,-1)).astype(float)

os.remove(folder+'/'+filetag+'.params')
os.remove(folder+'/'+filetag+'_0001.pdb')
os.remove(folder+'/mini_complex_'+filetag+'_lig1.sdf')
os.remove(folder+'/mini_complex_'+filetag+'_binana.out')
os.remove(folder+'/mini_complex_'+filetag+'_unb.pdbqt')
os.remove(folder+'/mini_complex_'+filetag+'_lig.pdbqt')
os.remove(folder+'/mini_complex_'+filetag+'_unb.pdb')
os.remove(folder+'/mini_complex_'+filetag+'_lig.pdb')
os.remove(folder+'/'+filetag+'.oeb.gz')
os.remove(folder+'/bou_'+filetag+'.log')
os.remove(folder+'/bou_'+filetag+'_out.oeb')
os.remove(folder+'/bou_'+filetag+'.status')
os.remove(folder+'/bou_'+filetag+'.param')
os.remove(folder+'/unb_'+filetag+'.log')
os.remove(folder+'/unb_'+filetag+'_out.oeb')
os.remove(folder+'/unb_'+filetag+'.status')
os.remove(folder+'/unb_'+filetag+'.param')

#predict vscreenml scrore

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
predict(clf, filetag, df)
