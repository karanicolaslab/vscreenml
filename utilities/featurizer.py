#! usr/bin/python
import os, subprocess, sys
from configuration import config as cfg

__author__ = "Yusuf Adeshina"
__email__ = "dryusufadeshina@gmail.com"

class VscreenMLFeaturizer(object):

    def __init__(self, folder, filetag):
	#initializing the class
        self.filetag = filetag
        self.folder = folder
        self.feature_list = []
        self.purpose = 'This featurizes protein-ligand complexes for vscreenml scoring'
        
    def prepare_protein_and_ligand_for_scoring(self):
        #Extract ligands and protein structures
        p=subprocess.Popen('grep -E \"ATOM|ZN|MG|MN\" '+self.folder+'/mini_complex_'+self.filetag+'.pdb > ' +self.folder+'/mini_complex_'+self.filetag+'_unb.pdb && grep HETATM '+self.folder+'/mini_complex_'+self.filetag+'.pdb|grep -v \'ZN\'|grep -v \'MG\'|grep -v \'MN\' > '+self.folder+'/mini_complex_'+self.filetag+'_lig.pdb',shell=True)

        p.wait()

        #reformatting the ligands
        p=subprocess.Popen(cfg.BABEL+'/bin/obabel '+self.folder+'/mini_complex_'+self.filetag+'_lig.pdb -O '+self.folder+'/mini_complex_'+self.filetag+'_lig1.sdf >/dev/null 2>&1',shell=True)
        p.wait()

        p=subprocess.Popen(cfg.CHEMAXON+'/structure-representation-toolkit/bin/structurecheck -c \"valenceerror\" -m fix -f sdf -o '+self.folder+'/mini_complex_'+self.filetag+'_lig.sdf '+self.folder+'/mini_complex_'+self.filetag+'_lig1.sdf',shell=True)
        p.wait()

    def generate_params_file(self):
        #generating params filetag
        p=subprocess.Popen(cfg.ROSETTA_DIR+'/main/source/scripts/python/public/molfile_to_params.py '+self.folder+'/mini_complex_'+self.filetag+'_lig1.sdf -p '+self.folder+'/'+self.filetag+' >/dev/null 2>&1' , shell=True)
        p.wait()



    def rosetta_minimization(self):
       #minimization of the complex
        try:
            p_rosetta = subprocess.Popen(cfg.ROSETTA_DIR+'/main/source/bin/minimize_ppi.'+cfg.ROSETTA_BUILD+' -s '+self.folder+'/mini_complex_'+self.filetag+'.pdb -extra_res_fa '+self.folder+'/'+self.filetag+'.params -database '+cfg.ROSETTA_DIR+'/main/database -ignore_zero_occupancy false -score_only |grep \'Interface_Scores:mini\' |awk \'{print $4,$5,$6,$7,$8}\' ', shell=True, stdout= subprocess.PIPE)
            rosetta = p_rosetta.communicate()[0].split()
	
        except:
            #use typical value (median) from the training set to impute missing value
            rosetta= ['-23.6933', '719.7050', '1.0000', '0.6410', '4.0000']
            
       # print(rosetta) 
       #Interface_Energy = rosetta[0], Total_BSA = rosetta[1], Interface_HB = rosetta[2], Total_packstats = rosetta[3], #Interface_unsat =rosetta[4]
        return rosetta[0], rosetta[1], rosetta[2], rosetta[3], rosetta[4]

    def calculate_sasa(self):
        #calculating sasa   
        try:
            p_sasa = subprocess.Popen(cfg.ROSETTA_DIR+'/main/source/bin/ragul_get_ligand_sasa.'+cfg.ROSETTA_BUILD+' -database '+cfg.ROSETTA_DIR+'/main/database -input_bound_pdb '+self.folder+'/mini_complex_'+self.filetag+'.pdb -input_unbound_pdb '+self.folder+'/mini_complex_'+self.filetag+'_unb.pdb -input_ligand_pdb '+self.folder+'/mini_complex_'+self.filetag+'_lig.pdb -extra_res_fa '+self.folder+'/'+self.filetag+'.params |grep \'Scores:\' |awk \'{print $2,$3,$4}\'', shell=True, stdout= subprocess.PIPE)
            sasa = p_sasa.communicate()[0].split()

        except:
            #use typical value (median) from the training set to impute missing value
            sasa = ['0.3304', '0.7138', '0.2867']

        #Total_pose_exposed_SASA = sasa[0], interface_hydrophobic_sasa = sasa[1], interface_polar_sasa = sasa[2]
        return sasa[0], sasa[1], sasa[2]


    def calculate_component_score(self):
        #calculating componenet scores    
        try:
            p_compo = subprocess.Popen(cfg.ROSETTA_DIR+'/main/source/bin/yusuf_interface_ddg.'+cfg.ROSETTA_BUILD+' -database '+cfg.ROSETTA_DIR+'/main/database -s '+self.folder+'/mini_complex_'+self.filetag+'.pdb -extra_res_fa '+self.folder+'/'+self.filetag+'.params -ignore_zero_occupancy false   |grep \'Component_score:\' |awk \'{print $4,$5,$6,$10,$11,$12}\'', shell=True, stdout= subprocess.PIPE)
            compo = p_compo.communicate()[0].split()

        except:
            #use typical value (median) from the training set to impute missing value
            compo = ['-28.0043', '6.45592', '7.77609', '-4.8185', '0.0000','0.0000']

        #fa_atr = compo[0], fa_rep = compo[1], fa_sol = compo[2], fa_elec = compo[3], hbond_bb_sc = compo[4], hbond_sc= compo[5]
        return compo[0], compo[1], compo[2], compo[3], compo[4], compo[5]

    def generate_dot_pdbqt_files(self):
        #generating .pdbqt filetags for cfg.BINANA features
        p=subprocess.Popen(cfg.MGLTOOLS+'/bin/pythonsh '+cfg.MGLTOOLS+'/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r '+self.folder+'/mini_complex_'+self.filetag+'_unb.pdb -A hydrogens -o '+self.folder+'/mini_complex_'+self.filetag+'_unb.pdbqt >/dev/null 2>&1',shell=True)
        p.wait()

        p=subprocess.Popen(cfg.MGLTOOLS+'/bin/pythonsh '+cfg.MGLTOOLS+'/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l '+self.folder+'/mini_complex_'+self.filetag+'_lig.pdb -A hydrogens -o '+self.folder+'/mini_complex_'+self.filetag+'_lig.pdbqt >/dev/null 2>&1',shell=True)
        p.wait()

    def generate_binana_features(self):
        #calculating cfg.BINANA features
        try:
            p=subprocess.Popen('python '+cfg.BINANA+'/binana_1_2_0.py -receptor '+self.folder+'/mini_complex_'+self.filetag+'_unb.pdbqt -ligand '+self.folder+'/mini_complex_'+self.filetag+'_lig.pdbqt > '+self.folder+'/mini_complex_'+self.filetag+'_binana.out',shell=True)
            p.wait()
            p_binana = subprocess.Popen('python '+cfg.BINANA+'/parse_binana.py '+self.folder+'/mini_complex_'+self.filetag+'_binana.out', shell=True, stdout=subprocess.PIPE)
            binana = p_binana.communicate()[0].split()
        except:
            #use typical value (median) from the training set to impute missing value
            binana = ['5', '22', '17', '0', '4', '5', '-56280.936', '1', '21', '0', '0', '0', '0']

        #SIDE_flex_ALPHA = binana[0], SIDE_flex_BETA= binana[1], SIDE_flex_OTHER=binana[2], BACK_flex_ALPHA= binana[3], BACK_flex_BETA= binana[4], BACK_flex_OTHER= binana[5], TotalElec= binana[6], TotalHbond= binana[7], Hphobe_contact= binana[8], Pi_Pi= binana[9], Tstack= binana[10], Cationpi= binana[11], Salt_Bridge= binana[12]
        return binana[0], binana[1], binana[2], binana[3], binana[4], binana[5], binana[6], binana[7], binana[8], binana[9], binana[10], binana[11], binana[12]

    def calculate_szybki(self):
        #generating cfg.OMEGA conformers for entropy calculations
        p=subprocess.Popen(cfg.OMEGA+' -in '+self.folder+'/mini_complex_'+self.filetag+'_lig1.sdf -out '+self.folder+'/'+self.filetag+'.oeb.gz -strictstereo false -strictatomtyping false -warts true -maxconfs 1000 -rms 0.1 >/dev/null 2>&1',shell=True)
        p.wait()

        try:
            p = subprocess.Popen(cfg.SZYBKI+' -entropy AN  -prefix '+self.folder+'/bou_'+self.filetag+' -complex '+self.folder+'/mini_complex_'+self.filetag+'.pdb -strict false >/dev/null 2>&1', shell=True, stdout=subprocess.PIPE)
            p.wait()
            p = subprocess.Popen(cfg.SZYBKI+' -entropy AN  -prefix '+self.folder+'/unb_'+self.filetag+' -sheffield -in '+self.folder+'/'+self.filetag+'.oeb.gz -strict false >/dev/null 2>&1', shell=True, stdout=subprocess.PIPE)
            p.wait()

            p_bound = subprocess.Popen('grep \'At T=298.15 K,\' '+self.folder+'/bou_'+self.filetag+'.log|awk \'{print $8}\'|tail -1', shell=True, stdout=subprocess.PIPE)
            bound = p_bound.communicate()[0].split()[0]

            p_unbound = subprocess.Popen('grep \'At T=298.15 K,\' '+self.folder+'/unb_'+self.filetag+'.log|awk \'{print $8}\'|tail -1', shell=True, stdout=subprocess.PIPE)
            unbound = p_unbound.communicate()[0].split()[0]
            entropy = float(bound) - float(unbound)
        except:
            #use typical value (median) from the training set to impute missing value
            entropy = '-1.4'
        return entropy

   
    def calculate_chemaxon_features(self): 
        #calculating ligand descriptors
   
        try:
            p_chemax = subprocess.Popen(cfg.CXCALC+' -N i logp pienergy acceptorcount donorcount fsp3 polarsurfacearea vdwsa '+self.folder+'/mini_complex_'+self.filetag+'_lig.sdf |tail -1' , shell=True, stdout= subprocess.PIPE)
            chemax = p_chemax.communicate()[0].split()
          
        except:
            #use typical value (median) from the training set to impute missing value
            chemax = ['FAIL', '39.46', 'FAIL', 'FAIL', '0.28', '79.90', '468.78']

        #logp= chemax[0], pienergy= chemax[1], acceptorcount= chemax[2], donorcount= chemax[3], fsp3= chemax[4], polarsurfacearea= chemax[5], vdwsa = chemax[6]
        return chemax[1], chemax[4], chemax[5], chemax[6]

    def calculate_rfscore_like_features(self):
        #calculating rfscore
        try:
            p_rfscore = subprocess.Popen('python  '+cfg.VSCREENML+'/rfscore_like_feature.py '+self.folder+'/mini_complex_'+self.filetag+'.pdb |tail -1 ', shell=True, stdout= subprocess.PIPE)
            rfscore = p_rfscore.communicate()[0].split()

        except:
            #use typical value (median) from the training set to impute missing value
            rfscore = ['2884', '743', '805', '27', '490', '128', '139', '4', '390', '99', '108', '3', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0']

        #order of features '6.6','7.6','8.6','16.6','6.7','7.7','8.7','16.7','6.8','7.8','8.8','16.8','6.9','7.9','8.9','16.9','6.15','7.15','8.15','16.15','6.16','7.16','8.16','16.16','6.17','7.17','8.17','16.17','6.35','7.35','8.35','16.35','6.53','7.53','8.53','16.53'
        return rfscore[0], rfscore[1], rfscore[2], rfscore[3], rfscore[4], rfscore[5], rfscore[6], rfscore[7], rfscore[8], rfscore[9], rfscore[10], rfscore[11], rfscore[12], rfscore[13], rfscore[14], rfscore[15], rfscore[16], rfscore[17], rfscore[18], rfscore[19], rfscore[20], rfscore[21], rfscore[22], rfscore[23], rfscore[24], rfscore[25], rfscore[26], rfscore[27], rfscore[28], rfscore[29], rfscore[30], rfscore[31], rfscore[32], rfscore[33], rfscore[34], rfscore[35]

    def clean_up(self):
        os.remove(self.folder+'/'+self.filetag+'.params')
        os.remove(self.folder+'/'+self.filetag+'_0001.pdb')
        os.remove(self.folder+'/mini_complex_'+self.filetag+'_lig1.sdf')
        os.remove(self.folder+'/mini_complex_'+self.filetag+'_binana.out')
        os.remove(self.folder+'/mini_complex_'+self.filetag+'_unb.pdbqt')
        os.remove(self.folder+'/mini_complex_'+self.filetag+'_lig.pdbqt')
        os.remove(self.folder+'/mini_complex_'+self.filetag+'_unb.pdb')
        os.remove(self.folder+'/mini_complex_'+self.filetag+'_lig.pdb')
        os.remove(self.folder+'/'+self.filetag+'.oeb.gz')
        os.remove(self.folder+'/bou_'+self.filetag+'.log')
        os.remove(self.folder+'/bou_'+self.filetag+'_out.oeb')
        os.remove(self.folder+'/bou_'+self.filetag+'.status')
        os.remove(self.folder+'/bou_'+self.filetag+'.param')
        os.remove(self.folder+'/unb_'+self.filetag+'.log')
        os.remove(self.folder+'/unb_'+self.filetag+'_out.oeb')
        os.remove(self.folder+'/unb_'+self.filetag+'.status')
        os.remove(self.folder+'/unb_'+self.filetag+'.param')



class VscreenMLAggregateFeatures(object):
    def __init__(self, folder_, filetag_):
        #initializing the class
        self.filetag_ = filetag_
        self.folder_ = folder_
        self.feature_list_ = []
        self.purpose_ = 'This featurizes protein-ligand complexes for vscreenml scoring'

    def vscreenml_featurize(self):
        """
        This function computes the overall vscreenml features
        """
        features_ = VscreenMLFeaturizer(self.folder_, self.filetag_)

	features_.prepare_protein_and_ligand_for_scoring()
        features_.generate_params_file()
        self.feature_list_.extend(features_.calculate_component_score())
        self.feature_list_.extend(features_.rosetta_minimization())
        self.feature_list_.extend(features_.calculate_sasa())
        self.feature_list_.extend(features_.calculate_rfscore_like_features())
        features_.generate_dot_pdbqt_files()
        self.feature_list_.extend(features_.generate_binana_features())
        self.feature_list_.append(features_.calculate_szybki())
        self.feature_list_.extend(features_.calculate_chemaxon_features())
        features_.clean_up()
        return self.feature_list_

#feature order['fa_atr','fa_rep','fa_sol','fa_elec','hbond_bb_sc','hbond_sc','Interface_Energy','Total_BSA','Interface_HB','Total_packstats','Interface_unsat','Total_pose_exposed_SASA','interface_hydrophobic_sasa','interface_polar_sasa','6.6','7.6','8.6','16.6','6.7','7.7','8.7','16.7','6.8','7.8','8.8','16.8','6.9','7.9','8.9','16.9','6.15','7.15','8.15','16.15','6.16','7.16','8.16','16.16','6.17','7.17','8.17','16.17','6.35','7.35','8.35','16.35','6.53','7.53','8.53','16.53','SIDE_flex_ALPHA','SIDE_flex_BETA','SIDE_flex_OTHER','BACK_flex_ALPHA','BACK_flex_BETA','BACK_flex_OTHER','TotalElec','TotalHbond','Hphobe_contact','Pi_Pi','T-stack','Cation-pi','Salt_Bridge','diff_entropy','pienergy','fsp3','polarsurfacearea','vdwsa']

