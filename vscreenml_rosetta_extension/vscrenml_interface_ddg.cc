// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief Prints the component energies scores without their weight of a pose
/// @author Yusuf Adeshina
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <string>
#include <core/io/pdb/pdb_writer.hh>




using namespace basic::options::OptionKeys;
using basic::options::option;
using core::Size;


static basic::Tracer TR( "apps.pilot.yusuf_interface_energy" );



////////////////////////////////////////////////////////////////////////////////
//                                    MAIN                                    //
////////////////////////////////////////////////////////////////////////////////
/*
int main( int argc, char * argv [] )
{
	devel::init(argc, argv);

	// load pose from pdb file
	core::pose::Pose ps;
	std::string const input_pdb_name( basic::options::start_file() );
	core::import_pose::pose_from_file( ps, input_pdb_name, core::import_pose::PDB_file ); */
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*	utility::vector1< std::pair< std::string,core::Real> > get_component_energy(core::pose::Pose const & protein_pose){
	
	core::scoring::ScoreFunctionOP scorefxn(
		core::scoring::get_score_function());

	(*scorefxn)(protein_pose);

	core::scoring::EnergyMap emap = protein_pose.energies().total_energies();
	core::scoring::EnergyMap weights = protein_pose.energies().weights();
	core::Real w_fa_sol = emap[core::scoring::fa_sol]*weights[core::scoring::fa_sol];

	std::string energy_terms[] = {"ref","fa_atr","fa_rep","fa_sol","fa_intra_rep","fa_elec","fa_pair","pro_close","hbond_bb_sc","hbond_sc","dslf_fa13","rama","omega","fa_dun","p_aa_pp","occ_sol_fitted","fa_plane"};
	utility::vector1< std::pair< std::string,core::Real> > component_vector[17];
	for (std::string::const_iterator i_e = energy_terms.begin(); i_e != energy_terms.end(); ++i_e) {
	core::Real "e_"+i_e = emap[core::scoring::i_e];
	std::pair< std::string,core::Real> energy_term_pair(i_e,"e_"+i_e);
	component_vector.push_back(energy_term_pair);
	}
	TR << component_vector << std::endl;
	return component_vector;
}

	utility::vector1< std::pair< std::string,core::Real> > component_vector;

	core::Real e_fa_atr = emap[core::scoring::fa_atr];
	std::pair< std::string,core::Real> energy_term_pair1("fa_atr",e_fa_atr);
	component_vector.push_back(energy_term_pair1);


	core::Real e_fa_rep = emap[core::scoring::fa_rep];
	std::pair< std::string,core::Real> energy_term_pair2("fa_rep",e_fa_rep);
	 component_vector.push_back(energy_term_pair2);


	core::Real e_fa_sol = emap[core::scoring::fa_sol];
	std::pair< std::string,core::Real> energy_term_pair3("fa_sol",e_fa_sol);
	component_vector.push_back(energy_term_pair3);


	core::Real e_fa_intra_rep = emap[core::scoring::fa_intra_rep];
	std::pair< std::string,core::Real> energy_term_pair4("fa_intra_rep",e_fa_intra);
   	component_vector.push_back(energy_term_pair4);


	core::Real e_fa_elec = emap[core::scoring::fa_elec];
	std::pair< std::string,core::Real> energy_term_pair5("fa_elec",e_fa_elec);
	component_vector.push_back(energy_term_pair5);


	core::Real e_fa_pair = emap[core::scoring::fa_pair];
	std::pair< std::string,core::Real> energy_term_pair6("fa_pair",e_fa_pair);
	component_vector.push_back(energy_term_pair6);

	
        core::Real e_pro_close = emap[core::scoring::pro_close];
	std::pair< std::string,core::Real> energy_term_pair7("pro_close",e_pro_close);
	component_vector.push_back(energy_term_pair7);


        core::Real e_hbond_bb_sc = emap[core::scoring::hbond_bb_sc];
	std::pair< std::string,core::Real> energy_term_pair8("hbond_bb_sc",e_hbond_bb_sc);
	component_vector.push_back(energy_term_pair8);


        core::Real e_hbond_sc = emap[core::scoring::hbond_sc];
	std::pair< std::string,core::Real> energy_term_pair9("hbond_sc",e_hbond_sc);
	 component_vector.push_back(energy_term_pair9);


        core::Real e_dslf_fa13 = emap[core::scoring::dslf_fa13];
	std::pair< std::string,core::Real> energy_term_pair10("dslf_fa13",e_dslf_fa13);
	 component_vector.push_back(energy_term_pair10);


        core::Real e_rama = emap[core::scoring::rama];
	std::pair< std::string,core::Real> energy_term_pair11("rama",e_rama);
	 component_vector.push_back(energy_term_pair11);

	core::Real e_omega = emap[core::scoring::omega];
	std::pair< std::string,core::Real> energy_term_pair12("omega",e_omega);
	 component_vector.push_back(energy_term_pair12);


        core::Real e_fa_dun = emap[core::scoring::fa_dun];
	std::pair< std::string,core::Real> energy_term_pair13("fa_dun",e_fa_dun);
	 component_vector.push_back(energy_term_pair13);


        core::Real e_p_aa_pp = emap[core::scoring::p_aa_pp];
	std::pair< std::string,core::Real> energy_term_pair14("p_aa_pp",e_p_aa_pp);
	 component_vector.push_back(energy_term_pair14;

        core::Real e_fa_occ_sol_fitted = emap[core::scoring::occ_sol_fitted];
	std::pair< std::string,core::Real> energy_term_pair15("occ_sol_fitted",e_occ_sol_fitted);
	 component_vector.push_back(energy_term_pair15);

        core::Real e_fa_plane = emap[core::scoring::fa_plane];
	std::pair< std::string,core::Real> energy_term_pair16("fa_plane",e_fa_plane);
	 component_vector.push_back(energy_term_pair16);
        
	core::Real e_ref = emap[core::scoring::ref];
        std::pair< std::string,core::Real> energy_term_pair17("ref",e_ref);
        component_vector.push_back(energy_term_pair17);


	// TR << "fa_sol: " << w_fa_sol << " " << "lk_polar+lk_nonpolar " << (w_lk_polar+w_lk_nonpolar) << std::endl;
	TR << "fa_sol: " << e_fa_sol << " " << "lk_nonpolar " << e_fa_rep << std::endl;

	return component_vector;
}

*/

//using namespace basic::options::OptionKeys;
//using basic::options::option;

OPT_KEY( Boolean, min_first )

//static basic::Tracer TR( "apps.pilot.bou-min-ubo-nrg-jump" );


////////////////////////////////////////////////////////////////////////////////
//                                    MAIN                                    //
////////////////////////////////////////////////////////////////////////////////

int main( int argc, char * argv [] )
{
try {
	NEW_OPT( min_first, "First minimize the protein+ligand complex", false);

	devel::init(argc, argv);

	// load bound pose from pdb file
	core::pose::Pose bou_ps;
	std::string const input_pdb_name( basic::options::start_file() );
	core::import_pose::pose_from_file( bou_ps, input_pdb_name , core::import_pose::PDB_file);

	//create tag for output filename
	int pfounddir = input_pdb_name.find_last_of("/\\");
	int pfounddot = input_pdb_name.find_last_of(".");
	std::string tag = input_pdb_name.substr((pfounddir+1),(pfounddot-(pfounddir+1)));

	// create score function
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

	// minimize bound pose, if requested
	if(option[min_first]) {

		core::optimization::MinimizerOptions minoptions("lbfgs_armijo_nonmonotone", 0.00001, true);
		minoptions.nblist_auto_update( true );

		core::kinematics::MoveMap mm;
		mm.set_bb(true);
		mm.set_chi(true);
		mm.set_jump(true);

		core::optimization::AtomTreeMinimizer minimizer;
		minimizer.run(bou_ps, mm, *scorefxn, minoptions);

		bou_ps.dump_pdb("min_bou.pdb");
	}

	// score bound pose
	core::Real bou_score = (*scorefxn)(bou_ps);
	TR << "bound pose's energy: " << bou_score << std::endl;
	core::scoring::EnergyMap emap = bou_ps.energies().total_energies();

	// create unbound pose
        core::pose::Pose ubo_ps = bou_ps;
        core::Size const rb_jump = ubo_ps.num_jump();
        core::Real const ubo_dist = 1000000.0;
        protocols::rigid::RigidBodyTransMover trans_mover( ubo_ps, rb_jump );
        trans_mover.step_size(ubo_dist);
        trans_mover.apply(ubo_ps);

        // score unboud pose
        core::Real ubo_score = (*scorefxn)(ubo_ps);
        TR << "unbound pose's energy: " << ubo_score << std::endl;
	TR << "energy difference: " << (bou_score - ubo_score) << std::endl;
	core::scoring::EnergyMap emap1 = ubo_ps.energies().total_energies();
	//ubo_ps.dump_pdb( "unbound.pdb" );

//	utility::vector1< std::pair< std::string,core::Real> > component_vector;

        core::Real e_fa_atr = emap[core::scoring::fa_atr] - emap1[core::scoring::fa_atr];
        //std::pair< std::string,core::Real> energy_term_pair1("fa_atr",e_fa_atr);
        //component_vector.push_back(energy_term_pair1);


        core::Real e_fa_rep = emap[core::scoring::fa_rep] - emap1[core::scoring::fa_rep];
        //std::pair< std::string,core::Real> energy_term_pair2("fa_rep",e_fa_rep);
        //component_vector.push_back(energy_term_pair2);
	
	core::Real e_lk_polar = emap[core::scoring::lk_polar] - emap1[core::scoring::lk_polar];

        core::Real e_lk_nonpolar = emap[core::scoring::lk_nonpolar] - emap1[core::scoring::lk_nonpolar];

      	core::Real e_fa_sol = emap[core::scoring::fa_sol] - emap1[core::scoring::fa_sol];
        //std::pair< std::string,core::Real> energy_term_pair3("fa_sol",e_fa_sol);
       // component_vector.push_back(energy_term_pair3);


        core::Real e_fa_elec = emap[core::scoring::fa_elec] - emap1[core::scoring::fa_elec];
        //std::pair< std::string,core::Real> energy_term_pair5("fa_elec",e_fa_elec);
        //component_vector.push_back(energy_term_pair5);


        core::Real e_hbond_bb_sc = emap[core::scoring::hbond_bb_sc] - emap1[core::scoring::hbond_bb_sc];
        //std::pair< std::string,core::Real> energy_term_pair8("hbond_bb_sc",e_hbond_bb_sc);
        //component_vector.push_back(energy_term_pair8);


        core::Real e_hbond_sc = emap[core::scoring::hbond_sc] - emap1[core::scoring::hbond_sc];
        //std::pair< std::string,core::Real> energy_term_pair9("hbond_sc",e_hbond_sc);
        // component_vector.push_back(energy_term_pair9);

        core::Real e_occ_sol_fitted = emap[core::scoring::occ_sol_fitted] - emap1[core::scoring::occ_sol_fitted];
       // std::pair< std::string,core::Real> energy_term_pair15("occ_sol_fitted",e_occ_sol_fitted);
         //component_vector.push_back(energy_term_pair15);

	TR << "|||||||||||||||||||||||COMPONENT SCORES||||||||||||||||||||||||||||||||"<< std::endl;
	//TR << "fa_atr"<<"\t"<<"fa_rep"<<"\t"<<"fa_sol"<<"\t"<<"fa_elec"<<"\t"<<"hbond_bb_sc"<<"\t"<<"hbond_sc"<< std::endl;
	//TR <<"Component_score"<<"\t"<< e_fa_atr<<"\t"<<e_fa_rep<<"\t"<<e_fa_sol<<"\t"<<e_fa_elec<<"\t"<<e_hbond_bb_sc<<"\t"<<e_hbond_sc<< std::endl;

	TR << "tag"<<"\t"<<"fa_atr"<<"\t"<<"fa_rep"<<"\t"<<"fa_sol"<<"\t"<<"lk_polar"<<"\t"<<"lk_nonpolar"<<"\t"<<"occ_sol_fitted"<<"\t"<<"fa_elec"<<"\t"<<"hbond_bb_sc"<<"\t"<<"hbond_sc"<< std::endl;
        TR <<"Component_score:"<<"\t"<<tag<<"\t"<< e_fa_atr<<"\t"<<e_fa_rep<<"\t"<<e_fa_sol<<"\t"<<e_lk_polar<<"\t"<<e_lk_nonpolar<<"\t"<<e_occ_sol_fitted<<"\t"<<e_fa_elec<<"\t"<<e_hbond_bb_sc<<"\t"<<e_hbond_sc<< std::endl;

}
catch ( utility::excn::EXCN_Base const & e ) {
  std::cerr << "caught exception " << e.msg() << std::endl;
  return -1;
}
}
