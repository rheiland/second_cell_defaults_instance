/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

Cell_Definition elmar_cell_def;

void create_default_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 


    // ------- 2nd cell default: copy what's done in P*_standard_models.cpp: initialize_default_cell_definition()

	elmar_cell_def.pMicroenvironment = BioFVM::get_default_microenvironment();
	elmar_cell_def.phenotype.secretion.sync_to_current_microenvironment();
	
	elmar_cell_def.type = 1; 
	elmar_cell_def.name = "elmar cell def";

	elmar_cell_def.parameters.pReference_live_phenotype = &(elmar_cell_def.phenotype); 
	
	elmar_cell_def.functions.cycle_model = Ki67_advanced; 
	
	elmar_cell_def.functions.volume_update_function = standard_volume_update_function;
	elmar_cell_def.functions.update_migration_bias = NULL; 
	
	elmar_cell_def.functions.update_phenotype = update_cell_and_death_parameters_O2_based; // NULL; 
	elmar_cell_def.functions.custom_cell_rule = NULL; 
	
	elmar_cell_def.functions.update_velocity = standard_update_cell_velocity;
	elmar_cell_def.functions.add_cell_basement_membrane_interactions = NULL; 
	elmar_cell_def.functions.calculate_distance_to_membrane = NULL; 
	
	elmar_cell_def.functions.set_orientation = NULL;
	
	elmar_cell_def.functions.plot_agent_SVG = standard_agent_SVG;
	elmar_cell_def.functions.plot_agent_legend = standard_agent_legend;
	
	elmar_cell_def.phenotype.death.add_death_model( 0.00319/60.0 , &apoptosis , apoptosis_parameters );
	elmar_cell_def.phenotype.death.add_death_model( 0.0 , &necrosis , necrosis_parameters );
	
	// set up the default phenotype (to be consistent with the default functions)
	elmar_cell_def.phenotype.cycle.sync_to_cycle_model( elmar_cell_def.functions.cycle_model ); 
	
	elmar_cell_def.phenotype.cell_interactions.sync_to_cell_definitions(); 
	elmar_cell_def.phenotype.cell_transformations.sync_to_cell_definitions(); 
	elmar_cell_def.phenotype.motility.sync_to_current_microenvironment(); 
	elmar_cell_def.phenotype.mechanics.sync_to_cell_definitions(); 


    //--------------------------------------------
	build_cell_definitions_maps(); 

	int m = microenvironment.number_of_densities(); 
	int n = cell_definition_indices_by_name.size(); 
    std::cout << __FUNCTION__ << ": # of densities= " << m << ", # of cell defs= " << n << std::endl;

	setup_signal_behavior_dictionaries(); 	

	display_cell_definitions( std::cout ); 
	
	return; 
}

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 


    // Cell_Definition MCF7 = cell_defaults;
    // MCF7 = cell_defaults;

    // ------------ Create another Cell_Definition 
    // cell_defaults_two.type = 1; 
    // cell_defaults_two.name = "elmar"; 
    // // cell_defaults_two.parameters.pReference_live_phenotype = &( cell_defaults_two.phenotype ); 
    // // cell_defaults_two.parameters.pReference_live_phenotype = &( cell_defaults.phenotype ); 
    // cell_defaults_two.phenotype.cycle.pCycle_Model = ( cell_defaults.phenotype.cycle.pCycle_Model ); 
    // cell_defaults_two.phenotype.death = ( cell_defaults.phenotype.death ); 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

    // ------------ Create another Cell_Definition 
    // initialize_another_cell_def( cell_defaults_two );

    // cell_defaults_two = cell_defaults;
    // cell_defaults_two.type = 1; 
    // cell_defaults_two.name = "elmar"; 
    // // cell_defaults_two.parameters.pReference_live_phenotype = &( cell_defaults_two.phenotype ); 
    // cell_defaults_two.parameters.pReference_live_phenotype = &( cell_defaults.phenotype ); 
    // // cell_defaults_two.phenotype.cycle.pCycle_Model = &( cell_defaults_two.phenotype.cycle.pCycle_Model ); 
    // // cell_defaults_two.phenotype.cycle.pCycle_Model = cell_defaults_two.phenotype.cycle.pCycle_Model; 
    // // cell_defaults_two.parameters.pReference_live_phenotype = &( cell_defaults.phenotype ); 
    // // cell_defaults_two.parameters = Cell_Parameters();

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 

	Cell_Definition* pCD = find_cell_definition( "cancer cell"); 
	pCD->functions.update_phenotype = tumor_cell_phenotype_with_oncoprotein; 

	pCD->parameters.o2_proliferation_saturation = 38; 
	pCD->parameters.o2_reference = 38; 


    // rwh
	// Cell_Definition* pCD = find_cell_definition( "cancer cell"); 
    // Cell_Definition cell_defaults_three;

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 

	// custom placement 

	Cell_Definition* pCD = find_cell_definition( "cancer cell"); 
	double cell_radius = pCD->phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double tumor_radius = parameters.doubles( "tumor_radius" ); // 250.0; 
	
	// Parameter<double> temp; 
	
	int i = parameters.doubles.find_index( "tumor_radius" ); 
	
	Cell* pCell = NULL; 
	
	double x = 0.0; 
	double x_outer = tumor_radius; 
	double y = 0.0; 
	
	double p_mean = parameters.doubles( "oncoprotein_mean" ); 
	double p_sd = parameters.doubles( "oncoprotein_sd" ); 
	double p_min = parameters.doubles( "oncoprotein_min" ); 
	double p_max = parameters.doubles( "oncoprotein_max" ); 
	
	int n = 0; 
	while( y < tumor_radius )
	{
		x = 0.0; 
		if( n % 2 == 1 )
		{ x = 0.5*cell_spacing; }
		x_outer = sqrt( tumor_radius*tumor_radius - y*y ); 
		
		while( x < x_outer )
		{
			pCell = create_cell( *pCD ); // tumor cell 
			pCell->assign_position( x , y , 0.0 );
			double p = NormalRandom( p_mean, p_sd );
			if( p < p_min )
			{ p = p_min; }
			if( p > p_max )
			{ p = p_max; }
			set_single_behavior( pCell, "custom:oncoprotein" , p ); 
			
			if( fabs( y ) > 0.01 )
			{
				pCell = create_cell(*pCD); // tumor cell 
				pCell->assign_position( x , -y , 0.0 );
				double p = NormalRandom( p_mean, p_sd );
				if( p < p_min )
				{ p = p_min; }
				if( p > p_max )
				{ p = p_max; }
				set_single_behavior( pCell, "custom:oncoprotein" , p ); 
			}
			
			if( fabs( x ) > 0.01 )
			{ 
				pCell = create_cell(*pCD); // tumor cell 
				pCell->assign_position( -x , y , 0.0 );
				double p = NormalRandom( p_mean, p_sd );
				if( p < p_min )
				{ p = p_min; }
				if( p > p_max )
				{ p = p_max; }
				set_single_behavior( pCell, "custom:oncoprotein" , p ); 
		
				if( fabs( y ) > 0.01 )
				{
					pCell = create_cell(*pCD); // tumor cell 
					pCell->assign_position( -x , -y , 0.0 );
					double p = NormalRandom( p_mean, p_sd );
					if( p < p_min )
					{ p = p_min; }
					if( p > p_max )
					{ p = p_max; }
					set_single_behavior( pCell, "custom:oncoprotein" , p ); 

				}
			}
			x += cell_spacing; 
			
		}
		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
	
	double sum = 0.0; 
	double min = 9e9; 
	double max = -9e9; 
	for( int i=0; i < all_cells->size() ; i++ )
	{
		double r = get_single_signal( (*all_cells)[i] , "custom:oncoprotein" ); 
		sum += r;
		if( r < min )
		{ min = r; } 
		if( r > max )
		{ max = r; }
	}
	double mean = sum / ( all_cells->size() + 1e-15 ); 
	// compute standard deviation 
	sum = 0.0; 
	for( int i=0; i < all_cells->size(); i++ )
	{
		double r = get_single_signal( (*all_cells)[i] , "custom:oncoprotein" ); 
		sum +=  ( r - mean )*( r - mean ); 
	}
	double standard_deviation = sqrt( sum / ( all_cells->size() - 1.0 + 1e-15 ) ); 
	
	std::cout << std::endl << "Oncoprotein summary: " << std::endl
			  << "===================" << std::endl; 
	std::cout << "mean: " << mean << std::endl; 
	std::cout << "standard deviation: " << standard_deviation << std::endl; 
	std::cout << "[min max]: [" << min << " " << max << "]" << std::endl << std::endl; 	
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

std::vector<std::string> heterogeneity_coloring_function( Cell* pCell )
{
	double p = get_single_signal( pCell, "custom:oncoprotein"); 
	
	static double p_min = parameters.doubles( "oncoprotein_min" ); 
	static double p_max = parameters.doubles( "oncoprotein_max" ); 
	
	// immune are black
	std::vector< std::string > output( 4, "black" ); 
	
	if( pCell->type == 1 )
	{ return output; } 
	
	// live cells are green, but shaded by oncoprotein value 
	if( pCell->phenotype.death.dead == false )
	{
		int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (p-p_min) * 255.0 ); 
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein );
		output[0].assign( szTempString );
		output[1].assign( szTempString );

		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/p_max) , (int)round(output[0][1]/p_max) , (int)round(output[0][2]/p_max) );
		output[2].assign( szTempString );
		
		return output; 
	}

	// if not, dead colors 
	
	if( get_single_signal( pCell, "apoptotic") > 0.5 )
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown
	if( get_single_signal(pCell, "necrotic") > 0.5 )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}	
	
	return output; 
}

void tumor_cell_phenotype_with_oncoprotein( Cell* pCell, Phenotype& phenotype, double dt )
{
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// if cell is dead, don't bother with future phenotype changes. 
	if( get_single_signal( pCell, "dead") > 0.5 )
	{
		pCell->functions.update_phenotype = NULL; 		
		return; 
	}

	// multiply proliferation rate by the oncoprotein 

	double cycle_rate = get_single_behavior( pCell, "cycle entry"); 
	cycle_rate *= get_single_signal( pCell , "custom:oncoprotein"); 
	set_single_behavior( pCell, "cycle entry" , cycle_rate ); 
	
	return; 
}