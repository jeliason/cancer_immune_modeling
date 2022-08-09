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
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

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

	Cell_Definition* pCD = find_cell_definition( "cancer" ); 
	pCD->functions.update_phenotype = cancer_phenotype; 

	pCD = find_cell_definition( "macrophage" ); 
	pCD->functions.update_phenotype = macrophage_phenotype; 

	pCD = find_cell_definition( "CD8+ T cell" ); 
	pCD->functions.update_phenotype = CD8_Tcell_phenotype; 

	
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
	/*
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
	*/


	// place tumor cells 
	double max_distance = parameters.doubles("max_initial_distance");

	Cell_Definition* pCD = find_cell_definition( "cancer" ); 
	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
	for( int k=0 ; k < parameters.ints( "number_of_cancer_cells" ); k++ )
	{
		std::vector<double> position = {0,0,0}; 
		double r = sqrt(UniformRandom())* max_distance; 
		double theta = UniformRandom() * 6.2831853; 
		position[0] = r*cos(theta); 
		position[1] = r*sin(theta); 
		
		pC = create_cell( *pCD ); 
		pC->assign_position( position );
	}

	// place macrophages
	pCD = find_cell_definition( "macrophage" ); 

	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
	for( int k=0 ; k < parameters.ints( "number_of_macrophages" ); k++ )
	{
		std::vector<double> position = UniformOnUnitCircle(); 
		position *= parameters.doubles("macrophage_distance"); 
		
		pC = create_cell( *pCD ); 
		pC->assign_position( position );
	}

	// place CD8+ T cells
	pCD = find_cell_definition( "CD8+ T cell" ); 

	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
	for( int k=0 ; k < parameters.ints( "number_of_CD8_Tcells" ); k++ )
	{
		std::vector<double> position = UniformOnUnitCircle(); 
		position *= parameters.doubles("CD8_Tcell_distance"); 
		
		pC = create_cell( *pCD ); 
		pC->assign_position( position );
	}

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

void cancer_phenotype( Cell* pCell, Phenotype& p, double dt)
{
	// if dead, set secretion/uptake zero, release (export) debris 

	bool dead = (bool) get_single_signal( pCell, "dead"); 
	if( dead )
	{
		double volume = get_single_signal( pCell, "volume"); 

		set_single_behavior( pCell, "oxygen uptake" , 0.0 ); 
		set_single_behavior( pCell, "debris export" , 1*volume ); 

		return; 
	}

	// sample O2
	double o2 = get_single_signal( pCell, "oxygen"); 
	
	// set birth rate -- scale value from cell definition
	// set necrosis rate -- scale max value

	// set birth rate -- scale value from cell definition

		// base birth rate
	double rate0 = get_single_base_behavior( pCell , "cycle entry"); 
		// scale based on o2
	double o2_sat = get_single_signal( pCell , "custom:pO2_proliferation_saturation"); 
	double o2_threshold = 
		get_single_signal( pCell , "custom:pO2_proliferation_threshold"); 
	double rate = rate0 * linear_response_function( o2 , o2_threshold , o2_sat ); 
	set_single_behavior( pCell , "cycle entry" , rate ); 

	// set necrosis rate -- scale max value

		// get max rate
	double rateMax = get_single_behavior( pCell, "custom:max_necrosis_rate"); 
		// scale by O2 
	o2_sat = get_single_behavior( pCell, "custom:pO2_necrosis_saturation"); 
	o2_threshold = get_single_behavior( pCell, "custom:pO2_necrosis_threshold"); 
	rate = rateMax * decreasing_linear_response_function( o2, o2_sat , o2_threshold ); 
	set_single_behavior( pCell, "necrosis" , rate ); 

	// set apoptosis rate based on damage 

	double damage = get_single_signal( pCell , "damage"); 
	rate0 = get_single_base_behavior( pCell , "apoptosis"); 
	double halfmax = get_single_behavior( pCell , "custom:damage_halfmax"); // 180 
	double hillpower = get_single_behavior( pCell , "custom:damage_hillpower"); // 2 
	rateMax = get_single_behavior( pCell , "custom:damage_relative_max_death") * rate0; 
	double hill = Hill_response_function( damage , halfmax , hillpower ); 
	rate = rate0 + (rateMax-rate0)*hill; 
	set_single_behavior( pCell , "apoptosis" , rate) ; 

	return;
}

std::vector<std::string> custom_coloring_function( Cell* pCell )
{
	// start with color-by-number

	std::vector<std::string> output = paint_by_number_cell_coloring(pCell);

	// dead cells: black if apoptotic, brown if necrotic 

	bool dead = (bool) get_single_signal( pCell, "dead"); 
	if( dead )
	{
		if( pCell->phenotype.cycle.model().name == "Apoptosis" )
		{ output[0] = "rgb(0,0,0)"; }
		else
		{ output[0] = "rgb(111,78,55)"; }
	}

	// live tumor cells: shade by proliferation rate 

	if( pCell->type_name == "cancer" && dead == false ) 
	{
		// get relative birth rate
		double s = 10 * get_single_behavior( pCell, "cycle entry" ) 
			/ get_single_base_behavior( pCell, "cycle entry" ); 
		if( s > 1 )
		{ s = 1; }

		// make color
		int color = (int) round( 255.0 * s ); 
		char szColor [1024];
		// interpolate from blue to yellow 
		sprintf( szColor, "rgb(%u,%u,%u)",color,color,255-color );

		// modify output
		output[0] = szColor; 
		output[2] = szColor; 
		output[3] = szColor; 
	}

	// if live CD8+ T cell, color dark green 

	if( pCell->type_name == "CD8+ T cell" && dead == false ) 
	{
		// get relative birth rate
		char szColor [1024] = "rgb(0,128,0)"; 

		// modify output
		output[0] = szColor; 
		output[2] = szColor; 
		output[3] = szColor; 
	}

	// live macrophage cells: shade by anti-inflammatory state

	if( pCell->type_name == "macrophage" && dead == false ) 
	{
		// get relative secretion rate
		double s = 1 * get_single_behavior( pCell, "anti-inflammatory secretion" ) 
			/ get_single_base_behavior( pCell, "anti-inflammatory secretion" ); 
		if( s > 1 )
		{ s = 1; }

		// make color
		int color = (int) round( 255.0 * s ); 
		char szColor [1024];
		// interpolate from red to white
		sprintf( szColor, "rgb(%u,%u,%u)",255,color,color );

		// modify output
		output[0] = szColor; 
		output[2] = szColor; 
		output[3] = szColor; 
	}	


	return output; 
}

void macrophage_phenotype( Cell* pCell, Phenotype& p, double dt)
{
	// secrete anti-inflammatory in low O2 

	double o2 = get_single_signal( pCell, "oxygen"); 
	double halfmax = get_single_behavior( pCell , "custom:M1M2_halfmax"); // 5
	double hillpower = get_single_behavior( pCell , "custom:M1M2_hillpower"); // 4 
	double hill = Hill_response_function( o2 , halfmax, hillpower ); 
	double rate = get_single_base_behavior( pCell, "anti-inflammatory secretion")*(1-hill); 
	set_single_behavior( pCell , "anti-inflammatory secretion" , rate); 

	// secrete pro-inflammatory except in low O2 

	rate = get_single_base_behavior( pCell, "pro-inflammatory secretion")*hill; 
	set_single_behavior( pCell , "pro-inflammatory secretion" , rate); 

	return; 
}

void CD8_Tcell_phenotype( Cell* pCell, Phenotype& p, double dt)
{
	// anti-inflammatory reduces motiltiy 

	double anti = get_single_signal(pCell, "anti-inflammatory"); 
	double halfmax = get_single_behavior( pCell , "custom:CD8_halfmax"); // 0.5 
	double hillpower = get_single_behavior( pCell , "custom:CD8_hillpower"); // 8 
	double hill = Hill_response_function( anti , halfmax , hillpower ); 

	double param = get_single_base_behavior( pCell, "migration speed") * (1-hill); 
	set_single_behavior( pCell , "migration speed" , param ); 

	// anti-inflammatory reduces cell killing 

	param = get_single_base_behavior( pCell, "attack cancer") * (1-hill); 
	set_single_behavior( pCell , "attack cancer" , param ); 

	return; 
}
 