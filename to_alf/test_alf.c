#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <X11/Xlib.h>
#include "t_o.h"
#include "epsilon.h"
#include "all.h"
#include "memfunc.h"

extern int t_o_domains_equal(t_o_domain* domain1, t_o_domain* domain2);

#define MARGIN        (30)
#define FRAME_WIDTH   (600)
#define FRAME_HEIGHT  (500)
#define WINDOW_WIDTH  (FRAME_WIDTH  + 2 * MARGIN)
#define WINDOW_HEIGHT (FRAME_HEIGHT + 3 * MARGIN)
#define RED           (0xFF0066)
#define BLUE          (0x0000FF)
#define GREEN         (0x00CC66)
#define GRAY          (0xB1B1B1) 
#define BLACK         (0x000000)
#define ONE_MINUTE    (60.0)
#define ONE_HOUR      (60.0 * ONE_MINUTE)
#define ONE_DAY       (24.0 * ONE_HOUR)

/* A link list struct stores depth, theta with depth increasing order. */
struct t_o_profile_link
{
  double depth;
  double theta;
  double pressure;
  struct t_o_profile_link *next;
};
typedef struct t_o_profile_link t_o_profile_link;

/* Insert a node with (depth, theta) into link list (head_pt), return error = TRUE if error occurs. */
int insert_t_o_profile_link(t_o_profile_link** head_pt, double depth, double theta, double pressure, int yes_groundater)
{
  int error = FALSE;
  assert(NULL != *head_pt);
  
  t_o_profile_link* pt_insert;
  t_o_profile_link* pt_rear;
  t_o_profile_link* pt_front;
  
  if (NULL == (pt_insert = (t_o_profile_link*)malloc(sizeof(t_o_profile_link))))
    {
      printf("malloc failed in creating t_o_profile_link node in function insert_t_o_profile_link(), exit.\n");
      error = TRUE;
      return error;
    }
  pt_insert->depth    = depth;
  pt_insert->theta    = theta;
  pt_insert->pressure = pressure;
  pt_insert->next     = NULL;
  
  pt_rear  = *head_pt;
  pt_front = (*head_pt)->next;
  // Head pointer *head_pt is initialized with domain->top_depth, and pass in value (depth) is larger than domain->top_depth.
  while (NULL != pt_front)
    {
      if (pt_insert->depth >= pt_front->depth)
        {
          pt_rear  = pt_front;
          pt_front = pt_front->next;     
        }
      else
        {
          break;
        }
    }
  
  // Insert pt_insert only when depth is strictly larger than pt_rear->depth, so that depth in link is increasing without same value.
  if (epsilon_greater(pt_insert->depth, pt_rear->depth))
    {
      pt_rear->next   = pt_insert;
      pt_insert->next = pt_front;
    }
  else if (!yes_groundater && epsilon_equal(pt_insert->depth, pt_rear->depth) && pt_rear->theta < pt_insert->theta)
    { // Add 10/09/14, to handle no groundwate cases,
     //printf("theta = %lf \n", pt_insert->theta); getchar();
      pt_rear->theta    = pt_insert->theta;
      pt_rear->pressure = pt_insert->pressure;
    }
  return error;
}

/* Map t_o domain into water content in 1D space discretization, water_content[num_elements] and effective_porosity are output values.
 *
 * Parameters:
 * domain             - A pointer to the t_o_domain struct.
 * num_elements       - Number of element in 1d soil column.
 * soil_depth_z       - A pointer to 1d array of size num_elements contains depth of each element's lower bound, in unit of [meters], positive downward.
 *
 * Output:
 * water_content      - A pointer to 1d array of size num_elements contains water content of each element.
 * pressure_head      - A pointer to 1d array of size num_elements contains pressure head of each element.
 * effective_porosity - Scalar passed by reference, unitless.
 */
int get_t_o_domain_profile(t_o_domain* domain, int num_elements, double* soil_depth_z, 
                           double* water_content,  double* pressure_head, double* effective_porosity)
{
  int error = FALSE;
  int ii, jj, jj_start;
  t_o_profile_link *head_pt = NULL, *tmp_pt = NULL, tmp;
  
  // Initialize head pointer node, fake node.
  tmp.depth    = domain->layer_top_depth;
  tmp.theta    = domain->parameters->bin_water_content[domain->parameters->num_bins];
  tmp.pressure = 0.0;
  tmp.next     = NULL;
  head_pt      = &tmp;
  
  assert(NULL != domain);
  
  // Find effective porosity.
  *effective_porosity = domain->parameters->bin_water_content[domain->parameters->num_bins];
  
  // Step 1, map the to_domain water content into link list t_o_profile_link_link.
  // Loop over bins, start from bin 2, since bin 1 is always saturated.
  for (ii = 2; ii <= domain->parameters->num_bins; ii++)
     {
       // Surface front.
       if (domain->surface_front[ii] > domain->layer_top_depth && domain->surface_front[ii] <= domain->layer_bottom_depth) // Add "<=" 10/09/14.
         {
           error = insert_t_o_profile_link(&head_pt, domain->surface_front[ii], domain->parameters->bin_water_content[ii], 
                                                                               -domain->parameters->bin_capillary_suction[ii], domain->yes_groundwater);
         }
        
       // Slug.
       slug* temp_slug = domain->top_slug[ii];
       while (NULL != temp_slug)
         {
           // slug top.
           error = insert_t_o_profile_link(&head_pt, temp_slug->top, domain->parameters->bin_water_content[ii - 1], 
                                                                    -domain->parameters->bin_capillary_suction[ii - 1], domain->yes_groundwater);
           
           // slug bottom.
           error = insert_t_o_profile_link(&head_pt, temp_slug->bot, domain->parameters->bin_water_content[ii], 
                                                                   -domain->parameters->bin_capillary_suction[ii], domain->yes_groundwater);
           
           temp_slug = temp_slug->next;
         }
        
       // groundwater front. 
       if (domain->yes_groundwater)
         {
           if (domain->groundwater_front[ii] > domain->layer_top_depth)
             {
              error = insert_t_o_profile_link(&head_pt, domain->groundwater_front[ii], domain->parameters->bin_water_content[ii - 1], 
                                                                                 -domain->parameters->bin_capillary_suction[ii - 1], domain->yes_groundwater);
             }
         }
     } // End of bin loop.
  
  // Insert fully saturated bin;
  if (domain->yes_groundwater)
    {
      error = insert_t_o_profile_link(&head_pt, domain->layer_bottom_depth, domain->parameters->bin_water_content[domain->parameters->num_bins], 
                                                          0.0, domain->yes_groundwater);
    }
  else
    {
      int first_bin = 2; // The leftmost bin that is not completely full of water.
 
      while (first_bin <= domain->parameters->num_bins && domain->parameters->bin_water_content[first_bin] <= domain->initial_water_content)
       {
         first_bin++;
       }
      error = insert_t_o_profile_link(&head_pt, domain->layer_bottom_depth, domain->parameters->bin_water_content[first_bin - 1], 
                                                                         -domain->parameters->bin_capillary_suction[first_bin - 1], domain->yes_groundwater);
    }

  // Step 2, fill 1D array water content using t_o_profile_link_link.
  tmp_pt   = head_pt;
  jj_start = 1;
  while (NULL != tmp_pt)
    {
      for (jj = jj_start; jj <= num_elements; jj++)
         {
           if (soil_depth_z[jj] <= tmp_pt->depth)
             {
               water_content[jj] = tmp_pt->theta;
               pressure_head[jj] = tmp_pt->pressure;
             }
           else
             {
               jj_start = jj;
               break;
             }
         }
      tmp_pt = tmp_pt->next;
    } 
    
  // Set head_pt to head_pt->next, as the first node is not malloc. Then free link-list.
  head_pt = head_pt->next;
  while (NULL != head_pt)
    {
       tmp_pt  = head_pt;
       head_pt = head_pt->next;
       free(tmp_pt);
    }
  
  return error;
} // End of get_t_o_domain_profile().

int main(void)
{
  t_o_parameters* parameters;
  t_o_domain*     domain;
	
		int    ii, jj, error = FALSE;
    	int    num_bins              = 100;//#0;           	// Number of bins.   
      	double conductivity          = 5.0/86400; 		// Meters per second.
      	double porosity              = 0.45;     			// Unitless fraction.
      	double residual_saturation   = 0.07;    			// Unitless fraction.
      	int    van_genutchen         = TRUE;     			// Yes, use Van Genuchten.
      	double vg_alpha              = 0.036*100;      	// One over meters.
      	double vg_n                  = 3.5;     			// Unitless.
      	double bc_lambda             = 5.5;    			// Unitless.
      	double bc_psib               = 0.37;     			// Meters.
      	double layer_top_depth       = 0.0;      			// Meters.
      	double layer_bottom_depth    = 1.50;      		// Meters.  
	    double initial_tension_top = 205;				//centimeters mbar; used for setting initial condition
	
	  int buf_alloc                = 1000;
		  
      int    yes_groundwater       = FALSE;     // SImulate Groundwater?
      int    yes_runoff            = FALSE;     // Remove excess surface water after infilt step.?
      double water_table           = layer_bottom_depth;      // Meters.
      
      double current_time          = 0.0;                                                     // Current time in seconds.
<<<<<<< HEAD
      double delta_time            = 3600;                                                     // The duration of the timestep in seconds.
      double max_time              = 100 * ONE_HOUR;                                         // How long to run the simulation in seconds.
      double runoff,frate,rech;
=======
      double delta_time            = 900;                                                     // The duration of the timestep in seconds.
      double max_time              = 24 * ONE_HOUR;                                         // How long to run the simulation in seconds.
      double runoff = 0.0;
	  //double frate = 0.0;
	  //double rech = 0.0;
>>>>>>> 42a0f9357b90890ffdd54a509014deaf81fd8ca9
      
      double infiltration_rate     = 0.0;
      double groundwater_inf_rate = 0.0;
      double rainfall_input_time[buf_alloc];
      double rainfall_input_intensity[buf_alloc];
      double potential_ET[buf_alloc];
      char string[buf_alloc];
	  double params[buf_alloc];
      
      double evaporated_water      = 0.0;
      double surfacewater_depth    = 0.0;                                                     // Meters.
      double total_water           = 0.0;                                                     // The should be value for total water in meters of water.
      double groundwater_recharge  = 0.0;                                                     // The total amount supplied to groundwater in meters of water.
      double accum_infil           = 0.0;
      double surfacewater_depth_old   = 0.0;
      double groundwater_recharge_old = 0.0;
      double domain_initial_water = 0.0;
      double accu_PET             = 0.0;
      double accu_rain            = 0.0;
      double PET                  = 0.0;
      double rainfall_rate        = 0.0;
      double rainfall             = 0.0;
  	  
      //Load parameters from flat file 
	  FILE*  param_fptr      = NULL;
      if (NULL == (param_fptr = fopen("TO.IN", "r")))
        {
          printf("Error reading TO.IN \n");
          exit(1);
        }
	
      fgets(string, 100, param_fptr);      // Ignore header.
	  for (ii = 1; ii <= 20 ; ii++) {
		  fscanf(param_fptr, "%lf", &params[ii]);  // Original data in mm /15 min.
		  printf("Param %d %lf\n",ii,params[ii]); 
	  }
      fclose(param_fptr);
	
	  num_bins              = params[1];           		// Number of bins.   
      conductivity          = params[2]; 				// Meters per second.
      porosity              = params[3];     			// Unitless fraction.
      residual_saturation   = params[4];    			// Unitless fraction.
      van_genutchen         = params[5];     			// Yes, use Van Genuchten.
      vg_alpha              = params[6];      			// One over meters.
      vg_n                  = params[7];     			// Unitless.
      bc_lambda             = params[8];    			// Unitless.
      bc_psib               = params[9];     			// Meters.
      layer_top_depth       = params[10];      			// Meters.
      layer_bottom_depth    = params[11];    			// Meters.
	  initial_tension_top 	= params[12];				// cm H2O. Greater than zero. Used to define initial water content in terms of vGM.
	
	
double initial_water_content = pow(1/pow((vg_alpha*initial_tension_top/100),vg_n),1-(1/vg_n))*porosity;
	
  int    test_id  = 1;  
  
  
<<<<<<< HEAD
  int    test_id  = 1;  // ###################################################################################################################
  /*
  test_id =  1, origianl panama test.
  test_id =  2, sand infiltration.
  test_id = 101 - 112, tests using 12 USAD soil type with van Genucthen parameters.
  */
  if ( 1 == test_id)
    {
      num_bins              = 300;            // Number of bins.
      conductivity          = 0.3 / 86400; // Meters per second.
      porosity              = 0.4;     // Unitless fraction.
      residual_saturation   = 0.027;    // Unitless fraction.
      van_genutchen         = TRUE;  // Yes, use Van Genuchten.
      vg_alpha              = 3.6;      // One over meters.
      vg_n                  = 22;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = FALSE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = FALSE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.1;      // Unitless fraction.
      water_table           = 1*layer_bottom_depth;      // Meters.
      
      delta_time            = 600;                                                     // The duration of the timestep in seconds.
      max_time              = 72 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 2 == test_id)
    { // Sand.
      num_bins              = 300;            // Number of bins.
      conductivity          = 29.7 / 360000.0; // Meters per second.
      porosity              = 0.43;     // Unitless fraction.
      residual_saturation   = 0.045;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 14.5;      // One over meters.
      vg_n                  = 2.68;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.

      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 1.0;                                                     // The duration of the timestep in seconds.
      max_time              = 6.0 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 3 == test_id)
    { // Silt
      num_bins              = 400;            // Number of bins.
      conductivity          = 0.25 / 360000.0; // Meters per second.
      porosity              = 0.46;     // Unitless fraction.
      residual_saturation   = 0.034;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 1.6;      // One over meters.
      vg_n                  = 1.37;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 1.0;                                                     // The duration of the timestep in seconds.
      max_time              = 6.0 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 101 == test_id)
    { // Panama test using Sand. 
      conductivity          = 29.7 / 360000.0; // Meters per second.
      porosity              = 0.43;     // Unitless fraction.
      residual_saturation   = 0.045;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 14.5;      // One over meters.
      vg_n                  = 2.68;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 102 == test_id)
    {
      conductivity          = 14.5917 / 360000.0; // Meters per second.
      porosity              = 0.41;     // Unitless fraction.
      residual_saturation   = 0.057;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 12.4;      // One over meters.
      vg_n                  = 2.28;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 103 == test_id)
    {
      conductivity          = 4.42083 / 360000.0; // Meters per second.
      porosity              = 0.41;     // Unitless fraction.
      residual_saturation   = 0.065;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 7.5;      // One over meters.
      vg_n                  = 1.89;     // Unitless.
      bc_lambda             = 5.5;    // Yeptless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 104 == test_id)
    {
      conductivity          = 1.04 / 360000.0; // Meters per second.
      porosity              = 0.43;     // Unitless fraction.
      residual_saturation   = 0.078;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 3.6;      // One over meters.
      vg_n                  = 1.56;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }  
  else if ( 105 == test_id)
    {
      conductivity          = 0.25 / 360000.0; // Meters per second.
      porosity              = 0.46;     // Unitless fraction.
      residual_saturation   = 0.034;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 1.6;      // One over meters.
      vg_n                  = 1.37;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 106 == test_id)
    {
      conductivity          = 0.45 / 360000.0; // Meters per second.
      porosity              = 0.45;     // Unitless fraction.
      residual_saturation   = 0.067;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 2.0;      // One over meters.
      vg_n                  = 1.41;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 107 == test_id)
    {
      conductivity          = 1.31 / 360000.0; // Meters per second.
      porosity              = 0.39;     // Unitless fraction.
      residual_saturation   = 0.1;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 5.9;      // One over meters.
      vg_n                  = 1.48;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 108 == test_id)
    {
      conductivity          = 0.26 / 360000.0; // Meters per second.
      porosity              = 0.41;     // Unitless fraction.
      residual_saturation   = 0.095;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 1.9;      // One over meters.
      vg_n                  = 1.31;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 109 == test_id)
    {
      conductivity          = 0.07 / 360000.0; // Meters per second.
      porosity              = 0.43;     // Unitless fraction.
      residual_saturation   = 0.089;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 1.0;      // One over meters.
      vg_n                  = 1.23;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 110 == test_id)
    {
      conductivity          = 0.12 / 360000.0; // Meters per second.
      porosity              = 0.38;     // Unitless fraction.
      residual_saturation   = 0.1;      // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 2.7;      // One over meters.
      vg_n                  = 1.23;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 111 == test_id)
    {
      conductivity          = 0.02 / 360000.0; // Meters per second.
      porosity              = 0.36;     // Unitless fraction.
      residual_saturation   = 0.07;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 0.5;      // One over meters.
      vg_n                  = 1.09;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  else if ( 112 == test_id)
    {
      conductivity          = 0.2 / 360000.0; // Meters per second.
      porosity              = 0.38;     // Unitless fraction.
      residual_saturation   = 0.068;    // Unitless fraction.
      van_genutchen         = TRUE;     // Yes, use Van Genuchten.
      vg_alpha              = 0.8;      // One over meters.
      vg_n                  = 1.09;     // Unitless.
      bc_lambda             = 5.5;    // Unitless.
      bc_psib               = 0.37;     // Meters.
      layer_top_depth       = 0.0;      // Meters.
      layer_bottom_depth    = 1.0;      // Meters.     ##################################
      yes_groundwater       = TRUE;     // Yes, simulate groundwater. #################################################################
      yes_runoff            = TRUE;     // Yes, remove excess surface water after infilt step.
      initial_water_content = 0.08;      // Unitless fraction.
      water_table           = layer_bottom_depth;      // Meters.
      
      delta_time            = 10.0;                                                     // The duration of the timestep in seconds.
      max_time              = 5750 * ONE_HOUR;                                         // How long to run the simulation in seconds.
    }
  
=======
>>>>>>> 42a0f9357b90890ffdd54a509014deaf81fd8ca9
  if (t_o_parameters_alloc(&parameters, num_bins, conductivity, porosity, residual_saturation, van_genutchen,  vg_alpha, vg_n, bc_lambda, bc_psib))
    {
      fprintf(stderr, "ERROR: Could not allocate t_o_parameters.\n");
      exit(1);
    }

  if (t_o_domain_alloc(&domain, parameters, layer_top_depth, layer_bottom_depth, yes_groundwater, initial_water_content, yes_groundwater, water_table))
    {
      fprintf(stderr, "ERROR: Could not allocate t_o_domain.\n");
      exit(1);
    }

  t_o_check_invariant(domain);

    /***********************/
   /* Run the simulation. */
  /***********************/
  if (1 == test_id)
    {
      // Read input rainfall file.
      FILE*  rain_fptr      = NULL;
  
      if (NULL == (rain_fptr = fopen("SURFACE_BC.IN", "r")))
        {
          printf("Error reading SURFACE_BC.IN \n");
          exit(1);
        }
        
      fgets(string, 100, rain_fptr);      // Ignore header.
      for (ii = 1; ii <= buf_alloc; ii++)
         {
           fscanf(rain_fptr, "%*d %lf %lf", &rainfall_input_intensity[ii], &potential_ET[ii]);  // Original data in mm /15 min.
           rainfall_input_time[ii]       = (ii - 1) * 15.0 * 60.0;  // In s.
           rainfall_input_intensity[ii] /= 900000.0;     // In m/s. 
           //rainfall_input_intensity[ii] *= 5.0;           // convert back to real values.
           potential_ET[ii]             /= 900000.0;     // In m/s.
      	   //printf("IN: %lf",rainfall_input_intensity[ii]);
	  }
      fclose(rain_fptr);
    }

  FILE*  f_fptr;                                                                          // File output of infiltration rate.
  FILE*  acc_depth_fptr;
  FILE* fptr_obsnode_60;
  FILE* fptr_obsnode_150;
  FILE* fptr_surfbc;
  //double acc_depth = 0.0;
  
  if (NULL == (fptr_obsnode_60 = fopen("alf_obsnode_60.txt", "w")) )
  {
	printf("ERROR: Could not open observation node output file.\n");
	exit(1);
  }  

  if (NULL == (fptr_obsnode_150 = fopen("alf_obsnode_150.txt", "w")) )
  {
	printf("ERROR: Could not open observation node output file.\n");
	exit(1);	
  }

  if (NULL == (fptr_surfbc = fopen("surfbc.out", "w")) )
  {
	printf("ERROR: Could not open observation node output file.\n");
	exit(1);	
  }
	
  if ((f_fptr = fopen("f.out", "w")) == NULL)
    {
     fprintf(stderr, "ERROR: Could not open infiltration output file.\n");
     exit(1);
    }
  if ((acc_depth_fptr = fopen("accum_depth.out", "w")) == NULL)
    {
     fprintf(stderr, "ERROR: Could not open accumulate depth output file.\n");
     exit(1);
    }
  
  int num_layers   = 1;
  int num_elements = 1000;
  double** soil_depth_z;   
  double** water_content;   
  double** pressure_head;   
  double effective_porosity;
  
  error =  dtwo_alloc(&soil_depth_z,  num_layers, num_elements);
  error =  dtwo_alloc(&water_content, num_layers, num_elements);
  error =  dtwo_alloc(&pressure_head, num_layers, num_elements);
	  
  for (ii = 1; ii <= num_layers ; ii++)
     {
       soil_depth_z[ii][0] = domain->layer_top_depth;
       for (jj = 1; jj <= num_elements; jj++)
          {
            soil_depth_z[ii][jj]  = soil_depth_z[ii][0] + jj * (domain->layer_bottom_depth  - domain->layer_top_depth) / num_elements;
          }
     }
  
  time_t time_start       = time(NULL); // Wall clock time.
  time_t time_end;                      // Wall clock time.
  domain_initial_water = t_o_total_water_in_domain(domain);
  accu_PET             = 0.0;
  accu_rain            = 0.0;
  PET                  = 0.0;
  rainfall_rate        = 0.0;
  rainfall             = 0.0;
  runoff               = 0.0;
  total_water          = surfacewater_depth + domain_initial_water + accu_PET + groundwater_recharge + accu_rain;
	
  while (current_time < max_time)
    {
      // Rainfall.  ###########################################################
      if (1 == test_id)
        {
          ii = 1;
          while (ii < buf_alloc && current_time >= rainfall_input_time[ii + 1])
            {
             ii++;
            }
          PET           = potential_ET[ii];                  // m/s.      ##################################### Set no PET.
          rainfall_rate = rainfall_input_intensity[ii];      // m/s.
          rainfall      = rainfall_rate * delta_time;        // Meters of water.
      	  //printf("Rainfall %lf\n",rainfall);
	  }
      //printf("Add %lf \n",rainfall);
	  surfacewater_depth  += rainfall;
      total_water         += rainfall;
      accu_rain           += rainfall;
      accu_PET            += PET * delta_time;
      
      surfacewater_depth_old   = surfacewater_depth;
      groundwater_recharge_old = groundwater_recharge;
      // Moving water table
     // water_table = 1.0 - 0.7 * fabs(sin(2.0 * 3.14 * current_time / (24.0 * 3600.0)));  //
      
      if (t_o_timestep(domain, delta_time, surfacewater_depth, &surfacewater_depth, water_table, &groundwater_recharge))
        {
          fprintf(stderr, "ERROR: t_o_timestep returned error.\n");
          exit(1);
        }
      
      // FIXME, Add ET as seperate function, Dec. 10, 2014. Better put them inside t_o_timestep.
      //if (1 == test_id)
      //  {
      //   // t_o_ET(domain, delta_time, root_depth, PET, field_capacity, wilting_ponint, 
                    // use_feddes, field_capacity_suction, wilting_ponint_suction, &surfacewater_depth, &evaporated_water);
      //    if (t_o_ET(domain, delta_time, 0.5,       PET, 0.32,           0.03,           TRUE, 0.27, 1527.68, &surfacewater_depth, &evaporated_water))
       //     {
        //      fprintf(stderr, "ERROR: t_o_timestep returned error.\n");
         //     exit(1);
         //   }
        //} End of ET for test 1.
	  
      infiltration_rate    = (surfacewater_depth_old - surfacewater_depth) / delta_time * 100.0 * ONE_HOUR; // cm/hr.
      groundwater_inf_rate = (groundwater_recharge - groundwater_recharge_old) / delta_time * 100 * ONE_HOUR; // cm/hr.
      //fprintf(f_fptr, "%lf %lf %lf %lf %lf %lf %lf %lf \n", current_time, rainfall_rate * 360000.0, accu_rain * 100, infiltration_rate, accum_infil * 100.0,
      //                                                     groundwater_inf_rate, groundwater_recharge * 100.0, evaporated_water * 100.0);

	  double frate, rech;
      frate = (surfacewater_depth_old - surfacewater_depth) / delta_time * 100.0 * ONE_HOUR; // cm/hr.
      rech  = (groundwater_recharge - groundwater_recharge_old) / delta_time * 100.0 * ONE_HOUR; // cm/hr.
      //printf("%lf %lf %lf %lf %lf \n", current_time/86400, rainfall_rate * 360000.0, frate, rech, surfacewater_depth);
  	  if(surfacewater_depth > 0) {
		  	printf("%lf",surfacewater_depth);
	  }
      current_time += delta_time;
      
      if(yes_runoff)  // new FLO 21 Dec. 2014 
        {
          runoff            += surfacewater_depth;
          surfacewater_depth = 0.0;
        }
     
     assert(epsilon_equal(total_water, evaporated_water + surfacewater_depth + groundwater_recharge + t_o_total_water_in_domain(domain) + runoff));
      
      fprintf(acc_depth_fptr, "%lf %lf\n", current_time, runoff * 100.0); // Time in s, runoff in cm.
      //fprintf(acc_depth_fptr,"%lf %lf\n",current_time/86400.0, t_o_total_water_in_domain(domain)); 
     fprintf(fptr_surfbc,"%lf\n",surfacewater_depth);
	  // obsERVATION NODE OUTPUT - AB
      get_t_o_domain_profile(domain, num_elements, soil_depth_z[1], water_content[1], pressure_head[1], &effective_porosity);
      for (jj = 1; jj <= num_elements; jj++) {
      	if(floor(soil_depth_z[1][jj]==0.6)) {
			fprintf(fptr_obsnode_60, "%lf %lf %lf %lf %lf %lf\n", current_time, soil_depth_z[1][jj], water_content[1][jj], pressure_head[1][jj],water_content[1][1],pressure_head[1][1]);
        } else if (floor(soil_depth_z[1][jj]==1.5))  {
			fprintf(fptr_obsnode_150, "%lf %lf %lf %lf\n", current_time, soil_depth_z[1][jj], water_content[1][jj], pressure_head[1][jj]);		    	
		}
	  }
	        
    } // End of time loop.
  time_end = time(NULL);
  FILE* fptr_simout = NULL;
  if (NULL == (fptr_simout = fopen("TO.OUT", "w")) )
  {
	printf("ERROR: Could not open simulation output file.\n");
	exit(1);
  }  
  printf("\nTotal simulation time  = %lf hours\n", max_time / ONE_HOUR);
  printf("\nElapsed time = %Lf seconds \n \n",(long double)(time_end - time_start));
  printf("Mass balance info: \n");
  printf("Initial water in domain  = %lf mm \n", domain_initial_water * 1000);
  printf("Accumulated rainfall     = %lf mm \n", accu_rain * 1000);
  printf("Accumulated PET          = %lf mm \n", accu_PET * 1000);
  printf("Accumulated AET          = %lf mm \n", evaporated_water * 1000);
  printf("Accumulated infiltration = %lf mm \n", accum_infil * 1000);
  printf("Groundwater recharge     = %lf mm \n", groundwater_recharge * 1000);
  printf("Final water in domain    = %lf mm \n", t_o_total_water_in_domain(domain) * 1000);
  printf("Final surface water      = %lf mm \n", surfacewater_depth);
  printf("Total surface runoff     = %lf mm \n", runoff*1000.0);
  printf("Mass error               = %8.5e mm \n", (domain_initial_water + accu_rain - evaporated_water - groundwater_recharge - 
                                                 t_o_total_water_in_domain(domain) - surfacewater_depth-runoff) * 1000);
  printf("Number of bins: %d \n", num_bins);
  printf("Initial tension: %lf \n",initial_tension_top);
  printf("Initial water content: %lf \n",initial_water_content);
  printf("Ksat: %lf meters/s \n",conductivity);
  printf("n: %lf \n",vg_n);
  printf("alpha: %lf \n",vg_alpha);
	
	
  fprintf(fptr_simout,"\nTotal simulation time  = %lf hours\n", max_time / ONE_HOUR);
  fprintf(fptr_simout,"\nElapsed time = %Lf seconds \n \n",(long double)(time_end - time_start));
  fprintf(fptr_simout,"Mass balance info: \n");
  fprintf(fptr_simout,"Initial water in domain  = %lf mm \n", domain_initial_water * 1000);
  fprintf(fptr_simout,"Accumulated rainfall     = %lf mm \n", accu_rain * 1000);
  fprintf(fptr_simout,"Accumulated PET          = %lf mm \n", accu_PET * 1000);
  fprintf(fptr_simout,"Accumulated AET          = %lf mm \n", evaporated_water * 1000);
  fprintf(fptr_simout,"Accumulated infiltration = %lf mm \n", accum_infil * 1000);
  fprintf(fptr_simout,"Groundwater recharge     = %lf mm \n", groundwater_recharge * 1000);
  fprintf(fptr_simout,"Final water in domain    = %lf mm \n", t_o_total_water_in_domain(domain) * 1000);
  fprintf(fptr_simout,"Final surface water      = %lf mm \n", surfacewater_depth);
  fprintf(fptr_simout,"Total surface runoff     = %lf mm \n", runoff*1000.0);
  fprintf(fptr_simout,"Mass error               = %8.5e mm \n", (domain_initial_water + accu_rain - evaporated_water - groundwater_recharge - 
                                                 t_o_total_water_in_domain(domain) - surfacewater_depth-runoff) * 1000);

    /*************
   /* Clean up. */
  /*************/
  t_o_domain_dealloc(&domain);
  t_o_parameters_dealloc(&parameters);
  fclose(fptr_obsnode_60);
  fclose(fptr_obsnode_150);
  fclose(f_fptr);
  fclose(acc_depth_fptr);
  fclose(fptr_simout);
  fclose(fptr_surfbc);
  return 0;
}
