////////////////////////////////////////////////////////////
//
//       Pseudocode for the Exigere software
//
//         Gives a succinct overview of what goes on
//         and points to the locations of the relavent
//         code
//
///////////////////////////////////////////////////////////


exigere()
{

  read_parameter_file();   //read in the simulation parameters from a .prm file
                           //readParamfile() in readParamfile.c

  build_world();           //set up everything to do with the 
                           //simulated observations that does not change
                           //from event to event
                           //buildWorld() in buildWorld.c
  
  write_event_header();    //write header info to output file
                           //writeEventHdr() in info.c
 
  initialize_blending();   //Load in blending profiles
                           //initBlending() in findBlending.c

  read_synthetic_galaxy(); //Load in the synthetic galaxy data
                           //readSynthGal() in readSynthGal.c

  for each galaxy event    //each event in the galaxy is simulated, 
    {                      //then weighted

      build_event();          //Set up the event parameters
                              //buildEvent() in buildEvent.c

      time_sequencer();       //Work out when each observation will be
                              //timeSequencer() in timeSequencer.c

      find_blending();        //Calculate the blending for the event
                              //findBlending() in findBlending.c

      generate_lightcurve();  //Calculate the lightcurve
                              //lightcurveGenerator() in lightcurveGenerator.c

      if(event_on_detector)   //only do anything if the event was seen
	{
	  fit_event_with_PSPL();    //fit with a PSPL model
	                            //lightcurveFitter() in lightcurveFitter.c

	  if(u0<2rs)                //if finite source effects likely
	    {
	      fit_event_with_FSPL();  //fit with a FSPL model
                                      //lightcurveFitter_FS() in 
                                      //   lightcurveFitter_FS.c
	    }

	  output_lightcurve();      //Output the lightcurve for later use
	                            //outputLightcurve() in outputLightcurve.c

	  write_event_parameters(); //Write intresting data to file
	                            //writeEventParams() in info.c
	}

      write_event_stats_to_logfile(); //Write out statistics for a Monte Carlo 
                                      //integration of the survey area
                                      //performed inline
    }
}


build_world()
{

  for each observatory            //its possible to have multiple observatories
    {
      read_observatory_parameters();  //read in all their parameters
                                      //readObservatoryfile() in
                                      //   readObservatoryfile.c
    }

  find_collecting_areas();        //work out the telescope's effective 
                                  //collecting areas
                                  //collectingArea() in buildWorld.c

  if there are ground based observatories
    {
      find_sunrise_and_set_times();   //Compute the sunrise and set times for
                                      //each observatory
                                      //findSunRiseSetTimes() in buildWorld.c
                                        //uses functions in astroFns.c
    }

  compute_cadence();              //Work out times of all possible 
                                  //observations given sunrise/set etc
                                  //computeCadence() in buildWorld.c

  apply_weather();                //Work out when the weather allows 
                                  //observations - reads in a weather 
                                  //profile, which specifies probability of
                                  //bad weather, or e.g. spacecraft
                                  //downtime due to pointing restrictions
                                  //applyWeather() in buildWorld.c
 
  build_fields();                 //Load in the field centres and calculate 
                                  //their corners based on detector geometry
                                  //buildFields() in buildWorld.c
                                    //uses functions in astroFns.c

  compute_saturation_limits();    //computeSatLimit() in buildWorld.c
                                    //uses functions in astroFns.c

}

build_event()
{
  get_synthetic_galaxy_parameters();  //get the parameters set by the galaxy
  generate_parameters();              //generate the random parameters
}

time_sequencer()
{
  galactic_to_radec();        //Convert event galactic coordinates to ra/dec
                              //eq2gal() in astroFns.c

  which_fields_is_event_in(); //Work out which field the event is in
                              //whichFields() in timeSequencer.c
                             
  determine_epochs();         //Work out when the event will be observed
                              //determineEpochs() in timeSequencer.c
                                //uses functions in astroFns.c
}




/////////////////////////////////////////////////////////////////////
//
//    Other files:
//
/////////////////////////////////////////////////////////////////////

//definitions.h:    Holds hard coded simulation parameters - those that 
                    //affect memory usage etc

//blendHist_s?d?.h: //Defines the blending histograms. Number after s refers to
                    //the seeing, and after d the stellar density with the key:
                    //  s=1 --> 2.1"      d=1 --> 65.5 stars per sq arcmin 
                    //  s=4 --> 1.05"     d=3 --> 131 stars per sq arcmin 
                    //  s=7 --> 0.7"      d=5 --> 196.5 stars per sq arcmin 
                    //see Smith et al, 2007, MNRAS, 380, 805 for details

//constdefs.h:      //Defines pi etc

//constants.h:      //Defines astro constants e.g. G, Rsol, etc.

//mathFns.h/c       //Defines a few useful maths related functions

//random.h/c        //Random number generator

//strfns.h/c        //Useful string related functions

//structures.h      //Defines all the data structures - those that hold:
                    //   filekeywords - simulation parameters
                    //   obsfilekeywords - observatory parameters
                    //   event - event parameters
                    //   galaxy - synthetic galaxy parameters
                    //etc.

//*.f & global.h    //Fortran files for the calculation of FS binary lens 
                    //lightcurves



