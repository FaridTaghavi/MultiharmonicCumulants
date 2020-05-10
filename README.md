# MultiharmonicCumulants_v1.1

Author: Seyed Farid Taghavi,

Technical University of Munich,
Dense and Strange Hadronic Matter Group,

Date: May 9, 2020,

The package is prepared by Mathematica 12.0 to generate all possible multi-harmonic cumulants applicable in heavy-ion physics by employing the generating function method.
It contains the following function:

	c[ cumulant_order _list, phase_list , harmonic_list , flow_amplitude _symbole , flow_phase _symbol ]            
	 --->    cumulant written in terms of symbolic moments.
	
	cCorr[ cumulant_order _list, phase_list , harmonic_list , correlation_symbol ]                            
	 --->    cumulant written in terms of correlations.
	
	cTable[ harmonic_list , maximum_order , flow_amplitude _symbole , flow_phase _symbol ]                                                
	--->    all cumulants with order between 2 to maximum_order in terms of symbolic moments.
	
	cTable[ harmonic_list , minimum_order , maximum_order , flow_amplitude _symbole , flow_phase _symbol ]                               
	 --->   all cumulants with order between minimum_order to maximum_order written in terms of symbolic moments.
	 
	Nsigma2Moment[ correlated_harmonics _list , flow_amplitude _symbole , flow_phase _symbole , multiplicity_symbol ]                      
	 --->   Number of events multiplied by statistical uncertainty square of correlator <<correlated_harmonic_list>> written in terms of symbolic moments.
	 
	Nsigma2Moment[ moment_order _list, phase_list , harmonic_list , flow_amplitude _symbole , flow_phase _symbol , multiplicity_symbol ]
	 --->   Number of events multiplied by statistical uncertainty square of correlator moment[moment_order _list, phase_list , harmonic_list] written in terms of symbolic moments.
	 
	Nsigma2MomentCorr[ correlated_harmonics _list , correlation_symbol , multiplicity_symbol ]                                            
	--->    Number of events multiplied by statistical uncertainty square of correlator <<correlated_harmonic_list>> written in terms of correlat
