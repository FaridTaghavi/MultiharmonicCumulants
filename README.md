MultiharmonicCumulants_v2 _ 0.m

Author: Seyed Farid Taghavi
Date: April 1, 2021

The package is prepared by Mathematica 12.0 to generate all possible multiharmonic cumulants applicable in heavy-ion physics by employing the generating function method.
The package contains the following functions:

	1. c,
	2. cCorr,
	3. cQvec,
	4. cMean,
	5. cTable,
	6. ncQvec,
	7. ncMean,
	8. ncCorr,
	9. Nsigma2,
	10. Nsigma2Mean,
	11. Nsigma2Qvec,
	12. Nsigma2P.
	
To run, save the package together with a Mathematica notebook (example.nb for instance) in the same folder. By running the following lines in the notebook, the functions in the package are accessible:

	SetDirectory[NotebookDirectory[]];
	<<MultiharmonicCumulants_v2_ 0.m
As an example:
	
	in[]:= cCorr[{3, 1, 1, 1}, {4, 5, -6}, {3, 4, 5, 6}, corr]
	out[]= -corr[-3, 3] (corr[-6, -3, 4, 5] + corr[-5, -4, 3, 6]) + 1/2 (corr[-6, -3, -3, 3, 4, 5] + corr[-5, -4, -3, 3, 3, 6])
	
An explanation for each function can be obtained by adding a "?" at the beginning of the function's name. For instance,

	?cQvec
returns	

	cQvec[cumulant_order _list , phase_list, harmonic_list , multiplicity_symbole , Q-vector_symbol]: cumulant written in terms of Q-vectors.
