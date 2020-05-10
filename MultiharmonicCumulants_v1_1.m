(* ::Package:: *)

(*  
MultiharmonicCumulants_v1_1.m

Author: Seyed Farid Taghavi
Date: May 8, 2020

The package is prepared by Mathematica 12.0 to generate all possible multi-harmonic cumulants applicable in heavy-ion physics by employing the generating function method.
It contains the following function:

	c[ cumulant_order_list, phase_list , harmonic_list , flow_amplitude_symbole , flow_phase_symbol ]            
	 --\[Rule]    cumulant written in terms of symbolic moments.
	
	cCorr[ cumulant_order_list, phase_list , harmonic_list , correlation_symbol ]                            
	 --\[Rule]    cumulant written in terms of correlations.
	
	cTable[ harmonic_list , maximum_order , flow_amplitude_symbole , flow_phase_symbol ]                                                
	--\[Rule]    all cumulants with order between 2 to maximum_order in terms of symbolic moments.
	
	cTable[ harmonic_list , minimum_order , maximum_order , flow_amplitude_symbole , flow_phase_symbol ]                               
	 --\[Rule]   all cumulants with order between minimum_order to maximum_order written in terms of symbolic moments.
	 
	Nsigma2Moment[ correlated_harmonics_list , flow_amplitude_symbole , flow_phase_symbole , multiplicity_symbol ]                      
	 --\[Rule]   Number of events multiplied by statistical uncertainty square of correlator <<correlated_harmonic_list>> written in terms of symbolic moments.
	 
	Nsigma2Moment[ moment_order_list, phase_list , harmonic_list , flow_amplitude_symbole , flow_phase_symbol , multiplicity_symbol ]
	 --\[Rule]   Number of events multiplied by statistical uncertainty square of correlator moment[moment_order_list, phase_list , harmonic_list] written in terms of symbolic moments.
	 
	Nsigma2MomentCorr[ correlated_harmonics_list , correlation_symbol , multiplicity_symbol ]                                            
	--\[Rule]    Number of events multiplied by statistical uncertainty square of correlator <<correlated_harmonic_list>> written in terms of correlations.
	
*)


BeginPackage["MultiharmonicCumulants`"]


Unprotect[
	c,
	cCorr,
	cTable,
	Nsigma2Moment,
	Nsigma2MomentCorr
]


ClearAll[
	c,
	cCorr,
	cTable,
	Nsigma2Moment,
	Nsigma2MomentCorr
]


c::usage=     
   "c[ cumulant_order_list, phase_list , harmonic_list , flow_amplitude_symbole , flow_phase_symbol]: cumulant written in terms of symbolic moments "
cCorr::usage=     
   "cCorr[ cumulant_order_list, phase_list , harmonic_list , correlation_symbol]: cumulant written in terms of correlation functions"
cTable::usage=     
   "cTable[ harmonic_list , maximum_order , flow_amplitude_symbole , flow_phase_symbol] or cTable[ harmonic_list , minimum_order , maximum_order , flow_amplitude_symbole , flow_phase_symbol]: table of all cumulants written in terms of symbolic moments"
Nsigma2Moment::usage= 
   "Nsigma2Moment[ correlated_harmonics_list , flow_amplitude_symbole , flow_phase_symbole , Multiplicity_Symbol ] or Nsigma2Moment[ moment_order_list, phase_list , harmonic_list , flow_amplitude_symbole , flow_phase_symbol , Multiplicity_Symbol ]: Number of events multiplied by statistical uncertainty square"
Nsigma2MomentCorr::usage=   
	"Nsigma2MomentCorr[ correlated_harmonics_list , correlation_symbol , Multiplicity_Symbol ]:  Number of events multiplied by statistical uncertainty square"


Begin["`Private`"]


(* Useful Functions For Generating Cumulant*)

integrateFunc[term_]:=Module[{output,faze,Dfaze,fazeAlpha,fazeBeta,coeff}, (* We use a simple substitution instead of integrating to increas the computation speed. *)
	faze=Exponent[term,E];
	Dfaze=D[faze,fi];
	If[faze==0 || Dfaze==0,
	Return[2\[Pi] term];
	Exit[];
	];
	coeff=term/E^faze//Simplify;
	fazeAlpha=(Dfaze)/.{fi->0};
	fazeBeta=faze/.{fi->0};
	output=coeff ((E^fazeBeta (-1+E^(2 \[Pi] fazeAlpha)))/fazeAlpha)//Simplify
]
reduceToAverage[term_,rList_,fazList_]:=Module[{fazExpList,rExpList,output},
	fazExpList=Exponent[term,fazList];
	rExpList=Exponent[term,rList];
	output=term/((Times@@(fazList^fazExpList))(Times@@(rList^rExpList))) av[rExpList,fazExpList]//Simplify;
	Return[output]
]
reduceToAverageNice[term_,rList_,fazList_,HarmonicList_,AmpSymbol_,phaseSymbol_]:=Module[{fazExpList,rExpList,output,amp,\[CapitalDelta]},(* This function replaces combination of flow amplitudes with a symbolic moment. The moment is presented with brackets and flow amplirude and phase symbols with inputs *)
	fazExpList=Exponent[term,fazList];
	rExpList=Exponent[term,rList];
	(*The following is just cosmetic*)
	amp=Product[Subscript[AmpSymbol, HarmonicList[[i]]]^rExpList[[i]],{i,1,HarmonicList//Length}];
	\[CapitalDelta]=Sum[fazExpList[[i]](Subscript[phaseSymbol, HarmonicList[[i+1]]]-Subscript[phaseSymbol, HarmonicList[[1]]]),{i,1,fazExpList//Length}]//Simplify;
	(*The above is just cosmetic*)
	output=term/((Times@@(fazList^fazExpList))(Times@@(rList^rExpList))) AngleBracket[amp E^(I \[CapitalDelta])  ]//Simplify;
Return[output]
]
NonvanishingCombination[nList_,maxord_]:=Module[{pList={},qList={},p,q,output={},outout={},appendList={},negativeCheck={}},  (* Finding all non-vanishing combinations*)
	pList=Tuples[Range[0,maxord],nList//Length];
	qList=pList;
	Do[
		If[Plus@@pList[[ip]]==0 || Plus@@qList[[iq]]==0,Continue[]];
		If[(Plus@@pList[[ip]])+( Plus@@qList[[iq]])>maxord,Continue[]];
		If[AnyTrue[pList[[ip]]+qList[[iq]],(#1==0)&],Continue[]];
		If[Plus@@(nList(pList[[ip]]-qList[[iq]]))==0,
			appendList=(pList[[ip]]+qList[[iq]])~Join~Flatten[{Drop[nList(pList[[ip]]-qList[[iq]]),1]}];
			negativeCheck=(pList[[ip]]+qList[[iq]])~Join~(-Flatten[{Drop[nList(pList[[ip]]-qList[[iq]]),1]}]);
			If[((Position[output,negativeCheck])//Length)!=0,Continue[]];
			output=Append[output,appendList];
		];
	,{ip,1,pList//Length},{iq,1,qList//Length}];
	outout=Sort[DeleteDuplicates[output],(Sum[#1[[ii]],{ii,1,nList//Length}]<=Sum[#2[[ii]],{ii,1,nList//Length}])&];
	Return[outout]
]
NonvanishingCombination[nList_,minord_,maxord_]:=Module[{pList={},qList={},p,q,output={},outout={},appendList={},negativeCheck={}}, (* Finding all non-vanishing combinations*)
	pList=Tuples[Range[0,maxord],nList//Length];
	qList=pList;
	Do[
		If[Plus@@pList[[ip]]==0 || Plus@@qList[[iq]]==0,Continue[]];
		If[(Plus@@pList[[ip]])+( Plus@@qList[[iq]])>maxord,Continue[]];
		If[(Plus@@pList[[ip]])+( Plus@@qList[[iq]])<minord,Continue[]];
		If[AnyTrue[pList[[ip]]+qList[[iq]],(#1==0)&],Continue[]];
		If[Plus@@(nList(pList[[ip]]-qList[[iq]]))==0,
			appendList=(pList[[ip]]+qList[[iq]])~Join~Flatten[{Drop[nList(pList[[ip]]-qList[[iq]]),1]}];
			negativeCheck=(pList[[ip]]+qList[[iq]])~Join~(-Flatten[{Drop[nList(pList[[ip]]-qList[[iq]]),1]}]);
			If[((Position[output,negativeCheck])//Length)!=0,Continue[]];
			output=Append[output,appendList];
		];
	,{ip,1,pList//Length},{iq,1,qList//Length}];
	outout=Sort[DeleteDuplicates[output],(Sum[#1[[ii]],{ii,1,nList//Length}]<=Sum[#2[[ii]],{ii,1,nList//Length}])&];
Return[outout]
]
sortingFunc[term_]:=Module[{output,out1,toList,mom},(* This function sorts the output moments, from the higherst rank moment to lower ranks *)
	ClearAttributes[Times,Orderless];
	out1=(term/.{\[LeftAngleBracket]Cos[A_.] B_.\[RightAngleBracket]:> \[LeftAngleBracket]B Cos[A] //Simplify\[RightAngleBracket]});
	ClearAttributes[Plus,Orderless];
	toList=List@@out1;
	If[(toList//Length)==1,
		output=ToString[out1,StandardForm];
	,
		output=ToString[((Plus@@(Sort[toList,(Exponent[(#1/.{\[LeftAngleBracket]xx_\[RightAngleBracket]:>mom}),mom]<=Exponent[(#2/.{\[LeftAngleBracket]xx_\[RightAngleBracket]:>mom}),mom])&]))),StandardForm];
	];
	SetAttributes[Plus,Orderless];
	SetAttributes[Times,Orderless];
	Return[output]
]
averageToCorrlation[mList_,delList_,nlist_,symb_]:=Module[{output={},out1={},extras,secondAndHigherHarmonics={},SeriesRange,firstHarmonic,numberOfHarmonics,finalHarmonics},
	numberOfHarmonics=nlist//Length;
	firstHarmonic=Table[-Sign[Plus@@delList]nlist[[1]],{i1,1,Abs[Plus@@delList]/nlist[[1]]}];
	Do[
		secondAndHigherHarmonics=Append[secondAndHigherHarmonics,Table[Sign[delList[[i1]]]nlist[[i1+1]],{i2,1,Abs[delList[[i1]]]/nlist[[i1+1]]}]];
	,{i1,1,numberOfHarmonics-1}];
	out1={firstHarmonic}~Join~secondAndHigherHarmonics;
	Do[
		finalHarmonics=mList[[i3]]-(out1[[i3]]//Length);
	If[finalHarmonics==0,
		output=Append[output,out1[[i3]]];
	,
		extras=(Table[{-nlist[[i3]],nlist[[i3]]},{i4,1,finalHarmonics/2}]//Flatten);
		out1[[i3]]=out1[[i3]]~Join~extras;
		output=Append[output,out1[[i3]]];
		Clear[extras];
	]
	,{i3,1,mList//Length}];
	SeriesRange=Sort[Flatten[output],#1<=#2&];
	output=symb[Sequence@@SeriesRange//Evaluate];
	Return[output]
]


(* Mathematica Error Messages *)

harmonicOrder::argerr = "Harmonic list should be in strictly ascending order.";
alphaListSize::argerr = "The phase list size should be one unit less than harmonic list size.";
orderListSize::argerr = "The cumulant order list size should be equal to the harmonic list size.";
correlationSum::argerr = "The total sum of the correlation list elements should be zero.";


(* Generating Cumulants  *)

c[MAmplitureList_,m\[Delta]List_,HarmonicList_,ampSymb_,phaseSymb_]:=Module[{output,w1,w2,w3,w4,w5,w6,HighestMoment,SeriesRange,generatingFuc,gExponent,coeff,fazList,\[Delta]kList,\[Delta]List,kList,rList,avZeroList1,avZeroList2,amp,\[CapitalDelta]},
(*   =====================================================================\[Equal]     *)
(*                                Errors                                        *)
	If[!OrderedQ[HarmonicList,(#1<#2)&] && ((HarmonicList)//Length)!=1,
		Message[harmonicOrder::argerr,1];
		Return[Null];
	];
	If[((MAmplitureList//Length)!=(HarmonicList//Length)),
		Message[orderListSize::argerr,1];
		Return[Null];
	];
	If[((m\[Delta]List//Length)+1!=(HarmonicList//Length)),
		Message[alphaListSize::argerr,1];
		Return[Null];
	];
(*   =====================================================================\[Equal]     *)
(*                     Preparing the initiation lists                           *)
	kList=Table[ki[ii],{ii,1,MAmplitureList//Length}];
	\[Delta]kList=Table[\[Delta]ki[ii],{ii,1,(m\[Delta]List//Length)}]; 
	rList=Table[ri[ii],{ii,1,MAmplitureList//Length}];
	\[Delta]List=Table[\[Delta]i[ii],{ii,1,m\[Delta]List//Length}];
	fazList=Table[faz[ii],{ii,1,m\[Delta]List//Length}];
(*   =====================================================================\[Equal]     *)
(*                     Preparing the Generating Function                        *)
	generatingFuc[kList_,\[Delta]kList_]:=Exp[I (rList[[1]] kList[[1]] Cos[HarmonicList[[1]] fi]+Sum[rList[[ii+1]] kList[[ii+1]]  Cos[HarmonicList[[ii+1]](fi+\[Delta]List[[ii]]-\[Delta]kList[[ii]])],{ii,1,m\[Delta]List//Length}])];
	SeriesRange=Table[{kList[[ii]],0,MAmplitureList[[ii]]},{ii,1,MAmplitureList//Length}];
(*   =====================================================================\[Equal]     *)
(*                       Expanding the Integrand                                *)
	w1=((((Series[generatingFuc[kList,\[Delta]kList],Sequence@@SeriesRange//Evaluate])//Normal)//TrigReduce)/.{Cos[\[Alpha]_]:>(E^(I \[Alpha])+E^(-I \[Alpha]))/2}/.{Sin[\[Alpha]_]:>(E^(I \[Alpha])-E^(-I \[Alpha]))/(2I)})//Expand;
(*   =====================================================================\[Equal]     *)
(*         To increase the speed, we use "integrateFun" for integration         *)
	w2=DeleteCases[(*Monitor[*)Table[integrateFunc[w1[[i]]],{i,1,w1//Length}](*,i/(w1//Length)//N]*),0];
	w3=w2/.Table[\[Delta]i[ii]->Log[faz[ii]]/I,{ii,1,5}]//Simplify;
	avZeroList1=Table[0,{ii,1,MAmplitureList//Length}];
	avZeroList2=Table[0,{ii,1,m\[Delta]List//Length}];
	w4=(Plus@@Table[reduceToAverageNice[w3[[i]],rList,fazList,HarmonicList,ampSymb,phaseSymb],{i,1,w3//Length}])/.{\[LeftAngleBracket]1\[RightAngleBracket]->1,\[LeftAngleBracket]0\[RightAngleBracket]->0};
	w5=SeriesCoefficient[Log[w4],Sequence@@SeriesRange//Evaluate];
(*  ==================\[Equal]  *)
	If[(m\[Delta]List//Length)==0,
		w6=w5;
	,
		w6=FourierCoefficient[w5,\[Delta]kList,-m\[Delta]List];
	];
	If[w6==0,
		Return[w6];
		Exit[];
	];
	w6=(w6/.{E^\[Alpha]_:>Cos[(\[Alpha]/I)//Expand]})(*//Simplify//Expand*);
(*   =====================================================================\[Equal]     *)
(*   In this part, we extract the numerical factor behind the highest moment.   *)
	amp=Product[Subscript[ampSymb, HarmonicList[[i]]]^MAmplitureList[[i]],{i,1,HarmonicList//Length}];
	\[CapitalDelta]=Sum[m\[Delta]List[[i]](Subscript[phaseSymb, HarmonicList[[i+1]]]-Subscript[phaseSymb, HarmonicList[[1]]]),{i,1,m\[Delta]List//Length}]//FullSimplify;
	HighestMoment=AngleBracket[Cos[\[CapitalDelta]//Expand] amp ];
	coeff=D[w6,HighestMoment];
(*   =====================================================================\[Equal]     *)
(*                               The output                                     *)
	output=sortingFunc[((w6/coeff))//Simplify//Expand];
	Return[output]
]
cCorr[MAmplitureList_,m\[Delta]List_,HarmonicList_,corrSymb_]:=Module[{output,w1,w2,w3,w4,w5,w6,SeriesRange,generatingFuc,gExponent,coeff,fazList,\[Delta]kList,\[Delta]List,kList,rList,avZeroList1,avZeroList2},
(*   =====================================================================\[Equal]     *)
(*                                Errors                                        *)
	If[!OrderedQ[HarmonicList,(#1<#2)&] && ((HarmonicList)//Length)!=1,
		Message[harmonicOrder::argerr,1];
		Return[Null];
	];
	If[((MAmplitureList//Length)!=(HarmonicList//Length)),
		Message[orderListSize::argerr,1];
		Return[Null];
	];
	If[((m\[Delta]List//Length)+1!=(HarmonicList//Length)),
		Message[alphaListSize::argerr,1];
		Return[Null];
	];
(*   =====================================================================\[Equal]     *)
(*                     Preparing the initiation lists                           *)	
	kList=Table[ki[ii],{ii,1,MAmplitureList//Length}];
	\[Delta]kList=Table[\[Delta]ki[ii],{ii,1,(m\[Delta]List//Length)}]; 
	rList=Table[ri[ii],{ii,1,MAmplitureList//Length}];
	\[Delta]List=Table[\[Delta]i[ii],{ii,1,m\[Delta]List//Length}];
	fazList=Table[faz[ii],{ii,1,m\[Delta]List//Length}];
(*===============*)
	generatingFuc[kList_,\[Delta]kList_]:=Exp[I (rList[[1]] kList[[1]] Cos[HarmonicList[[1]] fi]+Sum[rList[[ii+1]] kList[[ii+1]]  Cos[HarmonicList[[ii+1]](fi+\[Delta]List[[ii]]-\[Delta]kList[[ii]])],{ii,1,m\[Delta]List//Length}])];
	SeriesRange=Table[{kList[[ii]],0,MAmplitureList[[ii]]},{ii,1,MAmplitureList//Length}];
	w1=((((Series[generatingFuc[kList,\[Delta]kList],Sequence@@SeriesRange//Evaluate])//Normal)//TrigReduce)/.{Cos[\[Alpha]_]:>(E^(I \[Alpha])+E^(-I \[Alpha]))/2}/.{Sin[\[Alpha]_]:>(E^(I \[Alpha])-E^(-I \[Alpha]))/(2I)})//Expand;
	w2=DeleteCases[Table[integrateFunc[w1[[i]]],{i,1,w1//Length}],0];
	w3=w2/.Table[\[Delta]i[ii]->Log[faz[ii]]/I,{ii,1,5}]//Simplify;
	avZeroList1=Table[0,{ii,1,MAmplitureList//Length}];
	avZeroList2=Table[0,{ii,1,m\[Delta]List//Length}];
	w4=(Plus@@Table[reduceToAverage[w3[[i]],rList,fazList(*,HarmonicList,v,\[Psi]*)],{i,1,w3//Length}])/.{av[avZeroList1,avZeroList2]->1};
	w5=SeriesCoefficient[Log[w4],Sequence@@SeriesRange//Evaluate];
(*  ==================\[Equal]  *)
	If[(m\[Delta]List//Length)==0,
		w6=w5;
	,
		w6=FourierCoefficient[w5,\[Delta]kList,-m\[Delta]List];
	];
	If[w6==0,
		Return[w6];
		Exit[];
	];
	coeff=D[w6,av[MAmplitureList,m\[Delta]List]];
	output=((w6/coeff)/.{av[a_,b_]:>averageToCorrlation[a,b,HarmonicList,corrSymb]})//Expand;
	Return[output]
]
cTable[nList_,maxord_,ampSymb_,phaseSymb_]:=Module[{output={},out1={},cumul,kkList,\[Delta]kkList,allowedList={},str1,str2,headerString={}},
(*   =====================================================================\[Equal]     *)
(*                                Errors                                        *)
	If[!OrderedQ[nList,(#1<#2)&] && ((nList)//Length)!=1,
		Message[harmonicOrder::argerr,1];
		Return[Null];
	];
(*   =====================================================================\[Equal]     *)
	allowedList=NonvanishingCombination[nList,maxord];
	Monitor[Do[
		kkList=Drop[allowedList[[ii]],-((nList//Length)-1)];
		\[Delta]kkList=Drop[allowedList[[ii]],(nList//Length)];
		cumul=c[kkList,\[Delta]kkList,nList,ampSymb,phaseSymb];
		If[cumul==0 ,Continue[]];
		str1=ToString[kkList];
		str2=ToString[\[Delta]kkList];
		If[(\[Delta]kkList//Length)==0,
			out1=Append[out1,{str1,cumul}];
		,
			out1=Append[out1,{str1,str2,cumul}];
		];
	,{ii,1,allowedList//Length}],ToString[(ii*100)/(allowedList//Length)//N//IntegerPart]~StringJoin~"%"];
	If[(\[Delta]kkList//Length)==0,
		headerString={{"{m}",ToString[Subscript[c, nList[[1]]],StandardForm]~StringJoin~"{m}"}};
	,
		headerString={{ToString[Table[Subscript["m", ToString[i,StandardForm]],{i,1,kkList//Length}],StandardForm],ToString[Table[Subscript["\[Alpha]", ToString[i,StandardForm]],{i,1,(kkList//Length)-1}],StandardForm],ToString[Subscript[c, Row[nList,","]]^ToString[Table[Subscript["\[Alpha]", ToString[i,StandardForm]],{i,1,(kkList//Length)-1}],StandardForm],StandardForm]~StringJoin~ToString[Table[Subscript["m", ToString[i,StandardForm]],{i,1,kkList//Length}],StandardForm]}};
	];
	output=Grid[headerString~Join~out1,Dividers->{{},{1->{Thick,Black},2->{Thick,Black},((allowedList//Length)+2)->{Thick,Black}}}(*,ItemSize\[Rule]10*),Alignment->Left];
	Return[output];
]
cTable[nList_,minord_,maxord_,ampSymb_,phaseSymb_]:=Module[{output={},out1={},cumul,kkList,\[Delta]kkList,allowedList={},str1,str2,headerString={}},
(*   =====================================================================\[Equal]     *)
(*                                Errors                                        *)
	If[!OrderedQ[nList,(#1<#2)&] && ((nList)//Length)!=1,
		Message[harmonicOrder::argerr,1];
		Return[Null];
	];
(*   =====================================================================\[Equal]     *)
	allowedList=NonvanishingCombination[nList,minord,maxord];
	Monitor[Do[
		kkList=Drop[allowedList[[ii]],-((nList//Length)-1)];
		\[Delta]kkList=Drop[allowedList[[ii]],(nList//Length)];
		cumul=c[kkList,\[Delta]kkList,nList,ampSymb,phaseSymb];
		If[cumul==0 ,Continue[]];
		str1=ToString[kkList];
		str2=ToString[\[Delta]kkList];
		If[(\[Delta]kkList//Length)==0,
			out1=Append[out1,{str1,cumul}];
		,
			out1=Append[out1,{str1,str2,cumul}];
		];
	,{ii,1,allowedList//Length}],ToString[(ii*100)/(allowedList//Length)//N//IntegerPart]~StringJoin~"%"];
	If[(\[Delta]kkList//Length)==0,
		headerString={{"{m}",ToString[Subscript[c, nList[[1]]],StandardForm]~StringJoin~"{m}"}};
	,
		headerString={{ToString[Table[Subscript["m", ToString[i,StandardForm]],{i,1,kkList//Length}],StandardForm],ToString[Table[Subscript["\[Alpha]", ToString[i,StandardForm]],{i,1,(kkList//Length)-1}],StandardForm],ToString[Subscript[c, Row[nList,","]]^ToString[Table[Subscript["\[Alpha]", ToString[i,StandardForm]],{i,1,(kkList//Length)-1}],StandardForm],StandardForm]~StringJoin~ToString[Table[Subscript["m", ToString[i,StandardForm]],{i,1,kkList//Length}],StandardForm]}};
	];
	output=Grid[headerString~Join~out1,Dividers->{{},{1->{Thick,Black},2->{Thick,Black},((allowedList//Length)+2)->{Thick,Black}}}(*,ItemSize\[Rule]10*),Alignment->Left];
	Return[output];
]


(* Statistical error for moment/n-partile correlations *)

ZeroOneList[size_,level_]:=Module[{output={},zeroList={},oneList={},i1}, (* This function make all possible list with a given size and number of zeros and one. This will be used to produce the partition list. *)
	zeroList=Table[0,{i1,1,size-level}];
	oneList=Table[1,{i1,1,level}];
	output=Append[output,zeroList~Join~ oneList];
	Return[Flatten[Permutations/@output,1]]
]

splitList[list_,splitingBoolean_]:=Module[{output={},i1}, (* This function split a list consists of other two dimensional list into two or leave untouched based on a set made out of zeros and ones. *)
	Do[
		If[splitingBoolean[[i1]]==1,
			output=Append[output,{{list [[i1]][[1]]},{list [[i1]][[2]]}}];
		,
			output=Append[output,{list[[i1]]}];
		];
	,{i1,1,list//Length}];
	Return[Flatten[output,1]]
]

partitioningI[k_,level_,XX_]:=Module[{output={},splitingBoolean,indexList,n,X,ii1,ii2,i1,i2,i3},  (* It generates all possible partitions needed for statistical uncertainty computations. Using functions splitList and  ZeroOneList.*)
	splitingBoolean=ZeroOneList[k,level];
	X=Table[n[ii2],{ii2,1,2k}];
	indexList=Table[{Drop[X,-(X//Length)/2],(Permutations@Drop[X,(X//Length)/2])[[ii1]]}//Transpose,{ii1,1,((X//Length)/2)!}];
	Do[
	Do[
		output=Append[output,splitList[indexList[[i2]],splitingBoolean[[i1]]]];
	,{i1,1,splitingBoolean//Length}];
	,{i2,1,indexList//Length}];
	output=Sort/@output//DeleteDuplicates;
	Return[output/.{n[i3_]:>XX[[i3]]}]
]

partitioningN[nList_,level_]:=Module[{output={},ilist,nPartition,i1}, (* Convert partitioningI into partitoin of a set of harmonics. *)
	ilist=nList~Join~nList;
	nPartition=partitioningI[nList//Length,level,ilist];
	Do[
		output=Append[output,Plus@@@(nPartition[[i1]])];
	,{i1,1,nPartition//Length}];
	Return[output]
]

fallingFactorial[M_,k_]:=Product[(M-i1),{i1,0,k-1}]

momentToNiceMoment[mList_,delList_,nlist_,ampSymb_,phaseSymb_]:=Module[{output,amp,phase,i1,i2},
	amp=Product[Subscript[ampSymb, nlist[[i1]]]^mList[[i1]],{i1,1,nlist//Length}];
	phase=Sum[-delList[[i2-1]](Subscript[phaseSymb, nlist[[i2]]]-Subscript[phaseSymb, nlist[[1]]]),{i2,2,nlist//Length}]//Simplify;
	output=AngleBracket[amp Cos[phase]];
	Return[output]
]

niceMomentToCorr[moment_,symbCorr_,ampSymb_,phaseSymb_]:=Module[{output,dummy,amp,phase,phase2,amp2,i1,ampList,phaseList,av,harmonicList,mList,deltaList}, (* A function for presenting output in terms of correlators. *)
	amp=(moment/.{\[LeftAngleBracket]A_. Cos[B_.] \[RightAngleBracket]:>av[A,B]}/.{\[LeftAngleBracket]A_  \[RightAngleBracket]:>av[A,0]})[[1]];
	phase=(moment/.{\[LeftAngleBracket]A_. Cos[B_.] \[RightAngleBracket]:>av[A,B]}/.{\[LeftAngleBracket]A_  \[RightAngleBracket]:>av[A,0]})[[2]]//Expand;
	If[amp==1,Return[1]];
	amp2=Drop[List@@(dummy * amp/.{Subscript[ampSymb, n_.]^m_.:>ampList[n,m]}),1]; (* dummy is added for the cases there is only one term in < > (otherwise changeing produt to list becomes crasy) *)
	amp2=Table[amp2[[i1]],{i1,1,amp2//Length}]/.{ampList[nn_.,mm_.]:>{nn,mm}};
	harmonicList=(amp2//Transpose)[[1]];
	mList=(amp2//Transpose)[[2]];
	phase2=Table[Coefficient[phase,Subscript[phaseSymb, harmonicList[[i1]]]],{i1,1,harmonicList//Length}];
	deltaList=-Drop[phase2,1];
	output=averageToCorrlation[mList,deltaList,harmonicList,symbCorr];
	Return[output];
]

momentToCorr[moment_,symbCorr_,ampSymb_,phaseSymb_]:=Module[{output,dummy,amp,phase,phase2,amp2,i1,ampList,phaseList,av,harmonicList,mList,deltaList}, (* A function for presenting output in terms of correlators. *)
	amp=(moment/.{\[LeftAngleBracket]A_. Cos[B_.] \[RightAngleBracket]:>av[A,B]}/.{\[LeftAngleBracket]A_  \[RightAngleBracket]:>av[A,0]})[[1]];
	phase=(moment/.{\[LeftAngleBracket]A_. Cos[B_.] \[RightAngleBracket]:>av[A,B]}/.{\[LeftAngleBracket]A_  \[RightAngleBracket]:>av[A,0]})[[2]]//Expand;
	If[amp==1,Return[1]];
	amp2=Drop[List@@(dummy * amp/.{Subscript[ampSymb, n_.]^m_.:>ampList[n,m]}),1]; (* dummy is added for the cases there is only one term in < > (otherwise changeing produt to list becomes crasy) *)
	amp2=Table[amp2[[i1]],{i1,1,amp2//Length}]/.{ampList[nn_.,mm_.]:>{nn,mm}};
	harmonicList=(amp2//Transpose)[[1]];
	mList=(amp2//Transpose)[[2]];
	phase2=Table[Coefficient[phase,Subscript[phaseSymb, harmonicList[[i1]]]],{i1,1,harmonicList//Length}];
	deltaList=-Drop[phase2,1];
	output=averageToCorrlation[mList,deltaList,harmonicList,symbCorr];
	Return[output];
]

corrToNicemoment[corrList_,AmpSymbol_,phaseSymb_]:=Module[{output,amp,phase,i1},
	amp=Product[Subscript[AmpSymbol, Abs@corrList[[i1]]],{i1,1,(corrList//Length)}];
	phase=Sum [corrList[[i1]]Subscript[phaseSymb, Abs@corrList[[i1]]],{i1,1,(corrList//Length)}];	
	output=AngleBracket[amp Cos[phase]];
	Return[output]
]

Nsigma2Moment[mList_,delList_,nlist_,AmpSymbol_,phaseSymb_,M_]:=Module[{output=0,listOfPartitions,ampList,fazeList,i1,i2,corrList,corrSymb},
(*   =====================================================================\[Equal]     *)
(*                                Errors                                        *)
	If[!OrderedQ[nlist,(#1<#2)&] && ((nlist)//Length)!=1,
		Message[harmonicOrder::argerr,1];
		Return[Null];
	];
	If[((mList//Length)!=(nlist//Length)),
		Message[orderListSize::argerr,1];
		Return[Null];
	];
	If[((delList//Length)+1!=(nlist//Length)),
		Message[alphaListSize::argerr,1];
		Return[Null];
	];
	corrList=List@@averageToCorrlation[mList,delList,nlist,corrSymb];
	Do[
		listOfPartitions=partitioningN[corrList,level];
		Do[
			ampList=Product[Subscript[AmpSymbol, Abs@listOfPartitions[[i2]][[i1]]],{i1,1,(listOfPartitions[[i2]]//Length)}];
			fazeList=Sum [listOfPartitions[[i2]][[i1]]Subscript[phaseSymb, Abs@listOfPartitions[[i2]][[i1]]],{i1,1,(listOfPartitions[[i2]]//Length)}];
			output=output+fallingFactorial[M,listOfPartitions[[i2]]//Length]*AngleBracket[ampList Cos[fazeList]];
		,{i2,1,(listOfPartitions)//Length}];
	,{level,0,corrList//Length}];
	output=output/.{Subscript[AmpSymbol, 0]->1,Subscript[phaseSymb, 0]->0}/.{\[LeftAngleBracket]1\[RightAngleBracket]->1};
	output=(output/(fallingFactorial[M,corrList//Length])^2);
	output=(output-momentToNiceMoment[mList,delList,nlist,AmpSymbol,phaseSymb]^2)//Simplify;
Return[output]
]

Nsigma2Moment[corrList_,AmpSymbol_,phaseSymb_,M_]:=Module[{output=0,listOfPartitions,ampList,fazeList,i1,i2,corrSymb},
(*   =====================================================================\[Equal]     *)
(*                                Errors                                        *)
	If[Plus@@corrList!=0,
		Message[correlationSum::argerr,1];
		Return[Null];
	];
	Do[
		listOfPartitions=partitioningN[corrList,level];
		Do[
			ampList=Product[Subscript[AmpSymbol, Abs@listOfPartitions[[i2]][[i1]]],{i1,1(listOfPartitions[[i2]]//Length)}];
			fazeList=Sum [listOfPartitions[[i2]][[i1]]Subscript[phaseSymb, Abs@listOfPartitions[[i2]][[i1]]],{i1,1(listOfPartitions[[i2]]//Length)}];
			output=output+fallingFactorial[M,listOfPartitions[[i2]]//Length]*AngleBracket[ampList Cos[fazeList]];
		,{i2,1,(listOfPartitions)//Length}];
	,{level,0,corrList//Length}];
	output=output/.{Subscript[AmpSymbol, 0]->1,Subscript[phaseSymb, 0]->0}/.{\[LeftAngleBracket]1\[RightAngleBracket]->1};
	output=(output/(fallingFactorial[M,corrList//Length])^2);
	output=((output-corrToNicemoment[corrList,AmpSymbol,phaseSymb]^2))//Simplify;
Return[output]
]

Nsigma2MomentCorr[corrList_,symbCorr_,M_]:=Module[{output=0,listOfPartitions,ampList,fazeList,AmpSymbol,phaseSymb,i1,i2,nnn,niceMoment},
(*   =====================================================================\[Equal]     *)
(*                                Errors                                        *)
	If[Plus@@corrList!=0,
		Message[correlationSum::argerr,1];
		Return[Null];
	];
	Do[
		listOfPartitions=partitioningN[corrList,level];
		Do[
			ampList=Product[Subscript[AmpSymbol, Abs@listOfPartitions[[i2]][[i1]]],{i1,1(listOfPartitions[[i2]]//Length)}];
			fazeList=Sum [listOfPartitions[[i2]][[i1]]Subscript[phaseSymb, Abs@listOfPartitions[[i2]][[i1]]],{i1,1(listOfPartitions[[i2]]//Length)}];
			niceMoment=((AngleBracket[ampList Cos[fazeList]])/.{Subscript[AmpSymbol, 0]->1});
			output=output+fallingFactorial[M,listOfPartitions[[i2]]//Length]*niceMomentToCorr[niceMoment,symbCorr,AmpSymbol,phaseSymb];
		,{i2,1,(listOfPartitions)//Length}];
	,{level,0,corrList//Length}];
	output=(output/(fallingFactorial[M,corrList//Length])^2)//Simplify;
	output=output-symbCorr[Sequence@@corrList]^2;
Return[output]
]


Protect[
	c,
	cCorr,
	cTable,
	Nsigma2Moment,
	Nsigma2MomentCorr
]


End[ ]
EndPackage[ ]
