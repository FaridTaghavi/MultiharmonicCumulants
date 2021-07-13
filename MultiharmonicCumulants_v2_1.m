(* ::Package:: *)

(*  
MultiharmonicCumulants_v2_1.m

Author: Seyed Farid Taghavi
Date: July 13, 2021

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
	12. Nsigma2P,
	13. cppInterpreter.
	
To run, save the package together with a Mathematica notebook (example.nb for instance) in the same folder. By running the following lines in the notebook, the functions in the package are accessible:

	SetDirectory[NotebookDirectory[]];
	<<MultiharmonicCumulants_v2_1.m

As an example:
	
	in[]:= cCorr[{3, 1, 1, 1}, {4, 5, -6}, {3, 4, 5, 6}, corr]
	out[]= -corr[-3, 3] (corr[-6, -3, 4, 5] + corr[-5, -4, 3, 6]) + 1/2 (corr[-6, -3, -3, 3, 4, 5] + corr[-5, -4, -3, 3, 3, 6])	
	
An explanation for each function can be obtained by adding a "?" at the beginning of the function's name. For instance,

	?cQvec	
	
returns
	
	cQvec[cumulant_order_list , phase_list, harmonic_list , multiplicity_symbole , Q-vector_symbol]: cumulant written in terms of Q-vectors.
*)


BeginPackage["MultiharmonicCumulants`"]


Unprotect[Power];
Format[Power[E,a_],CForm]:=exp[a]
Format[Power[a_,1/2],CForm]:=sqrt[a]
Format[Power[a_,b_],CForm]:=pow[a,b]
Protect[Power];


Unprotect[
	c,
	cCorr,
	cQvec,
	cMean,
	cTable,
	ncQvec,
	ncMean,
	ncCorr,
	Nsigma2,
	Nsigma2Mean,
	Nsigma2Qvec,
	Nsigma2P,
	cppInterpreter
]


ClearAll[
	c,
	cCorr,
	cQvec,
	cMean,
	cTable,
	ncQvec,
	ncMean,
	ncCorr,
	Nsigma2,
	Nsigma2Mean,
	Nsigma2Qvec,
	Nsigma2P,
	cppInterpreter
]


c::usage=     
   "c[ cumulant_order_list, phase_list , harmonic_list , flow_amplitude_symbole , flow_phase_symbol]: cumulant written in terms of symbolic moments "
cCorr::usage=     
   "cCorr[ cumulant_order_list, phase_list , harmonic_list ]: cumulant written in terms of symbolic correlations OR cCorr[ cumulant_order_list, phase_list , harmonic_list , correlation_symbol]: cumulant written in terms of correlations"
cQvec::usage=
	"cQvec[ cumulant_order_list, phase_list , harmonic_list , multiplicity_symbole , Q-vector_symbol]: cumulant written in terms of Q-vectors "
cMean::usage= 
	   "cMean[ cumulant_order_list, phase_list , harmonic_list , flow_amplitude_symbole , flow_phase_symbol]: cumulant written in terms of Mean[...]"
cTable::usage=     
   "cTable[ harmonic_list , maximum_order , flow_amplitude_symbole , flow_phase_symbol] OR cTable[ harmonic_list , minimum_order , maximum_order , flow_amplitude_symbole , flow_phase_symbol]: table of all cumulants written in terms of symbolic moments"
ncQvec::usage=
	"ncQvec[ cumulant_order_list, phase_list , harmonic_list , multiplicity_symbole , Q-vector_symbol]: normalized version of cQvec"
ncMean::usage=
	"cMean[ cumulant_order_list, phase_list , harmonic_list , flow_amplitude_symbole , flow_phase_symbol]: normalized version of cMean"
ncCorr::usage=   
	"ncCorr[ cumulant_order_list, phase_list , harmonic_list , correlation_symbol]: normalized version of ncCorr"
Nsigma2::usage=
	"Nsigma2[ any function of corr[a1,...,an] ]: the statistica error of any function of correlations written in terms of symbolic correlators."
Nsigma2Mean::usage=
	"Nsigma2[ any function of corr[a1,...,an] , flow_amplitude_symbole , flow_phase_symbol]: similar to Nsigma2 witten in terms of Mean[...]"
Nsigma2Qvec::usage=
	"Nsigma2[ any function of corr[a1,...,an] , multiplicity_symbole , Q-vector_symbol]: similar to Nsigma2 witten in terms Q-vectors"
Nsigma2P::usage=
	"Nsigma2P[ any function of corr[a1,...,an] ]: similar to Nsigma2 with a direct algorithm explained in the appendix of Ref. [arXiv: 2005.04742]"
cppInterpreter::usage=
	"cppInterpreter[ function of \[LeftAngleBracket]k\!\(\*SubscriptBox[\(\[RightAngleBracket]\), \(\(-k\),  ... , k\)]\)\[RightAngleBracket] , path and name of the output file , correlation_symbol]: translates any function of a correlators to Cpp format"


Begin["`Private`"]


(* Useful Functions For Generating Cumulant*)

integrateFunc[term_]:=Module[{output,faze,Dfaze,fazeAlpha,fazeBeta,coeff}, 
	(* ====================================================================
	We use a simple substitution instead of integrating to increas the computation speed. 
	==============================================================================*)
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
NonvanishingCombination[nList_,minord_,maxord_]:=Module[{pList={},qList={},p,q,output={},outout={},appendList={},negativeCheck={}}, 
	(* =======================================
	Finding all non-vanishing combinations
	from order min to order max.
	======================================*)
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
sortingFunc[term_]:=Module[{output,out1,toList,mom},
	(* ==========================================================
	This function sorts the output moments, from the higherst rank 
	moment to lower ranks. 
	============================================================*)
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
averageToNice[rExpList_,fazExpList_,HarmonicList_,AmpSymbol_,phaseSymbol_]:=Module[{output,amp,\[CapitalDelta]},
	(*The following is just cosmetic*)
	amp=Product[Subscript[AmpSymbol, HarmonicList[[i]]]^rExpList[[i]],{i,1,HarmonicList//Length}];
	\[CapitalDelta]=Sum[fazExpList[[i]](Subscript[phaseSymbol, HarmonicList[[i+1]]]-Subscript[phaseSymbol, HarmonicList[[1]]]),{i,1,fazExpList//Length}]//Simplify;
	(*The above is just cosmetic*)
	output= AngleBracket[Cos[\[CapitalDelta]//Expand] amp ]//Simplify;
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
	(*output=symb[Sequence@@SeriesRange//Evaluate];*)
	(* We only keep the real part of the averages. However, correlators are are complex functions. So we need to keep only the real part*)
	output=(symb[Sequence@@SeriesRange//Evaluate]+symb[Sequence@@(Sort@(-SeriesRange))//Evaluate])/2;
	Return[output]
]

averageToNiceCorr[mList_,delList_,nlist_]:=Module[{output={},out1={},extras,secondAndHigherHarmonics={},SeriesRange,firstHarmonic,numberOfHarmonics,finalHarmonics},
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
	(*output=symb[Sequence@@SeriesRange//Evaluate];*)
	(* We only keep the real part of the averages. However, correlators are are complex functions. So we need to keep only the real part*)
	output=(AngleBracket[Subscript[AngleBracket[SeriesRange//Length],Sequence@@SeriesRange//Evaluate]]+AngleBracket[Subscript[AngleBracket[SeriesRange//Length],Sequence@@(Sort@(-SeriesRange))//Evaluate]])/2;
	Return[output]
]

(* correlations in terms of Q-vectors  *)
corrQ[nList2_,M_,Q_]:=Module[{output2,num,denom,NK,NKP,n,replaceList1,replaceList2}, 
	(* ===================================================================================
	 This function is written to find the correlations in terms of Q-vectors based on the algorithm explained in paper [arXiv: 1312.3572]. 
	 =========================================================================================*)
	NK[nListP_]:=Module[{output=0,nList,clist,reducedNList,k,q,m,nn},
		nList=Table[nn[i],{i,1,nListP//Length}];
		m=nList//Length;
		reducedNList=Drop[nList,-1];
		clist=Subsets[reducedNList];
		Do[
			k=clist[[cInd]]//Length;
			q=Plus@@Complement[nList,clist[[cInd]]];
			output+=(-1)^(m-k-1) (m-k-1)!NKP[clist[[cInd]]] Q[{q},m-k];
		,{cInd,1,clist//Length}];
		nn[i_]:=nListP[[i]]; 
		(* =======================================================================================
		This numerator function does not work properly if we set the list as we want. So here we first start
		 with nn[1],nn[2],... list and finally substitute the actual nList in the result
		 ======================================================================================= *)
		Return[output];
	];
	NKP[l_]:=Piecewise[{{1,(l//Length)==0},{Q[l,1],(l//Length)==1},{NK[l],True}}];
	num=NK[nList2];
	denom=(NK[Table[0,{i,1,nList2//Length}]])(*/.{Q[m_,n_]\[RuleDelayed]M}*);
	replaceList1={Q[n_,m_]:>Q[n]};replaceList2={Q[n_]/;((n//Length)==1 && n[[1]]==0):>M};
	output2=(num/denom)/.replaceList1/.replaceList2/.{Q[l_]:>Q[l[[1]]]}; (* Since we ignored the weight Subscript[w, k] the second index in Q[n,m] is not needed and also Q[{0}]=M. We can simply consider weights later if it is needed. *)
	Return[output2//Simplify]
]
averageToMean[rExpList_,fazExpList_,HarmonicList_,AmpSymbol_,phaseSymbol_]:=Module[{output,amp,\[CapitalDelta]},
	(*The following is just cosmetic*)
	amp=Product[AmpSymbol[HarmonicList[[i]]]^rExpList[[i]],{i,1,HarmonicList//Length}];
	\[CapitalDelta]=Sum[fazExpList[[i]](phaseSymbol[HarmonicList[[i+1]]]-phaseSymbol[HarmonicList[[1]]]),{i,1,fazExpList//Length}]//Simplify;
	(*The above is just cosmetic*)
	output=Mean[amp Cos[\[CapitalDelta]]];
Return[output]
]
timesFunc[list_]:=Module[{output={},elements},
Do[
	elements=list[[ii]];
	If[ListQ[elements],elements=Times@@(elements)];
	output=Append[output,elements];
,{ii,1,list//Length}];
Return[output]
]


(* Mathematica Error Messages *)

harmonicOrder::argerr = "Harmonic list should be in strictly ascending order.";
alphaListSize::argerr = "The phase list size should be one unit less than harmonic list size.";
orderListSize::argerr = "The cumulant order list size should be equal to the harmonic list size.";
correlationSum::argerr = "The total sum of the correlation list elements should be zero.";


(* Generating Cumulants  *)

cRaw[MAmplitureList_,m\[Delta]List_,HarmonicList_,avSymb_]:=Module[{output,w1,w2,w3,w4,w5,w6,SeriesRange,generatingFuc,gExponent,coeff,fazList,\[Delta]kList,\[Delta]List,kList,rList,avZeroList1,avZeroList2},
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
	output=((w6/coeff)/.{av[a__]:>avSymb[a]})//Expand;
	Return[output]
]
c[MAmplitureList_,m\[Delta]List_,HarmonicList_,AmpSymbol_,phaseSymbol_]:=(cRaw[MAmplitureList,m\[Delta]List,HarmonicList,av]/.{av[a_,b_]:>averageToNice[a,b,HarmonicList,AmpSymbol,phaseSymbol]})//sortingFunc
cCorr[MAmplitureList_,m\[Delta]List_,HarmonicList_,corrSymb_]:=cRaw[MAmplitureList,m\[Delta]List,HarmonicList,av]/.{av[a_,b_]:>averageToCorrlation[a,b,HarmonicList,corrSymb]}
cCorr[MAmplitureList_,m\[Delta]List_,HarmonicList_]:=cRaw[MAmplitureList,m\[Delta]List,HarmonicList,av]/.{av[a_,b_]:>averageToNiceCorr[a,b,HarmonicList]}

cQvec[MAmplitureList_,m\[Delta]List_,HarmonicList_,MultiSymb_,QvecSymb_]:=Module[{output,mSymb,QSymb,corrSymb},
	(* =============================================================================================================
	Similar to cMean function, we could do a simple replacement of corr[a] with CorrQ function. But because CorrQ function does some symbolic computations,
    the reulsts would mess up when the input "MultiSymb" and "QvecSymb" are numerical lists. So in the module, we first compute CorrQ and then we substitute the final
	"MultiSymb" and "QvecSymb". Then the function works well for both symbolic inputs and numerical lists.
	================================================================================================== *)
	output=cCorr[MAmplitureList,m\[Delta]List,HarmonicList,corrSymb]/.{corrSymb[a__]:>Mean[corrQ[{a},mSymb,QSymb]]};
	Return[(output/.{mSymb->MultiSymb,QSymb->QvecSymb})]
]

cMean[MAmplitureList_,m\[Delta]List_,HarmonicList_,ampSymb_,phaseSymb_]:=cRaw[MAmplitureList,m\[Delta]List,HarmonicList,av]/.{av[a_,b_]:>averageToMean[a,b,HarmonicList,ampSymb,phaseSymb]}

ncQvec[MAmplitureList_,m\[Delta]List_,HarmonicList_,MultiSymb_,QvecSymb_]:=(cQvec[MAmplitureList,m\[Delta]List,HarmonicList,MultiSymb,QvecSymb]/Product[cQvec[{2},{},{HarmonicList[[i]]},MultiSymb,QvecSymb]^(MAmplitureList[[i]]/2),{i,1,HarmonicList//Length}])
ncMean[MAmplitureList_,m\[Delta]List_,HarmonicList_,ampSymb_,phaseSymb_]:=(cMean[MAmplitureList,m\[Delta]List,HarmonicList,ampSymb,phaseSymb]/Product[cMean[{2},{},{HarmonicList[[i]]},ampSymb,phaseSymb]^(MAmplitureList[[i]]/2),{i,1,HarmonicList//Length}])
ncCorr[MAmplitureList_,m\[Delta]List_,HarmonicList_,corr_]:=(cCorr[MAmplitureList,m\[Delta]List,HarmonicList,corr]/Product[cCorr[{2},{},{HarmonicList[[i]]},corr]^(MAmplitureList[[i]]/2),{i,1,HarmonicList//Length}])
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
cTable[nList_,maxord_,ampSymb_,phaseSymb_]:=cTable[nList,1,maxord,ampSymb,phaseSymb]


corrToMean[corrList_,AmpSymbol_,phaseSymbol_]:=Module[{output,amp,\[CapitalDelta]},
	(*The following is just cosmetic*)
	amp=Product[AmpSymbol[Abs[corrList[[i]]]],{i,1,corrList//Length}]//Simplify;
	\[CapitalDelta]=Sum[corrList[[i]](phaseSymbol[Abs[corrList[[i]]]]),{i,1,corrList//Length}]//Simplify;
	(*The above is just cosmetic*)
	output=amp Cos[\[CapitalDelta]];
Return[output]
]

Nsigma2[func_]:=Module[{output=0,variables,variableList,corr},
	variables=Reduce`FreeVariables[func];
	variableList=List@@@variables;
	(* ==== Below: Extracting the symbol of the correlation from the input ========= *)
	corr=ToExpression[StringDelete[ToString[variables[[1]]],"["<>StringDelete[StringDelete[ToString[variableList[[1]]],"{"],"}"]<>"]"]];
	(*  ===============\[Equal]  *)
		output=Sum[
		((D[func,variables[[i1]]]*D[func,variables[[i2]]])/.{corr[a__]:>AngleBracket[Subscript[AngleBracket[{a}//Length],a]]})*
		(AngleBracket[Subscript[AngleBracket[variableList[[i1]]//Length],Sequence@@variableList[[i1]]]*Subscript[AngleBracket[variableList[[i2]]//Length],Sequence@@variableList[[i2]]]]
		-AngleBracket[Subscript[AngleBracket[variableList[[i1]]//Length],Sequence@@variableList[[i1]]]]*AngleBracket[Subscript[AngleBracket[variableList[[i2]]//Length],Sequence@@variableList[[i2]]]])
	,{i1,1,variables//Length},{i2,1,variables//Length}];
	Return[output]
]

Nsigma2[func_,correl_]:=Module[{output=0,variables,variableList,corr},
	variables=Reduce`FreeVariables[func];
	variableList=List@@@variables;
	(* ==== Below: Extracting the symbol of the correlation from the input ========= *)
	corr=ToExpression[StringDelete[ToString[variables[[1]]],"["<>StringDelete[StringDelete[ToString[variableList[[1]]],"{"],"}"]<>"]"]];
	(*  ===============\[Equal]  *)
		output=Sum[
		((D[func,variables[[i1]]]*D[func,variables[[i2]]])/.{corr[a__]:>correl[a]})(correl[Sequence@@variableList[[i1]]][Sequence@@variableList[[i2]]]
		-correl[Sequence@@variableList[[i1]]]*correl[Sequence@@variableList[[i2]]])
	,{i1,1,variables//Length},{i2,1,variables//Length}];
	Return[output]
]

Nsigma2Mean[func_,vSymb_,\[Psi]Symb_]:=Module[{output=0,variables,variableList,corr},
	variables=Reduce`FreeVariables[func];
	variableList=List@@@variables;
	(* ==== Below: Extracting the symbol of the correlation from the input ========= *)
	corr=ToExpression[StringDelete[ToString[variables[[1]]],"["<>StringDelete[StringDelete[ToString[variableList[[1]]],"{"],"}"]<>"]"]];
	(*  ===============\[Equal]  *)
		output=Sum[
		((D[func,variables[[i1]]]*D[func,variables[[i2]]])/.{corr[a__]:>Mean[corrToMean[{a},vSymb,\[Psi]Symb]]})*
		(Mean[corrToMean[variableList[[i1]],vSymb,\[Psi]Symb]*corrToMean[variableList[[i2]],vSymb,\[Psi]Symb]]-Mean[corrToMean[variableList[[i1]],vSymb,\[Psi]Symb]]*Mean[corrToMean[variableList[[i2]],vSymb,\[Psi]Symb]])
	,{i1,1,variables//Length},{i2,1,variables//Length}];
	Return[output]
]

Nsigma2Qvec[func_,MP_,QPSymb_]:=Module[{output=0,M,QSymb,variables,variableList,corr},
(*   =====================================================================\[Equal]     *)
	(* ==== Below: Extracting the symbol of the correlation from the input ========= *)
	variables=Reduce`FreeVariables[func];
	variableList=List@@@variables;
	corr=ToExpression[StringDelete[ToString[variables[[1]]],"["<>StringDelete[StringDelete[ToString[variableList[[1]]],"{"],"}"]<>"]"]];
	(*  ===============\[Equal]  *)
	output=Sum[
		((D[func,variables[[i1]]]*D[func,variables[[i2]]])/.{corr[a__]:>Mean[corr[a]]})(Mean[variables[[i1]]*variables[[i2]]]-Mean[variables[[i1]]]*Mean[variables[[i2]]])
	,{i1,1,variables//Length},{i2,1,variables//Length}];
	(* ==================   *)
	output=(output/.{corr[a__]:>corrQ[{a},M,QSymb]});
Return[output/.{M->MP,QSymb->QPSymb}]
]


(*   =============   *)


(* Direct approach to statistical error for moment/n-partile correlations *)

ZeroOneList[size_,level_]:=Module[{output={},zeroList={},oneList={},i1}, 
	(* ================================================================
	This function make all possible list with a given size and number 
	of zeros and one. This will be used to produce the partition list.
	================================================================== *)
	zeroList=Table[0,{i1,1,size-level}];
	oneList=Table[1,{i1,1,level}];
	output=Append[output,zeroList~Join~ oneList];
	Return[Flatten[Permutations/@output,1]]
]

splitList[list_,splitingBoolean_]:=Module[{output={},i1}, 
	(* ==============================================================
	This function splits a list consists of other two dimensional list into two or 
	leave untouched based on a set made out of zeros and ones. 
	===============================================================*)
	Do[
		If[splitingBoolean[[i1]]==1,
			output=Append[output,{{list [[i1]][[1]]},{list [[i1]][[2]]}}];
		,
			output=Append[output,{list[[i1]]}];
		];
	,{i1,1,list//Length}];
	Return[Flatten[output,1]]
]

partitioningI[kk_,ll_,level_,XX_]:=Module[{output={},splitingBoolean,extendingZeroList,indexList,n,X1,X2,ReducedPermX2,k,l,ii1,ii2,i1,i2,i3},  
	(* ======================================================================
	It generates all possible partitions needed for statistical uncertainty computations 
	using functions splitList and  ZeroOneList.
	=========================================================================*)
If[kk<=ll,
k=kk;l=ll;
,
l=kk;k=ll;
];
extendingZeroList=Table[ConstantArray[0,l-k],{ii,1,(ZeroOneList[k,level])//Length}]; (*Since k neq to l, we shoud extend ZeroOneList with l-k extra zeros. Because we do not need to split one-elemets lists appended at the end of the indexlist *)
splitingBoolean=MapThread[Join,{ZeroOneList[k,level],extendingZeroList}]; (* The above zeroList is added to ZeroOneList*)
     X1=Table[n[ii2],{ii2,1,k}];
X2=Table[n[ii2],{ii2,k+1,k+l}];
ReducedPermX2=Permutations@X2;
indexList=Table[({X1,Drop[ReducedPermX2[[ii1]],-(l-k)]}//Transpose)~Join~Partition[Drop[ReducedPermX2[[ii1]],k],1],{ii1,1,l!}];
 Do[
	Do[
		output=Append[output,splitList[indexList[[i2]],splitingBoolean[[i1]]]];
		,{i1,1,splitingBoolean//Length}];
	,{i2,1,indexList//Length}];
	output=Sort/@output//DeleteDuplicates;
      Return[output/.{n[i3_]:>XX[[i3]]}]
]

partitioningN[lis1_,lis2_,level_]:=Module[{output={},nList1={},nList2={},ilist,nPartition,i1}, 
	(* =======================================\[Equal]======  
	Convert partitioningI into partitoin of a set of harmonics.
	================================================ *)
If[(lis1//Length)>=(lis2//Length),
nList1=lis1;nList2=lis2;
,
nList2=lis1;nList1=lis2;
];
	ilist=nList1~Join~nList2;
	nPartition=partitioningI[nList1//Length,nList2//Length,level,ilist];
     Do[
           output=Append[output,Plus@@@(nPartition[[i1]])];
	,{i1,1,nPartition//Length}];
	Return[output]
]

fallingFactorial[M_,k_]:=Product[(M-i1),{i1,0,k-1}]


NcovarianceCorr[corrList1_,corrList2_,M_]:=Module[{output=0,listOfPartitions,i1,i2,argList},
(*   =====================================================================\[Equal]     *)
(*                                Errors                                        *)
	If[Plus@@corrList1!=0 || Plus@@corrList2!=0 ,
		Message[correlationSum::argerr,1];
		Return[Null];
	];
	Do[
		listOfPartitions=partitioningN[corrList1,corrList2,level];
		Do[
			argList=DeleteCases[Sort@listOfPartitions[[i2]],0];
			output=output+(fallingFactorial[M,listOfPartitions[[i2]]//Length]/((fallingFactorial[M,corrList1//Length])(fallingFactorial[M,corrList2//Length])))*
			Subscript[AngleBracket[argList//Length],Sequence@@argList];
		,{i2,1,(listOfPartitions)//Length}];
	,{level,0,Min[corrList1//Length,corrList2//Length]}];
	output=AngleBracket[output//Simplify];
	output=output/.{Subscript[AngleBracket[0]]->1};
	output=((output-AngleBracket[Subscript[AngleBracket[corrList1//Length],Sequence@@corrList1]]*AngleBracket[Subscript[AngleBracket[corrList1//Length],Sequence@@corrList1]]))//Simplify;
Return[output//Simplify]
]

Nsigma2P[func_,M_]:=Module[{output=0,variables,variableList,corr},
(*   =====================================================================\[Equal]     *)
	variables=Reduce`FreeVariables[func];
	variableList=List@@@variables;
	(* ==== Below: Extracting the symbol of the correlation from the input ========= *)
	corr=ToExpression[StringDelete[ToString[variables[[1]]],"["<>StringDelete[StringDelete[ToString[variableList[[1]]],"{"],"}"]<>"]"]];
	(*  ===============\[Equal]  *)
	output=Sum[
		((D[func,variables[[i1]]]*D[func,variables[[i2]]])/.{corr[a__]:>AngleBracket[Subscript[AngleBracket[a//Length],a]]})
		*NcovarianceCorr[variableList[[i1]],variableList[[i2]],M]
	,{i1,1,variables//Length},{i2,1,variables//Length}];
Return[output]
]

cppInterpreter[funcP_,path_,symb_]:=Module[{output,func,outputII={},internalList,internalListII,replaceList,strm,corr=symb},
     func=funcP//Expand;
	internalList=((func/.{Subscript[\[LeftAngleBracket]n_\[RightAngleBracket], m__]:>corr[m]}/.{Plus->List}/.{Times->List})//Flatten)/.{p_. \[LeftAngleBracket]n_\[RightAngleBracket]^m_.:>n};  
	If[!ListQ[internalList],internalList={internalList}];
	internalList=DeleteCases[DeleteDuplicates[timesFunc[internalList]],aa_/;NumericQ[aa]];
	replaceList=Table[\[LeftAngleBracket]internalList[[ii]]\[RightAngleBracket]->symb[[ii]],{ii,1,internalList//Length}]//Quiet;
	output=(func/.{Subscript[\[LeftAngleBracket]n_\[RightAngleBracket], m__]:>corr[m]}/.replaceList);
	output=((output/.(f:symb)[i_]:>Part[f,i])//CForm//Quiet);
	(*===*)
	internalListII={Table[ToString[symb[ii]]<>"    "<>"=>"<>"    ",{ii,1,internalList//Length}],internalList}//Transpose;
	(*===*)
	outputII={{"=============  Correlation Dictionary    ================",""},{"",""}}~Join~internalListII~Join~{{"",""},{"=============  Correlation in c++ Form     ================",""},{"",""},{output,""}};
	Export[path,outputII]
]


Protect[
	c,
	cCorr,
	cQvec,
	cMean,
	cTable,
	ncQvec,
	ncMean,
	ncCorr,
	Nsigma2,
	Nsigma2Mean,
	Nsigma2Qvec,
	Nsigma2P,
	cppInterpreter
]


End[ ]
EndPackage[ ]
