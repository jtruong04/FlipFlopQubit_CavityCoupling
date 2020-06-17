(* ::Package:: *)

(* ::Input:: *)
(*(* Non-degenerate perturbation theory*)
(*   -- Input Energy levels E in list of size N, and V in array of size NxN*)
(*    -- Output either energies (N list) or states (NxN array, columns are eigenstates)*)
(**)*)


(* ::Input::Initialization:: *)
\[Sigma][input__]:=
Module[{x,y},
x = List[input];
y =PauliMatrix[x[[1]]];
For[i=2,i<= Length[x],i++,
y=KroneckerProduct[y,PauliMatrix[x[[i]]]];
];
y
];


(* ::Input::Initialization:: *)
Energies[EN_,V_,ORDER_:1]:=Module[{\[Epsilon],N},
(*   EN is list of unperturbed energies
     V is perturbation Hamiltonian matrix
 ORDER is the desired perturbative order *)
\[Epsilon]=EN;
N=Length[EN];
If[ORDER>=1,\[Epsilon]+=Table[V[[n,n]],{n,N}]];
If[ORDER>=2,\[Epsilon]+=Table[Sum[If[PossibleZeroQ[V[[n,k]]*V[[k,n]]],0,(V[[n,k]]*V[[k,n]])/(EN[[n]]-EN[[k]])],{k,Complement[Range[N],{n}]}],{n,N}];
If[ORDER>=3,\[Epsilon]+=Table[
Sum[If[PossibleZeroQ[V[[n,Subscript[k, 3]]]*V[[Subscript[k, 3],Subscript[k, 2]]]*V[[Subscript[k, 2],n]]],0,(V[[n,Subscript[k, 3]]]*V[[Subscript[k, 3],Subscript[k, 2]]]*V[[Subscript[k, 2],n]])/((EN[[n]]-EN[[Subscript[k, 2]]])*(EN[[n]]-EN[[Subscript[k, 3]]]))],{Subscript[k, 2],Complement[Range[N],{n}]},{Subscript[k, 3],Complement[Range[N],{n}]}]
-Sum[If[PossibleZeroQ[V[[n,n]]*V[[n,Subscript[k, 3]]]*V[[Subscript[k, 3],n]]],0,(V[[n,n]]*V[[n,Subscript[k, 3]]]*V[[Subscript[k, 3],n]])/(EN[[n]]-EN[[Subscript[k, 3]]])^2],{Subscript[k, 3],Complement[Range[N],{n}]}],{n,N}]];
];
\[Epsilon]
];


(* ::Input::Initialization:: *)
States[EN_,V_,ORDER_:1]:=Module[{U,N},
N=Length[EN];
U=IdentityMatrix[N];
If[ORDER>=1,U+=Table[If[Subscript[k, 1]!=n,If[PossibleZeroQ[V[[Subscript[k, 1],n]]],0,V[[Subscript[k, 1],n]]/(EN[[n]]-EN[[Subscript[k, 1]]])],0],{Subscript[k, 1],N},{n,N}]];
If[ORDER>=2,U+=Table[If[Subscript[k, 1]!=n,Sum[If[PossibleZeroQ[V[
[Subscript[k, 1],Subscript[k, 2]]]*V[[Subscript[k, 2],n]]],0,(V[
[Subscript[k, 1],Subscript[k, 2]]]*V[[Subscript[k, 2],n]])/((EN[[n]]-EN[[Subscript[k, 1]]])*(EN[[n]]-EN[[Subscript[k, 2]]]))],{Subscript[k, 2],Complement[Range[N],{n}]}]-If[PossibleZeroQ[V[[n,n]]*V[[Subscript[k, 1],n]]],0,(V[[n,n]]*V[[Subscript[k, 1],n]])/(EN[[n]]-EN[[Subscript[k, 1]]])^2],0],{Subscript[k, 1],N},{n,N}]+DiagonalMatrix[Table[Sum[-(1/2)*If[PossibleZeroQ[V[[n,Subscript[k, 1]]]*V[[Subscript[k, 1],n]]],0,(V[[n,Subscript[k, 1]]]*V[[Subscript[k, 1],n]])/(EN[[Subscript[k, 1]]]-EN[[n]])^2],{Subscript[k, 1],Complement[Range[N],{n}]}],{n,N}]]];
If[ORDER>=3,U+=Table[If[Subscript[k, 1]!=n,Sum[-Sum[If[PossibleZeroQ[V[[Subscript[k, 1],Subscript[k, 2]]]*V[[Subscript[k, 2],Subscript[k, 3]]]*V[[Subscript[k, 3],n]]],0,(V[[Subscript[k, 1],Subscript[k, 2]]]*V[[Subscript[k, 2],Subscript[k, 3]]]*V[[Subscript[k, 3],n]])/((EN[[Subscript[k, 1]]]-EN[[n]])*(EN[[n]]-EN[[Subscript[k, 2]]])*(EN[[n]]-EN[[Subscript[k, 3]]]))],{Subscript[k, 3],Complement[Range[N],{n}]}]If[PossibleZeroQ[V[[n,Subscript[k, 2]]]*V[[Subscript[k, 2],n]]*V[[Subscript[k, 1],n]]],0,(V[[n,Subscript[k, 2]]]*V[[Subscript[k, 2],n]]*V[[Subscript[k, 1],n]])/((EN[[Subscript[k, 1]]]-EN[[n]])*(EN[[n]]-EN[[Subscript[k, 2]]]))*(1/(EN[[n]]-EN[[Subscript[k, 1]]])+1/(2(EN[[n]]-EN[[Subscript[k, 2]]])))]+If[PossibleZeroQ[V[[n,n]]*V[[Subscript[k, 1],Subscript[k, 2]]]*V[[Subscript[k, 2],n]]],0,(V[[n,n]]*V[[Subscript[k, 1],Subscript[k, 2]]]*V[[Subscript[k, 2],n]])/((EN[[Subscript[k, 1]]]-EN[[n]])*(EN[[n]]-EN[[Subscript[k, 2]]]))*(1/(EN[[n]]-EN[[Subscript[k, 1]]])+1/(EN[[n]]-EN[[Subscript[k, 2]]]))],{Subscript[k, 2],Complement[Range[N],{n}]}]-If[PossibleZeroQ[V[[n,n]]^2*V[[Subscript[k, 1],n]]],0,(V[[n,n]]^2*V[[Subscript[k, 1],n]])/(EN[[Subscript[k, 1]]]-EN[[n]])^3],0],{Subscript[k, 1],N},{n,N}]+DiagonalMatrix[Table[Sum[-Sum[If[PossibleZeroQ[V[[n,Subscript[k, 2]]]*V[[Subscript[k, 2],Subscript[k, 1]]]*V[[Subscript[k, 1],n]]+V[[Subscript[k, 2],n]]*V[[Subscript[k, 1],Subscript[k, 2]]]*V[[n,Subscript[k, 1]]]],0,(V[[n,Subscript[k, 2]]]*V[[Subscript[k, 2],Subscript[k, 1]]]*V[[Subscript[k, 1],n]]+V[[Subscript[k, 2],n]]*V[[Subscript[k, 1],Subscript[k, 2]]]*V[[n,Subscript[k, 1]]])/(2(EN[[n]]-EN[[Subscript[k, 2]]])^2*(EN[[n]]-EN[[Subscript[k, 1]]]))],{Subscript[k, 2],Complement[Range[N],{n}]}]+If[PossibleZeroQ[V[[n,Subscript[k, 1]]]*V[[Subscript[k, 1],n]]*V[[n,n]]],0,(V[[n,Subscript[k, 1]]]*V[[Subscript[k, 1],n]]*V[[n,n]])/(EN[[n]]-EN[[Subscript[k, 1]]])^3],{Subscript[k, 1],Complement[Range[N],{n}]
}],{n,N}]]];
U
];


(* ::Input::Initialization:: *)
PauliCoefficients[M_]:=(
(*Returns the coefficients for following terms:(1	Subscript[\[Sigma], s,x]	Subscript[\[Sigma], s,y]	Subscript[\[Sigma], s,z]
Subscript[\[Sigma], c,x]	Subscript[\[Sigma], c,x].Subscript[\[Sigma], s,x]	Subscript[\[Sigma], c,x].Subscript[\[Sigma], s,y]	Subscript[\[Sigma], c,x].Subscript[\[Sigma], s,z]
Subscript[\[Sigma], c,y]	Subscript[\[Sigma], c,y].Subscript[\[Sigma], s,x]	Subscript[\[Sigma], c,y].Subscript[\[Sigma], s,y]	Subscript[\[Sigma], c,y].Subscript[\[Sigma], s,z]
Subscript[\[Sigma], c,z]	Subscript[\[Sigma], c,z].Subscript[\[Sigma], s,x]	Subscript[\[Sigma], c,z].Subscript[\[Sigma], s,y]	Subscript[\[Sigma], c,z].Subscript[\[Sigma], s,z]

)*)
Return[({
 {1/4 (M[[1,1]]+M[[2,2]]+M[[3,3]]+M[[4,4]]), 1/4 (M[[1,2]]+M[[2,1]]+M[[3,4]]+M[[4,3]]), 1/4 I (M[[1,2]]-M[[2,1]]+M[[3,4]]-M[[4,3]]), 1/4 (M[[1,1]]-M[[2,2]]+M[[3,3]]-M[[4,4]])},
 {1/4 (M[[1,3]]+M[[2,4]]+M[[3,1]]+M[[4,2]]), 1/4 (M[[1,4]]+M[[2,3]]+M[[3,2]]+M[[4,1]]), 1/4 I (M[[1,4]]-M[[2,3]]+M[[3,2]]-M[[4,1]]), 1/4 (M[[1,3]]-M[[2,4]]+M[[3,1]]-M[[4,2]])},
 {1/4 I (M[[1,3]]+M[[2,4]]-M[[3,1]]-M[[4,2]]), 1/4 I (M[[1,4]]+M[[2,3]]-M[[3,2]]-M[[4,1]]), 1/4 (-M[[1,4]]+M[[2,3]]+M[[3,2]]-M[[4,1]]), 1/4 I (M[[1,3]]-M[[2,4]]-M[[3,1]]+M[[4,2]])},
 {1/4 (M[[1,1]]+M[[2,2]]-M[[3,3]]-M[[4,4]]), 1/4 (M[[1,2]]+M[[2,1]]-M[[3,4]]-M[[4,3]]), 1/4 I (M[[1,2]]-M[[2,1]]-M[[3,4]]+M[[4,3]]), 1/4 (M[[1,1]]-M[[2,2]]-M[[3,3]]+M[[4,4]])}
})];
)


(* ::Input:: *)
(**)


(* ::Input::Initialization:: *)
SWTransH[EN_,V_,M_,ORDER_:1]:=Module[{L,HTrans},
(* Requires bases to be order such that desired submatrix is in upper left (i.e. desired bases listed first) 
Arguments: EN is the list of unperturbed energies;
            V is the perturbation Hamiltonian;
			M is the dimension of the desired Hamiltonian;
			ORDER is the order*)
L=Range[M+1,Length[EN]];
HTrans=Table[0,{M},{M}];
(*HTrans=DiagonalMatrix[EN[[1;;M]]];*)
If[ORDER>=1,HTrans+=V[[1;;M,1;;M]]];
If[ORDER>=2,HTrans+=Table[1/2*Sum[V[[m,l]]*V[[l,m']]*(1/(EN[[m]]-EN[[l]])+1/(EN[[m']]-EN[[l]])),{l,L}],{m,M},{m',M}]];
If[ORDER>=3,HTrans+=Table[-(1/2)*Sum[(V[[m,l]]*V[[l,m'']]*V[[m'',m']])/((EN[[m']]-EN[[l]])*(EN[[m'']]-EN[[l]]))+(V[[m,m'']]*V[[m'',l]]*V[[l,m']])/((EN[[m]]-EN[[l]])*(EN[[m'']]-EN[[l]])),{l,L},{m'',M}]+1/2*Sum[V[[m,l]]*V[[l,l']]*V[[l',m']]*(1/((EN[[m]]-EN[[l]])*(EN[[m]]-EN[[l']]))+1/((EN[[m']]-EN[[l]])*(EN[[m']]-EN[[l']]))),{l,L},{l',L}],{m,M},{m',M}]];
HTrans
]


SWTransS[EN_,V_,M_,ORDER_:2]:=Module[{L},
L=Length[EN];
S=Table[0,L,L];
If[ORDER==1,S[[1;;M,M+1;;L]]+=Table[-(V[[m,l]]/(EN[[m]]-EN[[l]])),{m,1,M},{l,M+1,L}]];
If[ORDER==2,S[[1;;M,M+1;;L]]+=Table[1/(EN[[m]]-EN[[l]])*(Sum[(V[[m,m']]*V[[m',l]])/(EN[[m']]-EN[[l]]),{m',1,M}]-Sum[(V[[m,l']]*V[[l',l]])/(EN[[m]]-EN[[l']]),{l',M+1,L}]),{m,1,M},{l,M+1,L}]];
If[ORDER==3,S[[1;;M,M+1;;L]]+=Table[1/(EN[[m]]-EN[[l]])*(-Sum[(V[[m,m'']]*V[[m'',m']]*V[[m',l]])/((EN[[m'']]-EN[[l]])*(EN[[m']]-EN[[l]])),{m',1,M},{m'',1,M}]-Sum[(V[[m,l']]*V[[l',l'']]*V[[l'',l]])/((EN[[m]]-EN[[l'']])*(EN[[m]]-EN[[l']])),{l',M+1,L},{l'',M+1,L}]+Sum[(V[[m,m']]*V[[m',l']]*V[[l',l]])/((EN[[m']]-EN[[l]])*(EN[[m']]-EN[[l']])),{m',1,M},{l',M+1,L}]+Sum[(V[[m,m']]*V[[m',l']]*V[[l',l]])/((EN[[m]]-EN[[l']])*(EN[[m']]-EN[[l']])),{m',1,M},{l',M+1,L}]+1/3*Sum[(V[[m,l']]*V[[l',m']]*V[[m',l]])/((EN[[m']]-EN[[l']])*(EN[[m']]-EN[[l]])),{m',1,M},{l',M+1,L}]+1/3*Sum[(V[[m,l']]*V[[l',m']]*V[[m',l]])/((EN[[m]]-EN[[l']])*(EN[[m']]-EN[[l']])),{m',1,M},{l',M+1,L}]+2/3*Sum[(V[[m,l']]*V[[l',m']]*V[[m',l]])/((EN[[m]]-EN[[l']])*(EN[[m']]-EN[[l]])),{m',1,M},{l',M+1,L}]),{m,1,M},{l,M+1,L}]];
S=S-Transpose[S];
S
]


(* ::Input::Initialization:: *)
PercentError[actual_,measured_]:=(measured-actual)/actual*100;
