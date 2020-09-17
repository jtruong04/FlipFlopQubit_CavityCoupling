(* ::Package:: *)

(* ::Input::Initialization:: *)
\[Sigma][input__]:=
Module[{x,y},
x = List[input];
y =PauliMatrix[x[[1]]/.{"I"->0,"X"->1,"Y"->2,"Z"->3}];
For[i=2,i<= Length[x],i++,
y=KroneckerProduct[y,PauliMatrix[x[[i]]/.{"I"->0,"X"->1,"Y"->2,"Z"->3}]];
];
y
];
CoefficientTransformation[n_,inv_:False]:=CoefficientTransformation[n,inv]=Module[{A,B,p,row},
A=Table[0,4^n,4^n];
B=Table[0,n+1];
p=1;
row=1;
While[B[[n+1]]==0,
A[[row]]=Flatten[\[Sigma]@@Reverse[B[[-n-1;;-2]]]];
row++;
B[[1]]++;
While[B[[p]]==4,
B[[p]]=0;
p++;
B[[p]]++;
If[B[[p]]!=4,p=1];
];
];
A=Transpose[A];
If[inv,Return[Inverse[A]],Return[A]];
];
PauliToMatrix[x_]:=PauliToMatrix[x]=Module[{n},
n=If[SquareMatrixQ[x],Log[2,Dimensions[x][[1]]],Log[4,Length[x]]];
Return[ArrayReshape[CoefficientTransformation[n].Flatten[x],{2^n,2^n}]];
];
MatrixToPauli[x_]:=MatrixToPauli[x]=Module[{n},
n=If[SquareMatrixQ[x],Log[2,Dimensions[x][[1]]],Log[4,Length[x]]];
Return[CoefficientTransformation[n,True].Flatten[x]];
];


commutator[a_,b_]:=a.b-b.a
