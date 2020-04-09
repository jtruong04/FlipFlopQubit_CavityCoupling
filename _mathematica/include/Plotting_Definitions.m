(* ::Package:: *)

(* ::Input::Initialization:: *)
(* Plotting Def *)
<<MaTeX`
UBBlue=RGBColor["#005bbb"];
UBHayesHallWhite=RGBColor["#ffffff"];
UBLetchworthAutumn=RGBColor["#e56a54"];
UBSolarStrand=RGBColor["#ffc72c"];
UBGreinerGreen=RGBColor["#ebec00"];
UBLakeLaSalle=RGBColor["#00a69c"];
UBCapenBrick=RGBColor["#990000"];
UBBronzeBuffalo=RGBColor["#ad841f"];
UBOlmstedGreen=RGBColor["#6da04b"];
UBNiagaraWhirlpool=RGBColor["#006570"];
UBVictorEBull=RGBColor["#2f9fd0"];
UBHarrimanBlue=RGBColor["#002f56"];
UBBairdPoint=RGBColor["#e4e4e4"];
UBPutnamGray=RGBColor["#666666"];
UBColors={UBHarrimanBlue,UBCapenBrick,UBOlmstedGreen,UBBronzeBuffalo,UBBlue,UBLetchworthAutumn,UBLakeLaSalle,UBSolarStrand,UBVictorEBull};
UBDark=Blend[{UBCapenBrick,White,UBHarrimanBlue},#]&;
UBLight=Blend[{UBLetchworthAutumn,White,UBBlue},#]&;
texStyle={FontFamily->"Latin Modern Roman",FontSize->24};
plotStyle={Frame->True,
BaseStyle->texStyle,PlotStyle->Table[{UBColors[[i]],Thickness[0.005]},{i,9}],
FrameStyle->Directive[FontSize->18,Thickness[0.005]],
GridLines->Automatic,
FrameTicksStyle->Directive[FontSize->24,FontFamily->"Times New Roman"],
ImageSize->Large};
