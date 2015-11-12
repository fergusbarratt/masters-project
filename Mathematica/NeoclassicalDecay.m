(* ::Package:: *)

(* ::Input:: *)
(*(*Nonexponential Decay*)*)
(*Manipulate[Plot[Evaluate[y[x]/.NDSolve[{y'[x] ==-(g^2/k)(1-(y[x])^2),y[0]== 0.999},y,{x,  0,2}]], {x, 0, 2}, PlotRange -> Full, AxesLabel -> {"time", "zeta"} ], {g, 0.001, 10, 0.001}, {k, 0.001, 10, 0.001}]*)


(* ::Input:: *)
(*Plot[Sqrt[1-(2e/0.5)^2], {e, 0, (0.5/2)}, PlotRange -> Full, AxesLabel ->{"Epsilon", "Zeta"} ]*)


(* ::Text:: *)
(* *)


(* ::Input:: *)
(*(*Neoclassical Steady State Solutions to the Maxwell Bloch equations*)*)
(**)
(*NewSign[e_] :=*)
(*Piecewise[{{Sign[e], e!=0}, {1, e==0}}]*)
(*(*Piecewise sign function making zero positive*)*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*asymeqn = n== e^2/(kappa^2+(det+(g/(2*Sqrt[n]))))*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*alphaeqnnodet[e_, plusminus_] :=*)
(* alpha == -I*drive*If[plusminus=="+", *)
(*(kappa+I*(g/(2*Abs[alpha])))^-1,*)
(*(kappa-I*(g/(2*Abs[alpha])))^-1] /.*)
(*{drive ->e, kappa->0.5, g->25}*)


(* ::Input:: *)
(*alphaeqn[e_, detuning_, plusminus_] :=alpha == -I*drive*If[plusminus =="+", (kappa-I*(det+NewSign[det]*(g^2/(Sqrt[det^2+4*g^2*Abs[alpha]^2]))))^-1, (kappa-I*(det-NewSign[det]*(g^2/(Sqrt[det^2+4*g^2*Abs[alpha]^2]))))^-1] /. {drive->e, det->detuning, g-> 25, kappa->0.5}*)


(* ::Input:: *)
(*alphaplus[e_] := alpha /.{ToRules[Quiet@Reduce[alphaeqnnodet[e, "+"], alpha, Complexes]]}[[1]]*)
(*alphaminus[e_]:= alpha /.Reverse[{ToRules[Quiet@Reduce[alphaeqnnodet[e, "-"], alpha, Complexes]]}][[1]]*)


(* ::Input:: *)
(*alphaplusdet[e_, detuning_] := alpha /.{ToRules[Quiet@Reduce[alphaeqn[e, detuning, "-"], alpha, Complexes]]}[[1]]*)
(*alphaminusdet[e_, detuning_]:= alpha /.Reverse[{ToRules[Quiet@Reduce[alphaeqn[e, detuning,  "+"], alpha, Complexes]]}][[1]]*)


(* ::Input:: *)
(*(*Reverse is for consistency between det and no det*)*)


(* ::Input:: *)
(*betaplus[e_] := (alphaplus[e]/(2*Abs[alphaplus[e]]))*)
(*betaminus[e_] := -(alphaminus[e]/(2*Abs[alphaminus[e]]))*)


(* ::Input:: *)
(* betaplusdet[e_, det_] := NewSign[det]*g*alphaplusdet[e, det]/Sqrt[det^2 + 4*g^2*Abs[alphaplusdet[e, det]]^2]/. g-> 25*)


(* ::Input:: *)
(* betaminusdet[e_, det_] := -NewSign[det]*g*alphaminusdet[e, det]/Sqrt[det^2 + 4*g^2*Abs[alphaminusdet[e, det]]^2]/. g-> 25*)


(* ::Input:: *)
(*(* Plot[{Arg[betaplus[e]]/Pi, Arg[betaminus[e]]/Pi}, {e, 0.25001, 4} ,  PlotRange \[Rule] Full, AxesLabel \[Rule] { "Drive Strength", "Beta/Pi"}] *)*)


(* ::Input:: *)
(*(* Plot[{Arg[betaplusdet[e, 0]]/Pi, Arg[betaminusdet[e, 0]]/Pi}, {e, 0.25001, 2} ,  PlotRange \[Rule] Full, AxesLabel \[Rule] { "Drive Strength", "Beta/Pi"}]*)*)


(* ::Input:: *)
(*(*Functions for mapping along a line or a matrix*)*)


(* ::Input:: *)
(*Clear[drive, detuning]*)
(*func[drive_, detuning_] := Abs[alphaminusdet[drive, detuning]]^2*)
(*gfunc[drive_, detuning_] := Abs[alphaplusdet[drive, detuning]]^2*)


(* ::Input:: *)
(*qfunc[q_]:= Apply[gfunc, q]*)


(* ::Input:: *)
(*EVals[estart_,space_, estop_, det_]:=Map[qfunc, Partition[Riffle[Range[estart, estop, space], ConstantArray[det, Length[Range[estart, estop, space]]]], 2]]*)


(* ::Input:: *)
(*DVals[e_, detstart_, space_, detstop_] := Map[qfunc, Partition[Riffle[ConstantArray[e, Length[Range[detstart, detstop, space]]], Range[detstart, detstop, space]], 2]]*)


(* ::Input:: *)
(*(*Generate Lines*)*)


(* ::Input:: *)
(*LineGenE[estart_, space_, estop_, det_]:=Partition[Riffle[Range[estart, estop, space], EVals[estart, space, estop, det]], 2]*)


(* ::Input:: *)
(*LineGenDet[e_, detstart_, space_, detstop_]:= Partition[{ConstantArray[e, Length[Range[detstart, detstop, space]]],Range[detstart, detstop, space], DVals[e, detstart, space, detstop]} ~Flatten~{2, 1},  3]*)


(* ::Input:: *)
(*MatGen[estart_, espace_, estop_, detstart_,detspace_, detstop_] :=(q = {}; starttime = 0;Do[initialtime = SessionTime[];AppendTo[q, LineGenDet[e, detstart, detspace, detstop]]; Print[IntegerPart[((e-estart)/espace)]+1, "/", IntegerPart[(estop-estart)/espace]+1, "||", starttime+=SessionTime[]-initialtime,"/",(SessionTime[]-initialtime)*(estop-estart)/espace] ,{e, estart, estop, espace}];Return[Flatten[q, 1]])*)
(**)


(* ::Input:: *)
(*Clear[data]*)


(* ::Input:: *)
(*data = MatGen[10, 0.04, 15.0, 0,  0, 0];*)


(* ::Input:: *)
(*highndata = MatGen[20, 0.05, 25, 0, 0, 0];*)


(* ::Input:: *)
(*Show[ListPlot3D[data, AxesLabel -> {"Drive", "Detuning", "Photon Number"}, ImageSize->Large, PlotRange->All, Mesh->None, PlotStyle->Directive[LightBlue,Specularity[Gray, 50]]],*)
(*Graphics3D[{Gray, PointSize[0.04], Point[{12.5, 0, gfunc[12.5, 0]}]}]]*)


(* ::Input:: *)
(*kappa = 0.5*)
(*det = 0*)
(*g = 25*)


(* ::Input:: *)
(*eqn = n== e^2/(kappa^2+(det+(g/(2*Sqrt[n]))))*)


(* ::Input:: *)
(*(* Data for n:10\[Rule]15 from full meanfield |\alpha|^2 and approximation*)*)


(* ::Input:: *)
(*Show[Plot[n/.ToRules[Quiet@Reduce[eqn/.e->drive, n]], {drive, 10, 15}, PlotRange -> All], ListPlot[Take[data,All,{1, -1, 2}]]]*)


(* ::Input:: *)
(*(* Data for n:20\[Rule]25 from full meanfield |\alpha|^2 and approximation*)*)


(* ::Input:: *)
(*Show[Plot[n/.ToRules[Quiet@Reduce[eqn/.e->drive, n]], {drive, 20, 25}, PlotRange -> All], ListPlot[Take[highndata,All,{1, -1, 2}]]]*)


(* ::Input:: *)
(*(*Assymetric Lorentzians*)*)


(* ::Input:: *)
(*gamma[center_, gamma0_, skew_] := 2*gamma0/(1+Exp[skew*(x-center)])*)


(* ::Input:: *)
(*Lorentzian[strength_, gamma_, center_] := (strength*1/Pi)*(0.5*gamma)/((x-center)^2+((1/2)*gamma)^2)*)
(**)
(**)


(* ::Input:: *)
(*Assym[strength_,center_, gamma0_, skew_] := Lorentzian[strength, gamma[center, gamma0, skew], center]*)


(* ::Input:: *)
(*Plot[Assym[1, 1, 1, 1], {x, -1, 1}]*)


(* ::Text:: *)
(* *)


(* ::Input:: *)
(*(*Animation of Development of Phase Bistability w/ Increasing Drive*)*)


(* ::Input:: *)
(*up[e_] := PolarPlot[0.5, {theta, Arg[betaplus[0.25000001]], Arg[betaplus[0.25000002+e]]}, PlotRange->All, PlotStyle -> Orange, AxesLabel -> {"Re[beta]", "Im[beta]"}]*)


(* ::Input:: *)
(*down[e_]:= PolarPlot[0.5, {theta, Arg[betaminus[0.25000001]], Arg[betaminus[0.25000002+e]]}, PlotRange->All, PlotStyle->Blue, AxesLabel -> {"Re[beta]", "Im[beta]"}]*)


(* ::Input:: *)
(*Animate[Show[up[e], down[e]], {e, 0, 5}]*)


(* ::Text:: *)
(* *)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Length[Range[5]]*)


(* ::Input:: *)
(*Partition[Riffle[Range[5], ConstantArray[2, 5]],2]*)


(* ::Input:: *)
(*Manipulate[Plot[a*x, {x, 0, 5}], {a, 0, 5}]*)
