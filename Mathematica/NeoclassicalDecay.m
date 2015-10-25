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
(*alphaeqnnodet[e_, plusminus_] :=*)
(* alpha == -I*drive*If[plusminus=="+", *)
(*(kappa+I*(g/(2*Abs[alpha])))^-1,*)
(*(kappa-I*(g/(2*Abs[alpha])))^-1] /.*)
(*{drive ->e, kappa->0.01, g->0.5}*)


(* ::Input:: *)
(*alphaeqn[e_, detuning_, plusminus_] :=alpha == -I*drive*If[plusminus =="+", (kappa-I*(det+NewSign[det]*(g^2/(Sqrt[det^2+4*g^2*Abs[alpha]^2]))))^-1, (kappa-I*(det-NewSign[det]*(g^2/(Sqrt[det^2+4*g^2*Abs[alpha]^2]))))^-1] /. {drive->e, det->detuning, g-> 0.5, kappa->0.01}*)


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
(* betaplusdet[e_, det_] := NewSign[det]*g*alphaplusdet[e, det]/Sqrt[det^2 + 4*g^2*Abs[alphaplusdet[e, det]]^2]/. g-> 0.5*)


(* ::Input:: *)
(* betaminusdet[e_, det_] := -NewSign[det]*g*alphaminusdet[e, det]/Sqrt[det^2 + 4*g^2*Abs[alphaminusdet[e, det]]^2]/. g-> 0.5*)


(* ::Input:: *)
(*Plot[{Arg[betaplus[e]]/Pi, Arg[betaminus[e]]/Pi}, {e, 0.25001, 4} ,  PlotRange -> Full, AxesLabel -> { "Drive Strength", "Beta/Pi"}]*)


(* ::Input:: *)
(*Plot[{Arg[betaplusdet[e, 1]]/Pi, Arg[betaminusdet[e, 1]]/Pi}, {e, 0.25001, 2} ,  PlotRange -> Full, AxesLabel -> { "Drive Strength", "Beta/Pi"}]*)


(* ::Input:: *)
(*drive = 100000;*)
(*detuning=0.01;*)
(*Arg[betaminusdet[drive, detuning]]/Pi*)
(*Arg[betaplusdet[drive, detuning]]/Pi*)


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


(* ::Input:: *)
(**)
