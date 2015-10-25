(* ::Package:: *)

(* ::Input:: *)
(*Manipulate[Plot[Evaluate[y[x]/.NDSolve[{y'[x] ==-(g^2/k)(1-(y[x])^2),y[0]== 0.999},y,{x,  0,2}]], {x, 0, 2}, PlotRange -> Full, AxesLabel -> {"time", "zeta"} ], {g, 0.001, 10, 0.001}, {k, 0.001, 10, 0.001}]*)


(* ::Input:: *)
(*Plot[Sqrt[1-(2e/0.111)^2], {e, 0, (0.111/2)}, PlotRange -> Full, AxesLabel ->{"Epsilon", "Zeta"} ]*)


(* ::Input:: *)
(*Manipulate[{ToRules[Quiet@Reduce[alpha == -I*e(1/(kappa-I*(g/(2*Abs[alpha])))), alpha, Complexes]]}, {e, 0, 100}, {kappa, 0, 10, 0.01}, {g, 0,10, 0.01}]*)


(* ::Input:: *)
(*fplus[e_] :=alpha /.{ToRules[Quiet@Reduce[alpha == -I*e(1/(0.01+I*(0.5/(2*Abs[alpha])))), alpha, Complexes]]}[[1]]*)
(*fminus[e_]:= alpha /.{ToRules[Quiet@Reduce[alpha == -I*e(1/(0.01-I*(0.5/(2*Abs[alpha])))), alpha, Complexes]]}[[2]]*)


(* ::Input:: *)
(*bplusplus[e_] := (fplus[e]/(2*Abs[fplus[e]]))*)
(*bminusminus[e_] := -(fminus[e]/(2*Abs[fminus[e]]))*)


(* ::Input:: *)
(*Plot[{Arg[bplusplus[e]]/Pi, Arg[bminusminus[e]]/Pi}, {e, 0.250001, 5} ,  PlotRange -> Full, AxesLabel -> {"Beta/Pi", "Drive Strength"}]*)
(**)
(**)


(* ::Input:: *)
(*f[e_] := PolarPlot[0.5, {theta, Arg[bplusplus[0.25000001]], Arg[bplusplus[0.25000002+e]]}, PlotRange->All, PlotStyle -> Orange, AxesLabel -> {"Re[beta]", "Im[beta]"}]*)


(* ::Input:: *)
(* *)


(* ::Input:: *)
(*g[e_]:= PolarPlot[0.5, {theta, Arg[bminusminus[0.25000001]], Arg[bminusminus[0.25000002+e]]}, PlotRange->All, AxesLabel -> {"Re[beta]", "Im[beta]"}]*)


(* ::Input:: *)
(*q =Table[Show[g[e], f[e]], {e, 0, 4, 0.05}]*)


(* ::Input:: *)
(*Export["PhaseCurve.avi",q] *)


(* ::Input:: *)
(*SystemOpen["PhaseCurve.avi"]*)
