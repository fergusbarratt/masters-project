(* ::Package:: *)

(* ::Input:: *)
(*Manipulate[Plot3D[Exp[-Abs[((x+I*y)-(p+I*q))^2]],{p,-4,4}, {q,-4,4}, PlotRange -> Full, PlotPoints -> 50, ImageSize -> Large, AxesLabel -> {"X_1", "X_2", "Q(X1, X2)"}], {x, 0, 4, 1}, {y, 0, 4, 1}]*)


(* ::Input:: *)
(*Manipulate[Plot3D[ Exp[-(1/2) ((p-x)^2+(q-y)^2)]/\[Pi],{p,-10,10},{q,-10,10},PlotStyle->RGBColor[0.`,0.71`,0.48`],PlotRange->Full, PlotPoints-> 75, ImageSize -> Large, AxesLabel -> {"X_1", "X_2", "W(X1, X2)"}],{x,-10, 10, 1}, {y, -10, 10, 1} ]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Manipulate[Plot3D[((Abs[p+I*q]^(2*k)/(k!)^2))*Exp[-Abs[(p+I*q)]^2],{p,-4,4},{q,-4,4},PlotRange->Full, PlotPoints-> 75, ImageSize ->Large, AxesLabel -> {"X_1", "X_2", "Q(X1, X2)"}], {k, 0, 5, 1}]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Manipulate[Plot3D[(2/Pi)*(-1)^n * LaguerreL[n, 4*Sqrt[x^2+y^2]]*Exp[-2*(x^2+y^2)], {x, -2, 2}, {y, -2, 2}, PlotRange->Full, PlotPoints-> 80,ImageSize -> Large,PlotStyle->RGBColor[0.`,0.71`,0.48`], AxesLabel -> {"X_1", "X_2", "W(X1, X2)"}], {n, 0, 10, 1}]*)


(* ::Input:: *)
(**)
