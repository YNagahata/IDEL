(* ::Package:: *)

(* import energy data *)
ImportEQ[InFileEQ_] := Import[InFileEQ][[All, 1 ;; 2]];
ImportTS[InFileTS_] := Import[InFileTS][[All, 1 ;; 4]];

(* build a coefficient matrix of the 1st order rate equation *)
RateConstants[EQs_, TSs_, PhysConst_, acc_, Tc_] := Module[{kBTh, RT, MxD, MxM, fk,dk,Nk,Ngpe,id,j1,j2},
  kBTh = Simplify[kB*(Tc Kelvin+Tk0)/h/.PhysConst];
  RT = Simplify[R*(Tc Kelvin+Tk0)/.PhysConst];
  Nk = Length[EQs];
  (*fk and dk*)
  id = {}; Do[id = Append[id, EQs[[i, 1]] -> i], {i, 1, Nk}];
  Do[fk[j1, j2] := 0/Sec, {j1, 1, Nk}, {j2, 1, Nk}];
  Do[dk[j1] := 0/Sec, {j1, 1, Nk}];
  Do[
    j1 = Lookup[id, TSs[[i, 2]]];
    j2 = Lookup[id, TSs[[i, 3]]];
    If[j1==j2,Continue]
    If[j1<=Nk&&j2<=Nk,
      fk[j1, j2] = fk[j1, j2] +
                   kBTh*Exp[-SetPrecision[TSs[[i, 4]] - EQs[[j1, 2]], acc]*10^3 Joule/Mole/(RT)];
      fk[j2, j1] = fk[j2, j1] +
                   kBTh*Exp[-SetPrecision[TSs[[i, 4]] - EQs[[j2, 2]], acc]*10^3 Joule/Mole/(RT)];
      ,
      If[j1>Nk,
        dk[j2] = dk[j2] +
                 kBTh*Exp[-SetPrecision[TSs[[i, 4]] - EQs[[j2, 2]], acc]*10^3 Joule/Mole/(RT)];
        ,
        dk[j1] = dk[j1] +
                 kBTh*Exp[-SetPrecision[TSs[[i, 4]] - EQs[[j1, 2]], acc]*10^3 Joule/Mole/(RT)];
      ]
    ]
    ,{i, 1, Length[TSs]}];
  Clear[j1, j2];
  (*convert to list*)
  Do[MxD[i, j] = 0, {j, 1, Nk}, {i, 1, Nk}];
  Do[MxD[i, i] = Sum[fk[i, j], {j, 1, Nk}] + dk[i], {i, 1, Nk}];
  Do[MxM[i, j] = fk[i, j] - MxD[i, j], {j, 1, Nk}, {i, 1, Nk}];
  Transpose[Array[MxM, {Nk, Nk}]]
  ];

(* Pairwise epsilon-indistinguishability list *)
Indistinguish[\[Epsilon]_, lsM_, EQs_, EVal_, EVecR_, EVecL_,acc_] := Module[{Log1E,Nk,LInfLogP,k,lsRt},
  Log1E = Log[1+\[Epsilon]];
  Nk = Length[EVal];
  (* LInfLogP *)
  LInfLogP[j1_, j2_, Dt_] :=
  Max[Abs[
     Log[Sum[
       Exp[EVal[[k]] Dt]
       EVecR\[ConjugateTranspose][[k, All]]*EVecL[[k, j1]]
       ,{k, 1, Nk}]]
     - Log[Sum[
       Exp[EVal[[k]] Dt] EVecR\[ConjugateTranspose][[k, All]]*EVecL[[k, j2]]
       , {k, 1, Nk}]]
     ]];
  (* lsRt *)
  k = 1;
  lsRt = {};
  Do[
    AppendTo[lsRt,
      {EQs[[j1, 1]], EQs[[j2, 1]], Dt /. FindRoot[
        LInfLogP[j1, j2, Dt] == Log1E,
        {Dt, 10^(-14), 0, Infinity},
        PrecisionGoal -> 10, 
        AccuracyGoal -> 10,
        WorkingPrecision -> acc]
      }
    ];
    k = k + 1;
    , {j1, 1, Nk}, {j2, j1 + 1, Nk}];
  Clear[k, LInfLogP];
  Sort[lsRt, #1[[3]] < #2[[3]] &]
  ]

(* Lump a Matrix K by Nagahata's Exact Lumping *)
LumpExact[Blng_, K_, EVecR_, EVecL_] := Module[{CHat, FHat, LHat, Nk},
  Nk = Length[K];
  (* CHat *)
  CHat = Table[0, Max[Blng], Length[Blng]];
  Do[
    If[j == Blng[[i]], CHat[[j]][[i]] = 1;];
    ,{i, 1, Length[Blng]} ,{j, 1, Max[Blng]}];
   (* FHat *)
   FHat = Sum[
     KroneckerProduct[
      Transpose[EVecR][[k, All]], EVecL[[k, All]]
      ]
     , {k, Nk - Length[CHat] + 1, Nk}];
   LHat = CHat.FHat;
   LHat.K.PseudoInverse[LHat]
   ];

(* Lump a Matrix K by GPE *)
LumpGenPE[Blng_, EQs_, TSs_, PhysConst_, acc_, Tc_] := Module[{MxD, MxM, ZTS, ZEQ, kBTh, RT, Nk, Ngpe, id, j1, j2},
  kBTh = Simplify[kB*(Tc Kelvin+Tk0)/h/.PhysConst];
  RT = Simplify[R*(Tc Kelvin+Tk0)/.PhysConst];
  Nk = Length[EQs];
  Ngpe = Max[Blng];
  id = {}; Do[id = Append[id, EQs[[i, 1]] -> i], {i, 1, Nk}];
  (*Z_EQ*)
  Do[ZEQ[j1] := 0, {j1, 1, Ngpe}];
  Do[
    ZEQ[Blng[[j1]]] = ZEQ[Blng[[j1]]] + Exp[-SetPrecision[EQs[[j1, 2]], acc]*10^3 Joule/Mole/(RT)]
  , {j1, 1, Nk}];
  (*Z_TS*)
  Do[ZTS[j1, j2] := 0/Sec, {j1, 1, Ngpe}, {j2, 1, Ngpe}];
  Do[
    j1 = Lookup[id, TSs[[i, 2]]];
    j2 = Lookup[id, TSs[[i, 3]]];
    If[j1 != j2 && j1 <= Nk && j2 <= Nk && Blng[[j1]] != Blng[[j2]], 
      ZTS[Blng[[j1]], Blng[[j2]]] = ZTS[Blng[[j1]], Blng[[j2]]] + Exp[-SetPrecision[TSs[[i, 4]], acc]*10^3 Joule/Mole/(RT)];
      ZTS[Blng[[j2]], Blng[[j1]]] = ZTS[Blng[[j2]], Blng[[j1]]] + Exp[-SetPrecision[TSs[[i, 4]], acc]*10^3 Joule/Mole/(RT)];
      ];
    , {i, 1, Length[TSs]}];
  Clear[j1, j2];
  Do[MxD[i, j] = 0, {j, 1, Ngpe}, {i, 1, Ngpe}];
  Do[MxD[i, i] = Sum[ZTS[i, j]/ZEQ[i], {j, 1, Ngpe}], {i, 1, Ngpe}];
  Do[MxM[i, j] = kBTh*(ZTS[i, j]/ZEQ[i] - MxD[i, j]), {j, 1, Ngpe}, {i, 1, Ngpe}];
  Transpose[Array[MxM, {Ngpe, Ngpe}]]
  ];

EigenSystemRL[Mx_] := Module[{EVal, EVecR, EVecL},
  {EVal, EVec} = Eigensystem[Mx];
  EVecR = Transpose[EVec];
  EVecL = Inverse[EVecR];
  {EVal, EVecR, EVecL}
  ];

ArrheniusParameters[EVals_, Nk_] := Module[{EValArrhenius,Vals,Coef},
  EValArrhenius = {};
  Do[
    Vals = {};Do[
      AppendTo[Vals, {Joule/Mole/((Tc Kelvin +Tk0)*R)/.PhysConst,Log[-EVals[Tc][[-1-k]] ]}];
      , {Tc, Keys[EVals]}];
    Coef = CoefficientList[Fit[Vals, {1,x}, x], x];
    AppendTo[EValArrhenius, {Exp[ Coef[[1]] ], -Coef[[2]]/1000} ]
    , {k, 1, Nk-1}];
  EValArrhenius
];

ExportArrhenius[OutFileName_, Arrhenius_,ArrheniusFull_, ValStr_] := Module[{Outs},
  Outs = {{"#", ValStr, "A / [s]", "E_a / [kJ/mol]", "RE[A] / [s]", "RE[E_a] / [kJ/mol]"}};
  Do[
    AppendTo[Outs, Join[{StringJoin[ValStr, "_", ToString[k]]}, Arrhenius[[k]], Abs[(Arrhenius[[k]]-ArrheniusFull[[k]])/ArrheniusFull[[k]] ] ] ];
    , {k, 1, Length[Arrhenius]}];
  Export[OutFileName, Outs, "TextDelimiters" -> ""]
  ];

ExportRM[K_, EQs_, Dir_, Ext_] := Module[{ks},
  ks = {};
  Do[
    If[i == j, Continue[]];
    If[(K[[i, j]] Second) == 0.0, Continue[]];
    AppendTo[ks, {EQs[[i, 1]], EQs[[j, 1]], K[[i, j]] Second}];
    , {i, 1, Length[EQs]}, {j, 1, Length[EQs]}];
  Export[FileNameJoin[{Dir, StringJoin["ks",Ext]}], ks];
];

ExportESys[EVal_, EVecR_, EVecL_, Dir_, Ext_] := Module[{},
  Export[FileNameJoin[{Dir, StringJoin["EVal",Ext]}], EVal];
  Export[FileNameJoin[{Dir, StringJoin["EVecR",Ext]}], EVecR];
  Export[FileNameJoin[{Dir, StringJoin["EVecL",Ext]}], EVecL];
];

ExportRate[K_,OutFileName_] := Module[{Rate},
  Rates = {{"#","Tc", "i", "j", "k_{i,j}"}};
  Do[
    Do[
      AppendTo[Rates, {Tc, j1, j2, K[Tc][[j1]][[j2]] Second}]
      ,{j1,1,Length[K[Tc]]},{j2,1,Length[K[Tc]]}];
  ,{Tc,Keys[K]}]
  Export[OutFileName, Rates, "TextDelimiters" -> ""];
]

ExportEVal[EVals_,EValsFull_,FileName_] := Module[{Outs},
  Do[
    OutFileName = StringJoin[FileName[[1]], ToString[-k-1], FileName[[2]]];
    Outs = {{"#","Tc","k", "lambda_k"}};
    Do[
      AppendTo[Outs, {Tc, -k-1, EVals[Tc][[k]]}];
      , {Tc, Keys[EVals]}];
    Export[OutFileName, Outs, "TextDelimiters" -> ""];
    , {k, -2, -Length[EVals[0]],-1}];
  Do[
    OutFileName = StringJoin[FileName[[1]], ToString[-k-1], ".wRE", FileName[[2]]];
    Outs = {{"#","Tc","k", "lambda_k", "RE"}};
    Do[
      AppendTo[Outs, {Tc, -k-1, EVals[Tc][[k]], Abs[(EVals[Tc][[k]]-EValsFull[Tc][[k]])/EValsFull[Tc][[k]] ]}];
      , {Tc, Keys[EVals]}];
    Export[OutFileName, Outs, "TextDelimiters" -> ""];
    , {k, -2, -Length[EVals[0]],-1}];
]