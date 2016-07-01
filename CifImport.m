(* ::Package:: *)

(* :Context: CifImport` *)

(* :Author: 
    Bianca Eifert and Christian Heiliger, 
    Theoretical Solid State Physics, 
    Institute for Theoretical Physics, 
    Justus Liebig University Giessen, Germany *)

(* :Copyright: Copyright (c) 2016 Bianca Eifert and Christian Heiliger *)

(* :Summary: The CifImport package contains an import function for CIF files. *)

(* :Discussion: 
    CifImport comes "as is", with no warranties whatsoever. Please do your own research on the CIF format.
    In particular, read this: http://oldwww.iucr.org/iucr-top/cif/faq/#ipr
    and this: http://www.iucr.org/__data/assets/pdf_file/0019/22618/cifguide.pdf
    This package is distributed under the terms of the GNU General Public License (Version 3, 29 June 2007) *)

(* :Package Version: 1.0.0 *)

(* :Mathematica Version: 10.0 *)


If[TrueQ[$VersionNumber<10],Print["Sorry, CifImport only works with version 10+ of the Wolfram Language."];Abort[]];


BeginPackage["CifImport`"];
CifImport::usage="CifImport[file] imports a crystal structures from a CIF file, returning multiple versions of the structure as explained below. \
(CifImport can only handle files with one structure, but files with multiple structures can easily be split with any text editor.)

CIF files do not store all atoms in a structure explicitly. Instead, they contain only those atoms that are not connected to one another by \
symmetry operations, or the atoms of the asymmetric unit cell. The file then also contains all symmetry operations of the cell, so that the \
remaining atoms can be constructed. This is what CifImport will do. Additionally, the atom types stored in a CIF file may not be actual chemical \
types, they can also be a concatenation of the chemical type and an integer unique to the atom. Once the other atoms have been added by symmetry, \
the types will be unique to the atom's symmetry position. Sometimes these \"non-chemical\" types are helpful to have, although in most cases you \
will want the chemical types. To allow access to all of this data, CifImport returns four structure datasets for each crystal structure: the \
asymmetric and the complete cell, and the chemical and non-chemical types for each. (These four datasets may or may not actually be different \
from one another depending on the input file, but you always get all four of them.) The comment in the last entry of each structure will let you \
know which version of the structure you're dealing with. The final structure in the list is always the complete cell with chemical atom types, \
so Last[CifImport[file]] will give you the most physically meaningful version of the structure.";
Begin["`Private`"];


(*construct lattice vectors from lengths and angles*)
LatticeVectors[lengths_,angles_]:=
Module[{a,b,c,\[Alpha],\[Beta],\[Gamma],veca,vecb,vecc,c2,c3},
{a,b,c}=N[lengths];
{\[Alpha],\[Beta],\[Gamma]}=N[angles]/180*\[Pi];
veca=a*{1,0,0};
vecb=b*{Cos[\[Gamma]],Sin[\[Gamma]],0};
vecc=c*{Cos[\[Beta]],c2,c3};
c2=(Cos[\[Alpha]]-Cos[\[Gamma]]*Cos[\[Beta]])/Sin[\[Gamma]];
c3=Sqrt[1-(Cos[\[Beta]])^2-(c2)^2];
Chop[{veca,vecb,vecc}]
];


(*find all atoms symmetrically equivalent to one specific atom*)
SymEquiv[sym_,coord_,conf_,index_]:=
Module[{xtemp,ytemp,ztemp,symmetries,unique,symout},
(*apply all symmetries found for the system to the atom, returning all symmetrically equivalent atoms*)
symmetries=sym/.{CifImport`Private`x->coord[[index,1]],CifImport`Private`y->coord[[index,2]],CifImport`Private`z->coord[[index,3]]};
(*delete duplicate atoms*)
unique=Round[#,.000001]-Floor[Round[#,.000001]]&/@symmetries;
unique=DeleteDuplicates[unique,(Abs[#1[[1]]-#2[[1]]]<.0001&&Abs[#1[[2]]-#2[[2]]]<.0001&&Abs[#1[[3]]-#2[[3]]]<.0001)&];
(*return atoms and their types*)
symout={#,conf[[index]]}&/@unique;
symout
];


(*find all unique atoms from the pool of symmetrical equivalents*)
SymAllAtoms[sym_,coord_,conf_]:=
Module[{symequiv},
(*run SymEquiv over all values of 'index', then delete duplicate atoms*)
symequiv=Flatten[SymEquiv[sym,coord,conf,#]&/@Range[Length[conf]],1];
DeleteDuplicates[symequiv,(Abs[#1[[1,1]]-#2[[1,1]]]<.0001&&Abs[#1[[1,2]]-#2[[1,2]]]<.0001&&Abs[#1[[1,3]]-#2[[1,3]]]<.0001&&TrueQ[#1[[2]]==#2[[2]]])&]
];


(*read a CIF file: *)
CifImport[file_String/;(FileExistsQ[file]&&TrueQ[FileExtension[file]=="cif"])]:=
Module[{data,unstring,lengths,angles,lattvec,coord,coordall,confall,confnonchem,chemnonchem,conf,knowntypes,chem,
symopstrings,splitstrings,sym,symequiv,explicitinput,fullinput,symequivnonchem,explicitinputnonchem,fullinputnonchem,coordnonchem,whichlabel,label},

data=Association[Import[file]];

(*string to expression, and remove uncertain digits: *)
unstring[x_]:=If[MatchQ[x,_String],ToExpression[StringSplit[x,{"(",")","[","]","{","}",","}][[1]]],x];

(*lattice vectors*)
lengths={data["_cell_length_a"],data["_cell_length_b"],data["_cell_length_c"]};
angles={data["_cell_angle_alpha"],data["_cell_angle_beta"],data["_cell_angle_gamma"]};
lattvec=LatticeVectors[unstring/@lengths,unstring/@angles];

(*atoms and their types*)
coord=Transpose[{unstring/@data["_atom_site_fract_x"],unstring/@data["_atom_site_fract_y"],unstring/@data["_atom_site_fract_z"]}];
conf=data["_atom_site_label"];
confnonchem=conf;
chemnonchem=DeleteDuplicates[confnonchem];
confnonchem=ReplacePart[confnonchem,Flatten@Table[#->ii&/@(Flatten@Position[confnonchem,DeleteDuplicates[confnonchem][[ii]]]),{ii,1,Length@DeleteDuplicates[confnonchem]}]];
(*separate proper elements from appendices: *)
knowntypes={"Wat","Uuu","Uut","Uus","Uuq","Uup","Uuo","Uuh","Uub","O-H","Zr","Zn","Yb","Xe","Tm","Tl","Ti","Th","Te","Tc","Tb","Ta","Sr","Sn","Sm","Si","Sg","Se",
"Sc","Sb","Ru","Rn","Rh","Rg","Rf","Re","Rb","Ra","Pu","Pt","Pr","Po","Pm","Pd","Pb","Pa","Os","Np","No","Ni","Ne","Nd","Nb","Na","Mt","Mo","Mn","Mg","Md","Lu",
"Lr","Li","La","Kr","Ir","In","Hs","Ho","Hg","Hf","He","Ge","Gd","Ga","Fr","Fm","Fe","Eu","Es","Er","Dy","Ds","Db","Cu","Cs","Cr","Co","Cn","Cm","Cl","Cf","Ce",
"Cd","Ca","Br","Bk","Bi","Bh","Be","Ba","Au","At","As","Ar","Am","Al","Ag","Ac","Y","W","V","U","T","S","P","O","N","K","I","H","F","D","C","B",""};
Table[(conf[[ii]]=Select[knowntypes,StringMatchQ[conf[[ii]],#<>"*"]&][[1]]),{ii,Range[Length[conf]]}];
chem=DeleteDuplicates[conf];
conf=ReplacePart[conf,Flatten@Table[#->ii&/@(Flatten@Position[conf,DeleteDuplicates[conf][[ii]]]),{ii,1,Length@DeleteDuplicates[conf]}]];

(*substance name: *)
whichlabel=Select[{"_chemical_name_mineral","_amcsd_formula_title","_pd_phase_name"},KeyExistsQ[data,#]&];
label=If[TrueQ[whichlabel=={}],FileBaseName[file],data[whichlabel[[1]]]];

(*asymmetric units: *)
explicitinput=<|"lattice"->lattvec,"atomcoords"->coord,"atomtypes"->conf,"chemical"->chem,"name"->label,"file"->AbsoluteFileName[file],"comment"->"atoms of the asymmetric unit with chemical types"|>;
explicitinputnonchem=<|"lattice"->lattvec,"atomcoords"->coord,"atomtypes"->confnonchem,"chemical"->chemnonchem,"name"->label,"file"->AbsoluteFileName[file],"comment"->"atoms of the asymmetric unit with NONchemical types"|>;

(*read symmetries: *)
symopstrings=data[Select[{"_space_group_symop_operation_xyz","_symmetry_equiv_pos_as_xyz"},KeyExistsQ[data,#]&][[1]]];
splitstrings=StringSplit[#,","]&/@symopstrings;
splitstrings=StringReplace[#,{"x"->"CifImport`Private`x","y"->"CifImport`Private`y","z"->"CifImport`Private`z"}]&/@splitstrings;
sym=ToExpression[splitstrings];

(*add symmetrically equivalent atoms for chemical structures: *)
symequiv=SymAllAtoms[sym,coord,conf];
coordall=symequiv[[All,1]];
confall=symequiv[[All,2]];
fullinput=<|"lattice"->lattvec,"atomcoords"->coordall,"atomtypes"->confall,"chemical"->chem,"name"->label,"file"->AbsoluteFileName[file],"comment"->"complete cell with chemical types"|>;
(*add symmetrically equivalent atoms for NONchemical structures: *)
If[TrueQ[chem==chemnonchem],
fullinputnonchem=fullinput;,
symequivnonchem=SymAllAtoms[sym,coord,confnonchem];
coordnonchem=symequivnonchem[[All,1]];
confnonchem=symequivnonchem[[All,2]];
fullinputnonchem=<|"lattice"->lattvec,"atomcoords"->coordnonchem,"atomtypes"->confnonchem,"chemical"->chemnonchem,"name"->label,"file"->AbsoluteFileName[file],"comment"->"complete cell with NONchemical types"|>;
];

{explicitinputnonchem,fullinputnonchem,explicitinput,fullinput}
];


End[];
EndPackage[];
