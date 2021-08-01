///////////////////////////////////////////////////////////////////////////////
/////////    Copyright 2021 Jaume Pujol and Merc� Villanueva           ///////
/////////                                                               ///////
/////////    This program is distributed under the terms of GNU         ///////
/////////               General Public License                          ///////
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

/************************************************************/
/*                                                          */
/* Project name: Qary nonlinear codes in MAGMA              */
/* File name: QaryCodes_ConstructionsBiplop.m               */
/*                                                          */
/* Comment: Package developed within the CCSG group         */
/*                                                          */
/* Authors: Biplop Day and Mercè Villanueva                 */
/*                                                          */
/* Revision version and last date: version 1.0 18/03/2021   */
/*                                         1.1 02/06/2021   */
/*                                                          */
/************************************************************/
//Uncomment freeze when package finished
//freeze;

/* PACKAGE VERSION */
intrinsic QaryCodes_Extension_version() -> SeqEnum
{Return the current version of this package.}
    version := [1,1];
    return version;
end intrinsic;

import "QaryCodes_Core.m": UpdateMinimumDistance;
import "QaryCodes_Core.m": UpdateMinimumDistanceLowerBound;
import "QaryCodes_Core.m": UpdateMinimumDistanceUpperBound;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////      CONSTRUCTION OF QARY CODES                                 ///////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/******************************************************************/
/*                                                                */
/* Function name: zeroPositions                                   */
/* Parameters: seqVectors, seqCoordPositions                      */
/* Function description: Given a sequence of vectors, seqVectors, */
/*   and a sequence of coordinate positions, seqCoordPositions,   */
/*   return a sequence with the coordinate positions where all    */
/*   the vectors in seqVectors have a zero.                       */                                   
/* Input parameters description:                                  */
/*   - seqVectors: A sequence of vectors of the same length       */
/*   - seqCoordPositions: A sequence of coordinate positions      */
/* Output parameters description:                                 */
/*   - A sequence of coordinate positions,                        */
/*     if there is not any, returns an empty sequence             */
/*                                                                */
/******************************************************************/
zeroPositions := function(seqVectors, seqCoordPositions)
    for c in seqVectors do 
        for i in seqCoordPositions do
            if c[i] ne 0 then
                seqCoordPositions := Exclude(seqCoordPositions, i);
            end if;
        end for;

        if #seqCoordPositions eq 0 then 
            break;
        end if;
    end for;
    return seqCoordPositions;
end function;

zeroPositionsBF := function(C)
    listLength := [1..(C`Length)];
    for c in Set(C) do 
        for i in listLength do
            if c[i] ne 0 then
                listLength := Exclude(listLength, i);
            end if;
        end for;

        if #listLength eq 0 then 
            break;
        end if;
    end for;
    return listLength;
end function;

/****************************************************************/
/*                                                              */
/* Function name: PadCode                                       */
/* Parameters: C, n                                             */
/* Function description: Add n zeros to the end of each         */
/*   codeword of C.                                             */
/* Input parameters description:                                */
/*   - C: A qary code                                           */
/*   - n: Integer with the number of zeros to be added          */
/* Output parameters description:                               */
/*   - The q-ary code                                           */
/*                                                              */
/* Signature: (<CodeFld> C, <RngIntElt> n) -> CodeFld           */
/*                                                              */
/****************************************************************/
intrinsic PadCodeBF(C::CodeFld, n::RngIntElt) -> CodeFld
{Add n zeros to the end of each codeword of C.}
    require  n ge 0: "n must be non negative";

    if n eq 0 then 
        return C;
    end if;

    V := VectorSpace(C`BaseField, C`Length + n);
    padCode := QaryCode([V!(Eltseq(c) cat [0^^n]) : c in Set(C)]);

    padCode`PAutSubgroup := DirectProduct(C`PAutSubgroup, Sym(n));
    padCode`IsPAutGroup := C`IsPAutGroup and (#zeroPositionsBF(C) eq 0);
   
    UpdateMinimumDistanceLowerBound(~padCode, C`MinimumDistanceLowerBound);
    UpdateMinimumDistanceUpperBound(~padCode, C`MinimumDistanceUpperBound);

    return padCode;

end intrinsic;

intrinsic PadCode(C::CodeFld, n::RngIntElt) -> CodeFld
{Add n zeros to the end of each codeword of C.}
    require  n ge 0: "n must be non negative";

    if n eq 0 then 
        return C;
    end if;

    V := VectorSpace(C`BaseField, C`Length + n);
    kernelPadCode := PadCode(C`Kernel, n);
    cosetRepPadCode := [ V!(Eltseq(c) cat [0^^n]) : c in C`CosetRepresentatives];
    padCode := QaryCode(kernelPadCode, cosetRepPadCode : IsFinalKernel := true);
    
    padCode`PAutSubgroup := DirectProduct(C`PAutSubgroup, Sym(n));
    if C`IsPAutGroup eq false then
        padCode`IsPAutGroup := false;
    elif C`IsLinear then                //C = C`Kernel
        padCode`IsPAutGroup := (#zeroPositions(Rows(GeneratorMatrix(C`Kernel)), [1..(C`Length)]) eq 0);
    elif Dimension(C`Kernel) eq 0 then  //C = C`CosetRepresentatives
        padCode`IsPAutGroup := (#zeroPositions(C`CosetRepresentatives, [1..(C`Length)]) eq 0);
    else                                // General case
        zeroCoordinatesKernel := zeroPositions(Rows(GeneratorMatrix(C`Kernel)), [1..(C`Length)]);
        padCode`IsPAutGroup := (#zeroCoordinatesKernel eq 0) and 
                               (#zeroPositions(C`CosetRepresentatives, zeroCoordinatesKernel) eq 0);
    end if;       
   
    UpdateMinimumDistanceLowerBound(~padCode, C`MinimumDistanceLowerBound);
    UpdateMinimumDistanceUpperBound(~padCode, C`MinimumDistanceUpperBound);
    
    return padCode;

end intrinsic;

/******************************************************************/
/*                                                                */
/* Function name: minPositive                                     */
/* Parameters: L                                                  */
/* Function description: Given a sequence of numbers, return the  */
/*   minimum positive number in the sequence.                     */                                     
/* Input parameters description:                                  */
/*   - L: A sequence of numbers                                   */
/* Output parameters description:                                 */
/*   - the minimum positive number                                */
/*                                                                */
/******************************************************************/
minPositive := function(L)
    m := 0;
    positiveL := [n : n in L | n gt 0];
    if #positiveL gt 0 then
        m := Min(positiveL);
    end if;
    return m;
end function;

/*******************************************************************/
/*                                                                 */
/* Function name: PlotkinSum                                       */
/* Parameters: C, D                                                */
/* Function description: Given q-ary codes C and D, both of the    */
/*   same length and over the same base field, construct the       */
/*   Plotkin sum of C and D. The Plotkin sum is a q-ary code that  */
/*   consists of all vectors of the form (u,u+v), where u in C and */
/*   v in D.                                                       */
/* Input parameters description:                                   */
/*   - C: A q-ary code                                             */
/*   - D: A q-ary code                                             */
/* Output parameters description:                                  */
/*   - The Plotkin sum code                                        */
/*                                                                 */
/* Remark: This function is based on Proposition 3.2.15            */
/*         in page 51 of Fanxuan Zeng's dissertation,              */
/*         "Nonlinear codes: representation, constructions,        */
/*         minimum distance computation and decoding"              */
/*                                                                 */
/* Signature: (<CodeFld> C, <CodeFld> D) -> CodeFld                */
/*                                                                 */
/*******************************************************************/
intrinsic PlotkinSumBF(C::CodeFld, D::CodeFld) -> CodeFld
{Given q-ary codes C and D, both of the same length and over the same base field, 
construct the Plotkin sum of C and D. The Plotkin sum is a q-ary code that 
consists of all vectors of the form (u,u+v), where u in C and v in D.}
    require not IsNull(C): "Argument 1 can not be an empty code";
    require not IsNull(D): "Argument 2 can not be an empty code";
    require (C`BaseField cmpeq D`BaseField): "Arguments must have same base field";
    require (C`Length eq D`Length): "Arguments must have the same length";

    V := VectorSpace(C`BaseField, 2*(C`Length)); 
    plotkinSumCodewords := [V!(Eltseq(c) cat Eltseq(c+d)) : c in Set(C), d in Set(D)];
    plotkinSumCode := QaryCode(plotkinSumCodewords);

    plotkinSumCode`PAutSubgroup := DirectProduct(C`PAutSubgroup, 
                                                 C`PAutSubgroup meet D`PAutSubgroup);
    plotkinSumCode`IsPAutGroup := false;

    UpdateMinimumDistanceLowerBound(~plotkinSumCode, 
          minPositive([2*C`MinimumDistanceLowerBound, D`MinimumDistanceLowerBound]));
    UpdateMinimumDistanceUpperBound(~plotkinSumCode, 
          minPositive([2*C`MinimumDistanceUpperBound, D`MinimumDistanceUpperBound]));

    return plotkinSumCode;

end intrinsic;

intrinsic PlotkinSum(C::CodeFld, D::CodeFld) -> CodeFld
{Given q-ary codes C and D, both of the same length and over the same base field, 
construct the Plotkin sum of C and D. The Plotkin sum is a q-ary code that 
consists of all vectors of the form (u,u+v), where u in C and v in D.}
    require not IsNull(C): "Argument 1 can not be an empty code";
    require not IsNull(D): "Argument 2 can not be an empty code";
    require (C`BaseField cmpeq D`BaseField): "Arguments must have same base field";
    require (C`Length eq D`Length): "Arguments must have the same length";
    
    V := VectorSpace(C`BaseField, 2*(C`Length));
    plotkinSumRep:= [V!(Eltseq(c) cat Eltseq(c+d)): c in CosetRepresentatives(C),
                                                    d in CosetRepresentatives(D)];                                
    plotkinSumKernel := PlotkinSum(C`Kernel, D`Kernel);                              
    plotkinSumCode := QaryCode(plotkinSumKernel, plotkinSumRep : IsFinalKernel := true);
    
    plotkinSumCode`PAutSubgroup := DirectProduct(C`PAutSubgroup,
                                                 C`PAutSubgroup meet D`PAutSubgroup);   
    plotkinSumCode`IsPAutGroup := false;

    UpdateMinimumDistanceLowerBound(~plotkinSumCode, 
          minPositive([2*C`MinimumDistanceLowerBound, D`MinimumDistanceLowerBound]));
    UpdateMinimumDistanceUpperBound(~plotkinSumCode, 
          minPositive([2*C`MinimumDistanceUpperBound, D`MinimumDistanceUpperBound]));

    return plotkinSumCode;

end intrinsic;

/****************************************************************/
/*                                                              */
/* Function name: DirectSum                                     */
/* Parameters: C, D                                             */
/* Function description: Given q-ary codes C and D, both over   */
/*   the same base field, construct the direct sum of C and D.  */
/*   The direct sum is a q-ary code that consists of all vectors*/
/*   of the form (u,v), where u in C and v in D.                */
/* Input parameters description:                                */
/*   - C: A q-ary code                                          */
/*   - D: A q-ary code                                          */
/* Output parameters description:                               */
/*   - The direct sum code                                      */
/*                                                              */
/* Signature: (<FldFin> C, <FldFin> D) -> CodeFld               */
/*                                                              */
/****************************************************************/
intrinsic DirectSumBF(C::CodeFld, D::CodeFld) -> CodeFld
{Given q-ary codes C and D, both over the same base field, construct the direct 
sum of C and D. The direct sum is a q-ary code that consists of all vectors of
the form (u,v), where u in C and v in D.}
    require not IsNull(C): "Argument 1 can not be an empty code";
    require not IsNull(D): "Argument 2 can not be an empty code";
    require (C`BaseField cmpeq D`BaseField): "Arguments must have same base field";

    V := VectorSpace(C`BaseField, C`Length + D`Length);
    directSumCodewords := [V!(Eltseq(c) cat Eltseq(d)) : c in Set(C), d in Set(D)];
    directSumCode := QaryCode(directSumCodewords);

    directSumCode`PAutSubgroup := DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
    directSumCode`IsPAutGroup := false; 

    UpdateMinimumDistanceLowerBound(~directSumCode, 
          minPositive([C`MinimumDistanceLowerBound, D`MinimumDistanceLowerBound]));
    UpdateMinimumDistanceUpperBound(~directSumCode, 
          minPositive([C`MinimumDistanceUpperBound, D`MinimumDistanceUpperBound]));

    return directSumCode;

end intrinsic;

intrinsic DirectSum(C::CodeFld, D::CodeFld) -> CodeFld
{Given q-ary codes C and D, both over the same base field, construct the direct 
sum of C and D. The direct sum is a q-ary code that consists of all vectors of
the form (u,v), where u in C and v in D.}  
    require not IsNull(C): "Argument 1 can not be an empty code";
    require not IsNull(D): "Argument 2 can not be an empty code";
    require (C`BaseField cmpeq D`BaseField): "Arguments must have same base field";

    V := VectorSpace(C`BaseField, C`Length + D`Length);
    directSumRep:= [V!(Eltseq(c) cat Eltseq(d)) : c in CosetRepresentatives(C),
                                                  d in CosetRepresentatives(D)];                               
    directSumKernel := DirectSum(C`Kernel, D`Kernel);                      
    directSumCode := QaryCode(directSumKernel, directSumRep: IsFinalKernel := true);
    
    directSumCode`PAutSubgroup := DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
    directSumCode`IsPAutGroup := false;

    UpdateMinimumDistanceLowerBound(~directSumCode, 
          minPositive([C`MinimumDistanceLowerBound, D`MinimumDistanceLowerBound]));
    UpdateMinimumDistanceUpperBound(~directSumCode, 
          minPositive([C`MinimumDistanceUpperBound, D`MinimumDistanceUpperBound]));

    return directSumCode;
    
end intrinsic;

/*******************************************************************/
/*                                                                 */
/* Function name: DirectSum                                        */
/* Parameters: Q                                                   */
/* Function description: Given a sequence of q-ary codes Q = [C_1, */
/*   ..., C_r], all defined over the same base field, construct    */
/*   the direct sum of all these q-ary codes C_i, 1 <= i <= r. The */
/*   direct sum is a q-ary code that consists of all vectors of    */
/*   the form (u_1,u_2,...,u_r), where u_i in C_i.                 */
/* Input parameters description:                                   */
/*   - Q: A sequence of q-ary codes                                */
/* Output parameters description:                                  */
/*   - The direct sum code                                         */
/*                                                                 */
/* Signature: (<[CodeFld]> Q) -> CodeFld                           */
/*                                                                 */
/*******************************************************************/
intrinsic DirectSumBF(Q::[CodeFld]) -> CodeFld
{Given a sequence of q-ary codes Q = [C_1,..., C_r], all defined over the same
base field, construct the direct sum of all these q-ary codes C_i, 1 <= i <= r. 
The direct sum is a q-ary code that consists of all vectors of the form 
(u_1,u_2,...,u_r), where u_i in C_i.}
    r := #Q;
    require r gt 0: "The sequence of codes has to have at least one element";
    if r eq 1 then 
        return Q[1];
    end if;
    baseField := Q[1]`BaseField;
    for i in [2..r] do 
        require Q[i]`BaseField cmpeq baseField: "Codes must have the same base field";
    end for;

    codeLength := &+[C`Length : C in Q];
    V := VectorSpace(Q[1]`BaseField, codeLength);
    SetQ := <Set(C) : C in Q>;
    carProdArray := [Flat([Eltseq(xi) : xi in c]) : c in CartesianProduct(SetQ)];
    directSumCode := QaryCode([V!c : c in carProdArray]);
        
    upperBoundList := [C`MinimumDistanceUpperBound : C in Q];
    lowerBoundList := [C`MinimumDistanceLowerBound : C in Q];
    UpdateMinimumDistanceLowerBound(~directSumCode, minPositive(lowerBoundList));
    UpdateMinimumDistanceUpperBound(~directSumCode, minPositive(upperBoundList));
    
    directSumCode`PAutSubgroup := DirectProduct([C`PAutSubgroup : C in Q]);
    directSumCode`IsPAutGroup := false;

    return directSumCode;
    
end intrinsic;

intrinsic DirectSum(Q::[CodeFld]) -> CodeFld
{Given a sequence of q-ary codes Q = [C_1,..., C_r], all defined over the same
base field, construct the direct sum of all these q-ary codes C_i, 1 <= i <= r. 
The direct sum is a q-ary code that consists of all vectors of the form 
(u_1,u_2,...,u_r), where u_i in C_i.}
    r := #Q;
    require r gt 0: "The sequence of codes has to have at least one element";
    if r eq 1 then 
        return Q[1];
    end if;
    baseField := Q[1]`BaseField;
    for i in [2..r] do 
        require Q[i]`BaseField cmpeq baseField: "Codes must have the same base field";
    end for;

    codeLength := &+[C`Length : C in Q];
    V := VectorSpace(Q[1]`BaseField, codeLength);
    SetQ := <Set(CosetRepresentatives(C)) : C in Q>;
    carProdArray := [Flat([Eltseq(xi) : xi in c]) : c in CartesianProduct(SetQ)];
    directSumRep := [V!c : c in carProdArray];    
    directSumKernel := DirectSum([C`Kernel : C in Q]);                     
    directSumCode := QaryCode(directSumKernel, directSumRep: IsFinalKernel := true);

    upperBoundList := [C`MinimumDistanceUpperBound : C in Q];
    lowerBoundList := [C`MinimumDistanceLowerBound : C in Q];
    UpdateMinimumDistanceLowerBound(~directSumCode, minPositive(lowerBoundList));
    UpdateMinimumDistanceUpperBound(~directSumCode, minPositive(upperBoundList));
    
    directSumCode`PAutSubgroup := DirectProduct([C`PAutSubgroup : C in Q]);
    directSumCode`IsPAutGroup := false;

    return directSumCode;
    
end intrinsic;

intrinsic DirectSumI(Q::[CodeFld]) -> CodeFld
{Given a sequence of q-ary codes Q = [C_1,..., C_r], all defined over the same
base field, construct the direct sum of all these q-ary codes C_i, 1 <= i <= r. 
The direct sum is a q-ary code that consists of all vectors of the form 
(u_1,u_2,...,u_r), where u_i in C_i.}
    r := #Q;
    require r gt 0: "The sequence of codes has to have at least one element";
    if r eq 1 then 
        return Q[1];
    end if;
   
    directSumCode := Q[1];
    for i in [2..r] do
        directSumCode := DirectSum(directSumCode, Q[i]);
    end for;

    return directSumCode;
    
end intrinsic;

