/************************************************************/
/*                                                          */
/* Project name: Qary nonlinear codes in MAGMA              */
/* Test file name: QaryPadCode_BB_test.m                    */
/*                                                          */
/* Comments: Black-box tests for the function PadCode       */
/*           included in the QaryCodes_Constructions.m file */
/*                                                          */
/* Authors: Biplop Day and Mercè Villanueva                 */
/*                                                          */
/* Revision version and last date: 1.0, 2021/03/18          */
/*                                                          */
/************************************************************/

SetAssertions(true);
Alarm(30*60);

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
print "test 1: Zero code of length 8 over GF(2) and pad 3";
zeroCode := ZeroCode(GF(2), 8);
C := QaryCode(zeroCode);
n := C`Length;
pad := 3;

expectedOutput1 := ZeroQaryCode(GF(2), n+pad);
Output1 := PadCode(C, pad);
// DELETE Output2, since function BF won't be in the package!!!  
Output2 := PadCodeBF(C, pad);

assert expectedOutput1 eq Output1;
assert expectedOutput1 eq Output2;
assert Output1`MinimumDistance eq 0;
assert Output2`MinimumDistance eq 0;
assert Output1`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert not Output1`IsPAutGroup;
assert not Output2`IsPAutGroup;
assert Output1`IsLinear eq C`IsLinear;
assert Output2`IsLinear eq C`IsLinear;

/****************************************************************/
print "test 2: Universe code of lenght 10 over GF(3) and pad 2";
universeCode := UniverseCode(GF(3), 10);
C := QaryCode(universeCode);
n := C`Length;
pad := 2;

V := VectorSpace(GF(3), n+pad);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat [0,0]) : c in Set(C)]);
Output1 := PadCode(C, pad);
Output2 := PadCodeBF(C, pad);

assert expectedOutput1 eq Output1;
assert expectedOutput1 eq Output2;
assert Output1`MinimumDistance eq 1;
assert Output2`MinimumDistance eq 1;
assert Output1`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert not Output1`IsPAutGroup;
assert not Output2`IsPAutGroup;
assert Output1`IsLinear eq C`IsLinear;
assert Output2`IsLinear eq C`IsLinear;

/****************************************************************/
print "test 3: Repetition code of lenght 12 over GF(4) and pad 0";
repetitionCode := RepetitionCode(GF(4), 12);
C := QaryCode(repetitionCode);
n := C`Length;
pad := 0;

V := VectorSpace(GF(4), n+pad);
expectedOutput1 := QaryCode( [V!(Eltseq(c)) : c in Set(C)]);
Output1 := PadCode(C, pad);
Output2 := PadCodeBF(C, pad);

assert expectedOutput1 eq Output1;
assert expectedOutput1 eq Output2;
assert Output1`MinimumDistance eq 12;
assert Output2`MinimumDistance eq 12;
assert Output1`PAutSubgroup eq C`PAutSubgroup;
assert Output2`PAutSubgroup eq C`PAutSubgroup;
assert Output1`IsPAutGroup eq C`IsPAutGroup;
assert Output2`IsPAutGroup eq C`IsPAutGroup;
assert Output1`IsLinear eq C`IsLinear;
assert Output2`IsLinear eq C`IsLinear;

/****************************************************************/
print "test 4: q=4, n=4, kerDim=2, #C=16, #rep=1, linear, pad=5"; 
print "        d not assigned, 1<= d <= 3";
CHadq4n4ker2_seq := [ VectorSpace(GF(2, 2), 4) |
    [ GF(2, 2) | 0, 0, 0, 0 ] where a := GF(2, 2).1,
    [ GF(2, 2) | 1, 1, 1, 1 ] where a := GF(2, 2).1,
    [ GF(2, 2) | a, a, a, a ] where a := GF(2, 2).1,
    [ GF(2, 2) | a^2, a^2, a^2, a^2 ] where a := GF(2, 2).1,
    [ GF(2, 2) | 0, 1, a, a^2 ] where a := GF(2, 2).1,
    [ GF(2, 2) | 1, 0, a^2, a ] where a := GF(2, 2).1,
    [ GF(2, 2) | a, a^2, 0, 1 ] where a := GF(2, 2).1,
    [ GF(2, 2) | a^2, a, 1, 0 ] where a := GF(2, 2).1,
    [ GF(2, 2) | 0, a, a^2, 1 ] where a := GF(2, 2).1,
    [ GF(2, 2) | 1, a^2, a, 0 ] where a := GF(2, 2).1,
    [ GF(2, 2) | a, 0, 1, a^2 ] where a := GF(2, 2).1,
    [ GF(2, 2) | a^2, 1, 0, a ] where a := GF(2, 2).1,
    [ GF(2, 2) | 0, a^2, 1, a ] where a := GF(2, 2).1,
    [ GF(2, 2) | 1, a, 0, a^2 ] where a := GF(2, 2).1,
    [ GF(2, 2) | a, 1, a^2, 0 ] where a := GF(2, 2).1,
    [ GF(2, 2) | a^2, 0, a, 1 ] where a := GF(2, 2).1];
L := CHadq4n4ker2_seq;
C := QaryCode(L);
n := C`Length; 
pad := 5;

V := VectorSpace(GF(4), n+pad);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat [0^^pad]) : c in Set(C)]);
Output1 := PadCode(C, pad);
Output2 := PadCodeBF(C, pad);

assert expectedOutput1 eq Output1;
assert expectedOutput1 eq Output2;
assert Output1`MinimumDistanceLowerBound eq 1;
assert Output1`MinimumDistanceUpperBound eq 3;
assert Output1`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output1`IsPAutGroup eq C`IsPAutGroup;
assert Output2`IsPAutGroup eq C`IsPAutGroup;
assert Output1`IsLinear eq C`IsLinear;
assert Output2`IsLinear eq C`IsLinear;

/****************************************************************/
print "test 5: q=3, n=13, kerDim=9, #C=59049, #rep=3, nonlinear, pad=3";
print "        d not assigned, 1<= d <= 4"; 
kernel := LinearCode<GF(3), 13 |
                     [ GF(3) | 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2 ],
                     [ GF(3) | 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2 ],
                     [ GF(3) | 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1 ],
                     [ GF(3) | 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 2, 0 ],
                     [ GF(3) | 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2 ],
                     [ GF(3) | 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 2, 1 ],
                     [ GF(3) | 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 2 ],
                     [ GF(3) | 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 2, 1 ],
                     [ GF(3) | 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 ] >;
cosetRep := [ VectorSpace(GF(3), 13) | 
       Matrix(SparseMatrix(GF(3), 1, 13, \[ 0])),
       Matrix(GF(3), 1, 13, \[ 0, 2, 2, 0, 2, 2, 0, 2, 0, 0, 0, 0, 0 ]),
       Matrix(GF(3), 1, 13, \[ 2, 0, 2, 0, 2, 0, 0, 2, 2, 0, 0, 0, 0 ])];
C := QaryCode(kernel, cosetRep : IsFinalKernel := true);
n := C`Length; 
pad := 3;

V := VectorSpace(GF(3), n+pad);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat [0^^pad]) : c in Set(C)]);
Output1 := PadCode(C, pad);
Output2 := PadCodeBF(C, pad);

assert expectedOutput1 eq Output1;
assert expectedOutput1 eq Output2;
assert Output1`MinimumDistanceLowerBound eq 1;
assert Output1`MinimumDistanceUpperBound eq 4;
assert Output1`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output1`IsPAutGroup eq C`IsPAutGroup;
assert Output2`IsPAutGroup eq C`IsPAutGroup;
assert Output1`IsLinear eq C`IsLinear;
assert Output2`IsLinear eq C`IsLinear;

/****************************************************************/
print "test 6: q=8, n=20, kerDim=3, #C=2560, #rep=5, nonlinear, pad=6";
print "        d not assigned, 1<= d <= 17";
kernel := LinearCode<GF(2, 3), 20 |
     [ GF(2, 3) | 1, 0, 0, 0, 0, w^3, w, w^4, w^4, w^4, w^5, w^6, w, w^2, 
                  1, w^2, w^6, w, 1, 0 ] where w := GF(2, 3).1,
     [ GF(2, 3) | 0, 1, 0, 0, 0, 1, w^5, w^3, 1, w^2, w^3, w^5, w^6, w^3, 
                  w, 0, w, w^3, w^5, w^3 ] where w := GF(2, 3).1,
     [ GF(2, 3) | 0, 0, 1, 0, 0, w^5, w^5, w^5, w^6, 0, w^5, w^2, w^6, w^5, 
                  w^3, w^3, 0, w^4, 0, 0 ] where w := GF(2, 3).1 >;
cosetRep := [ VectorSpace(GF(2, 3), 20) |
     [ GF(2, 3) | 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ] 
                  where W := GF(2, 3).1,
     [ GF(2, 3) | W, W^3, W^4, W^2, W, 0, 1, W, W^4, 0, W, W^2, W, W^4, W^2, 0, 
                  W^4, W^4, W^4, W^4 ] where W := GF(2, 3).1,
     [ GF(2, 3) | W^6, W^3, W^3, W, W^5, 0, W^5, W, W, 0, W^4, W^4, 1, W^4, W, W, 
                  W^4, W^4, 0, W^2 ] where W := GF(2, 3).1,
     [ GF(2, 3) | W^3, W^6, W^2, W, 1, W^6, W^2, W^4, 0, 0, W^4, W^4, W^2, W, W^3, 
                  W^2, W^2, W^4, W, W^2 ] where W := GF(2, 3).1,
     [ GF(2, 3) | W^2, W^6, 1, W^2, W^2, W^6, W, W^4, W^2, 0, W, W^2, W^5, W, W^6, 
                  W^4, W^2, W^4, W^2, W^4 ] where W := GF(2, 3).1] ;
C := QaryCode(kernel, cosetRep : IsFinalKernel := true);
n := C`Length; 
pad := 6;

V := VectorSpace(GF(8), n+pad);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat [0^^pad]) : c in Set(C)]);
Output1 := PadCode(C, pad);
Output2 := PadCodeBF(C, pad);

assert expectedOutput1 eq Output1;
assert expectedOutput1 eq Output2;
assert Output1`MinimumDistanceLowerBound eq 1;
assert Output1`MinimumDistanceUpperBound eq 17;
assert Output1`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output1`IsPAutGroup eq C`IsPAutGroup;
assert Output2`IsPAutGroup eq C`IsPAutGroup;
assert Output1`IsLinear eq C`IsLinear;
assert Output2`IsLinear eq C`IsLinear;

/****************************************************************/
print "test 7: q=3, n=13, kerDim=8, #C=59049, #rep=9, nonlinear, pad=3";
print "        d not assigned, 1<= d <= 4";
load "./data/Cq3n13ker8c.m";
load "./data/Cq3n13ker8c_seqgen.m";
L := Cq3n13ker8c_seq;
kernel := LinearCode(Matrix(Cq3n13ker8c`Kernel));
C := QaryCode(kernel, Cq3n13ker8c`CosetRepresentatives : IsFinalKernel := true);
n := C`Length; 
pad := 3;

V := VectorSpace(GF(3), n+pad);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat [0^^pad]) : c in Set(C)]);
Output1 := PadCode(C, pad);
Output2 := PadCodeBF(C, pad);

assert expectedOutput1 eq Output1;
assert expectedOutput1 eq Output2;
assert Output1`MinimumDistanceLowerBound eq 1;
assert Output1`MinimumDistanceUpperBound eq 4;
assert Output1`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output1`IsPAutGroup eq C`IsPAutGroup;
assert Output2`IsPAutGroup eq C`IsPAutGroup;
assert Output1`IsLinear eq C`IsLinear;
assert Output2`IsLinear eq C`IsLinear;

/****************************************************************/
print "test 8: q=2, n=20, kerDim=0, #C=127, #rep=127, nonlinear, 0 notin C, pad=2";
print "        d not assigned, 1<= d <= 10";
load "./data/Cq2n20ker0.m";
load "./data/Cq2n20ker0_seqgen.m";
L := Cq2n20ker0_seq;
kernel := LinearCode(Matrix(Cq2n20ker0`Kernel));
C := QaryCode(kernel, Cq2n20ker0`CosetRepresentatives : IsFinalKernel := true);
C`MinimumDistanceUpperBound := 10;
n := C`Length; 
pad := 2;

V := VectorSpace(GF(2), n+pad);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat [0^^pad]) : c in Set(C)]);
Output1 := PadCode(C, pad);
Output2 := PadCodeBF(C, pad);

assert expectedOutput1 eq Output1;
assert expectedOutput1 eq Output2;
assert Output1`MinimumDistanceLowerBound eq 1;
assert Output1`MinimumDistanceUpperBound eq 10;
assert Output1`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output1`IsPAutGroup eq C`IsPAutGroup;
assert Output2`IsPAutGroup eq C`IsPAutGroup;
assert Output1`IsLinear eq C`IsLinear;
assert Output2`IsLinear eq C`IsLinear;

/****************************************************************/
print "test 9: q=2, n=16, kerDim=4, #C=2048, #rep=128, nonlinear, 0 notin C, pad=16";
print "        d not assigned, 1<= d <= 6";
load "./data/CHamq2n16ker4t.m";
load "./data/CHamq2n16ker4t_seqgen.m";
L := CHamq2n16ker4t_seq;
kernel := LinearCode(Matrix(CHamq2n16ker4t`Kernel));
C := QaryCode(kernel, CHamq2n16ker4t`CosetRepresentatives : IsFinalKernel := true);
n := C`Length; 
pad := 16;

V := VectorSpace(GF(2), n+pad);
expectedOutput1 := QaryCode([ V!(Eltseq(c) cat [0^^pad]) : c in Set(C)]);
Output1 := PadCode(C, pad);
Output2 := PadCodeBF(C, pad);

assert expectedOutput1 eq Output1;
assert expectedOutput1 eq Output2;
assert Output1`MinimumDistanceLowerBound eq 1;
assert Output1`MinimumDistanceUpperBound eq 6;
assert Output1`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, Sym(pad));
assert Output1`IsPAutGroup eq C`IsPAutGroup;
assert Output2`IsPAutGroup eq C`IsPAutGroup;
assert Output1`IsLinear eq C`IsLinear;
assert Output2`IsLinear eq C`IsLinear;









