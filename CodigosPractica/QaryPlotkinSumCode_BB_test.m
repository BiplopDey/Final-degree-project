/************************************************************/
/*                                                          */
/* Project name: Qary nonlinear codes in MAGMA              */
/* Test file name: QaryPlotkinSumCode_BB_test.m             */
/*                                                          */
/* Comments: Black-box tests for the function PlotkinSum    */
/*           included in the QaryCodes_Constructions.m file */
/*                                                          */
/* Authors: Biplop Dey and Mercè Villanueva                 */
/*                                                          */
/* Revision version and last date: 1.0, 2021/03/18          */
/*                                 1.1, 2021/06/03          */
/*                                                          */
/************************************************************/

SetAssertions(true);
Alarm(30*60);

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
print "test 1: Two zero codes of length 8 over GF(2)";
zeroCode := ZeroCode(GF(2), 8);
C := QaryCode(zeroCode);
D := C;

n1 := C`Length;
n2 := D`Length;
expectedOutput1 := ZeroQaryCode(GF(2), n1+n2);
OutputKC := PlotkinSum(C, D);
OutputBF := PlotkinSumBF(C, D);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq OutputBF;
assert OutputKC`MinimumDistance eq 0;
assert OutputBF`MinimumDistance eq 0;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert OutputBF`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert not OutputBF`IsPAutGroup;
assert OutputKC`IsLinear;
assert OutputBF`IsLinear;

/****************************************************************/
print "test 2: Two universe codes of lenght 5 over GF(3)";
universeCode := UniverseCode(GF(3), 5);
C := QaryCode(universeCode);
D := C;

n1 := C`Length;
n2 := D`Length;
V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(c+d)) : c in Set(C), d in Set(D)]);
OutputKC := PlotkinSum(C, D);
OutputBF := PlotkinSumBF(C, D);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq OutputBF;
assert OutputKC`MinimumDistance eq 1;
assert OutputBF`MinimumDistance eq 1;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert OutputBF`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert OutputKC`IsLinear;
assert OutputBF`IsLinear;

/****************************************************************/
print "test 3: Universe and zero code of lenght 8 over GF(3)";
universeCode := UniverseCode(GF(3), 8);
C := QaryCode(universeCode);
D := QaryCode(ZeroCode(GF(3), 8));

n1 := C`Length;
n2 := D`Length;
V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(c+d)) : c in Set(C), d in Set(D)]);
OutputKC := PlotkinSum(C, D);
OutputBF := PlotkinSumBF(C, D);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq OutputBF;
assert OutputKC`MinimumDistance eq 2;
assert OutputBF`MinimumDistance eq 2;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert OutputBF`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert OutputKC`IsLinear;
assert OutputBF`IsLinear;

/****************************************************************/
print "test 4: q=2, Codes with kernel of dimension 0";
print "Code C: q=2, n=20, kerDim=0, MinDistance=4, notLinear, #C=127, #costeRep=127, 0 notin C, d not assigned, 1<=d<=14";
load "./data/Cq2n20ker0.m";
load "./data/Cq2n20ker0_seqgen.m";
L := Cq2n20ker0;
kernel := LinearCode(Matrix(L`Kernel));
C := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);

print "Code D: q=2, n=20, kerDim=0, MinDistance=1, notLinear, #D=100, #costeRep=100, 0 notin D, d not assigned, 1<=d<=4";
load "./data/Cq2n10ker0.m";
load "./data/Cq2n10ker0_seqgen.m";
L := Cq2n10ker0;
kernel := LinearCode(Matrix(L`Kernel));
D := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);
D := PadCode(D, 10);

n1 := C`Length;
n2 := D`Length;
V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(c+d)) : c in Set(C), d in Set(D)]);
OutputKC := PlotkinSum(C, D);
OutputBF := PlotkinSumBF(C, D);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq OutputBF;
assert OutputKC`MinimumDistanceLowerBound eq 1;
assert OutputBF`MinimumDistanceLowerBound eq 1;
assert OutputKC`MinimumDistanceUpperBound eq 4;
assert OutputBF`MinimumDistanceUpperBound eq 4;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert OutputBF`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert not OutputKC`IsLinear;
assert not OutputBF`IsLinear;

/****************************************************************/
print "test 5: q=2, Same code with kernel of dimension 0";
print "Code C: q=2, n=10, kerDim=0, MinDistance=1, notLinear, #C=100, #costeRep=100, 0 in C";
load "./data/Cq2n10ker0.m";
load "./data/Cq2n10ker0_seqgen.m";
L := Cq2n10ker0;
kernel := LinearCode(Matrix(L`Kernel));
C := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);
D := C;

d1 := MinimumDistance(C);

n1 := C`Length;
n2 := n1;
V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(c+d)): c in Set(C), d in Set(C)]);
OutputKC := PlotkinSum(C, D);
OutputBF := PlotkinSumBF(C, D);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq OutputBF;
assert OutputKC`MinimumDistance eq 1;
assert OutputBF`MinimumDistance eq 1;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert OutputBF`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert not OutputKC`IsLinear;
assert not OutputBF`IsLinear;

/****************************************************************/
print "test 6: q=4, codes with kernel of dimension 2 and 3";
print "Code C: q=4, n=16, kerDim=2, MinDistance=12, notLinear, #C=64, #costeRep=4, 0 in C, d not assigned, 2<=d<=14";
load "./data/CHadq4n16ker2a.m";
load "./data/CHadq4n16ker2a_seqgen.m";
L := CHadq4n16ker2a;
kernel := LinearCode(Matrix(L`Kernel));
C := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);
C`MinimumDistanceLowerBound := 2;

print "Code D: q=4,n=16,kerDim=3,MinDistance=12,notLinear, #D=64, #costeRep=1, 0 notin D ";
load "./data/CHadq4n16ker3t.m";
load "./data/CHadq4n16ker3t_seqgen.m";
L := CHadq4n16ker3t;
kernel := LinearCode(Matrix(L`Kernel));
D := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);

d2 := MinimumDistance(D);

n1 := C`Length;
n2 := D`Length;
V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(c+d)): c in Set(C), d in Set(D)]);
OutputKC := PlotkinSum(C, D);
OutputBF := PlotkinSumBF(C, D);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq OutputBF;
assert OutputKC`MinimumDistanceLowerBound eq 4;
assert OutputBF`MinimumDistanceLowerBound eq 4;
assert OutputKC`MinimumDistanceUpperBound eq 12;
assert OutputBF`MinimumDistanceUpperBound eq 12;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert OutputBF`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert not OutputKC`IsLinear;
assert not OutputBF`IsLinear;

/****************************************************************/
print "test 7: q=4, Codes with kernel of dimension 1 and 0";
print "Code C: q=4, n=16, kerDim=1, MinDistance= 12, notLinear, #C=64, #costeRep=16, 0 in C, d not assigned, 1<=d<=14";
load "./data/CHadq4n16ker1b.m";
load "./data/CHadq4n16ker1b_seqgen.m";
L := CHadq4n16ker1b;
kernel := LinearCode(Matrix(L`Kernel));
C := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);

print "Code D: q=4, n=16, kerDim=0, MinDistance= 6, notLinear, #C=16, #costeRep=16, 0 in C, d not assigned, 1<=d<=15";
load "./data/Cq4n16ker0.m";
load "./data/Cq4n16ker0_seqgen.m";
L := Cq4n16ker0;
kernel := LinearCode(Matrix(L`Kernel));
D := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);

n1 := C`Length;
n2 := D`Length;
V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(c+d)): c in Set(C), d in Set(D)]);
OutputKC := PlotkinSum(C, D);
OutputBF := PlotkinSumBF(C, D);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq OutputBF;
assert OutputKC`MinimumDistanceLowerBound eq 1;
assert OutputBF`MinimumDistanceLowerBound eq 1;
assert OutputKC`MinimumDistanceUpperBound eq 15;
assert OutputBF`MinimumDistanceUpperBound eq 15;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert OutputBF`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert not OutputKC`IsLinear;
assert not OutputBF`IsLinear;

/****************************************************************/
print "test 8: q=4, Codes with kernel of dimension 1 and 2";
print "Code C: q=4, n=16, kerDim=1, MinDistance=12, notLinear, #C=64, #costeRep=16, 0 in C, d not assigned, 1<=d<=13";
load "./data/CHadq4n16ker1c.m";
load "./data/CHadq4n16ker1c_seqgen.m";
L := CHadq4n16ker1c;
kernel := LinearCode(Matrix(L`Kernel));
C := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);
C`MinimumDistanceUpperBound := 13;

print "Code D: q=4, n=16, kerDim=2, MinDistance=1, notLinear, #D=32, #costeRep=2, 0 in D, d not assigned, 1<=d<=2";
load "./data/Cq4n4ker2.m";
load "./data/Cq4n4ker2_seqgen.m";
L := Cq4n4ker2;
kernel := LinearCode(Matrix(L`Kernel));
D := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);
D := PadCode(D, 12);

n1 := C`Length;
n2 := D`Length;
V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(c+d)): c in Set(C), d in Set(D)]);
OutputKC := PlotkinSum(C, D);
OutputBF := PlotkinSumBF(C, D);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq OutputBF;

assert OutputKC`MinimumDistanceLowerBound eq 1;
assert OutputBF`MinimumDistanceLowerBound eq 1;
assert OutputKC`MinimumDistanceUpperBound eq 2;
assert OutputBF`MinimumDistanceUpperBound eq 2;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert OutputBF`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert not OutputKC`IsLinear;
assert not OutputBF`IsLinear;

/****************************************************************/
print "test 9: Repetition code of lenght 20 over GF(8) and code with dim kernel 3";
repetitionCode := RepetitionCode(GF(8), 20); 
C := QaryCode(repetitionCode);

print "Code D: q=8, n=20, kerDim=3, lower=1, upper=17, notLinear, #D=2560, #costeRep=5, 0 in D, d not assigned, 1<=d<=17";
load "./data/Cq8n20ker3.m";
load "./data/Cq8n20ker3_seqgen.m";
L := Cq8n20ker3;
kernel := LinearCode(Matrix(L`Kernel));
D := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);

n1 := C`Length;
n2 := D`Length;
V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(c+d)): c in Set(C), d in Set(D)]);
OutputKC := PlotkinSum(C, D);
OutputBF := PlotkinSumBF(C, D);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq OutputBF;
assert OutputKC`MinimumDistanceLowerBound eq 1;
assert OutputBF`MinimumDistanceLowerBound eq 1;
assert OutputKC`MinimumDistanceUpperBound eq 17;
assert OutputBF`MinimumDistanceUpperBound eq 17;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert OutputBF`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert not OutputKC`IsLinear;
assert not OutputBF`IsLinear;

/****************************************************************/
print "test 10: linear codes with length 7 over GF(3)";
M := KMatrixSpace(GF(3), 4, 7);
G := M![0,0,0,2,0,1,1,1,2,0,0,0,0,0,1,1,2,0,0,1,1,0,2,1,0,1,0,0];
print "Code C: q=3, n=7, kerDim=4, Linear, #C=81, d not assigned, 1<=d<=4";
C := QaryCode(LinearCode(G));

M := KMatrixSpace(GF(3), 3, 7);
G := M![1,1,2,0,3,1,1,0,0,1,0,1,0,2,1,1,1,0,2,0,0];
print "Code D: q=3, n=7, kerDim=3, Linear, #D=27, d not assigned, 1<=d<=5";
D := QaryCode(LinearCode(G));

n1 := C`Length;
n2 := D`Length;
V := VectorSpace(GF(3),n1+n2 );
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(c+d)): c in Set(C), d in Set(D)]);
OutputKC := PlotkinSum(C, D);
OutputBF := PlotkinSumBF(C, D);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq OutputBF;
assert OutputKC`MinimumDistanceLowerBound eq 1;
assert OutputBF`MinimumDistanceLowerBound eq 1;
assert OutputKC`MinimumDistanceUpperBound eq 5;
assert OutputBF`MinimumDistanceUpperBound eq 5;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert OutputBF`PAutSubgroup eq DirectProduct(C`PAutSubgroup, 
                                C`PAutSubgroup meet D`PAutSubgroup);
assert OutputKC`IsLinear;
assert OutputBF`IsLinear;