/************************************************************/
/*                                                          */
/* Project name: Qary nonlinear codes in MAGMA              */
/* Test file name: QaryDirectSumQCode_BB_test.m             */
/*                                                          */
/* Comments: Black-box tests for the function DirectCode    */
/*           included in the QaryCodes_Constructions.m file */
/*                                                          */
/* Authors: Biplop Dey and Mercè Villanueva                 */
/*                                                          */
/* Revision version and last date: 1.0, 2021/03/18          */
/*                                 1.1, 2021/06/21          */
/*                                                          */
/************************************************************/

SetAssertions(true);
Alarm(30*60);

/****************************************************************/
print "test 1: Two zero codes of length 8 over GF(2)";
zeroCode := ZeroCode(GF(2), 8);
C := QaryCode(zeroCode);
D := C;
n1 := C`Length;
n2 := D`Length;

expectedOutput1 := ZeroQaryCode(GF(2), n1+n2);
OutputKC := DirectSum([C, D]);
Output2 := DirectSumI([C, D]);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistance eq 0;
assert Output2`MinimumDistance eq 0;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert OutputKC`IsLinear;
assert Output2`IsLinear;

/****************************************************************/
print "test 2: Two Universe codes of lenght 5 over GF(3)";
universeCode := UniverseCode(GF(3), 5);
C := QaryCode(universeCode);
D := C;
n1 := C`Length;
n2 := D`Length;

V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(d)): c in Set(C), d in Set(D)]);
OutputKC := DirectSum([C, D]);
Output2 := DirectSumI([C, D]);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistance eq 1;
assert Output2`MinimumDistance eq 1;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert OutputKC`IsLinear;
assert Output2`IsLinear;

/****************************************************************/
print "test 3: Universe code of lenght 8 over GF(3) and zero code of length 8 over GF(3)";
universeCode := UniverseCode(GF(3), 10);
C := QaryCode(universeCode);
D := QaryCode(ZeroCode(GF(3), 8));
n1 := C`Length;
n2 := D`Length;

V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(d)): c in Set(C), d in Set(D)]);
OutputKC := DirectSum([C, D]);
Output2 := DirectSumI([C, D]);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistance eq 1;
assert Output2`MinimumDistance eq 1;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert OutputKC`IsLinear;
assert Output2`IsLinear;

/****************************************************************/
print "test 4: q=2, Codes with kernel dim 0";
print "Code C: q=2, n=20, kerDim=0, MinDistance=4, notLinear, #C=127, #costeRep=127, 0 notin C, d not assigned, 1<= d <=14";
load "./data/Cq2n20ker0.m";
load "./data/Cq2n20ker0_seqgen.m";
L := Cq2n20ker0;
kernel := LinearCode(Matrix(L`Kernel));
C := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);

print "Code D: q=2, n=10, kerDim=0, MinDistance=1, notLinear, #D=100, #costeRep=100, 0 notin D, d not assigned, 1<= d <=4";
load "./data/Cq2n10ker0.m";
load "./data/Cq2n10ker0_seqgen.m";
L := Cq2n10ker0;
kernel := LinearCode(Matrix(L`Kernel));
D := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);

n1 := C`Length;
n2 := D`Length;
V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(d)): c in Set(C), d in Set(D)]);
OutputKC := DirectSum([C, D]);
Output2 := DirectSumI([C, D]);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistanceLowerBound eq 1;
assert Output2`MinimumDistanceLowerBound eq 1;
assert OutputKC`MinimumDistanceUpperBound eq 4;
assert Output2`MinimumDistanceUpperBound eq 4;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert not OutputKC`IsLinear;
assert not Output2`IsLinear;

/****************************************************************/
print "test 5: q=2, Same code";
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
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(d)): c in Set(C), d in Set(C)]);
OutputKC := DirectSum([C, D]);
Output2 := DirectSumI([C, D]);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistance eq 1;
assert Output2`MinimumDistance eq 1;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert not OutputKC`IsLinear;
assert not Output2`IsLinear;

/****************************************************************/
print "test 6: q=4, Codes with ker dim 2 and 3";
print "Code C: q=4, n=16, kerDim=2, MinDistance=12, notLinear, #C=64, #costeRep=4, 0 in C, d not assigned, 2<=d<=14";
load "./data/CHadq4n16ker2a.m";
load "./data/CHadq4n16ker2a_seqgen.m";
L := CHadq4n16ker2a;
kernel := LinearCode(Matrix(L`Kernel));
C := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);
C`MinimumDistanceLowerBound:=2;

print "Code D: q=4,n=16,kerDim=3,MinDistance=12,notLinear, #D=64, #costeRep=1, 0 notin D ";
load "./data/CHadq4n16ker3t.m";
load "./data/CHadq4n16ker3t_seqgen.m";
L := CHadq4n16ker3t;
kernel := LinearCode(Matrix(L`Kernel));
D := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);
d2 :=MinimumDistance(D);

n1 := C`Length;
n2 := D`Length;
V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(d)): c in Set(C), d in Set(D)]);
OutputKC := DirectSum([C, D]);
Output2 := DirectSumI([C, D]);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistanceLowerBound eq 2;
assert Output2`MinimumDistanceLowerBound eq 2;
assert OutputKC`MinimumDistanceUpperBound eq 12;
assert Output2`MinimumDistanceUpperBound eq 12;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert not OutputKC`IsLinear;
assert not Output2`IsLinear;

/****************************************************************/
print "test 7: q=4, Codes with ker dim 1 and 0";
print "Code C: q=4, n=16, kerDim=1, MinDistance=12, notLinear, #C=64, #costeRep=16, 0 in C, d not assigned, 1<=d<=14";
load "./data/CHadq4n16ker1b.m";
load "./data/CHadq4n16ker1b_seqgen.m";
L := CHadq4n16ker1b;
kernel := LinearCode(Matrix(L`Kernel));
C := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);

print "Code D: q=4, n=16, kerDim=0, MinDistance=6, notLinear, #C=16, #costeRep=16, 0 in C, d not assigned, 1<=d<=15";
load "./data/Cq4n16ker0.m";
load "./data/Cq4n16ker0_seqgen.m";
L := Cq4n16ker0;
kernel := LinearCode(Matrix(L`Kernel));
D := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);

n1 := C`Length;
n2 := D`Length;
V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(d)): c in Set(C), d in Set(D)]);
OutputKC := DirectSum([C, D]);
Output2 := DirectSumI([C, D]);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistanceLowerBound eq 1;
assert Output2`MinimumDistanceLowerBound eq 1;
assert OutputKC`MinimumDistanceUpperBound eq 14;
assert Output2`MinimumDistanceUpperBound eq 14;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert not OutputKC`IsLinear;
assert not Output2`IsLinear;

/****************************************************************/
print "test 8: q=4, Codes with ker dim 1 and 2";
print "Code C: q=4, n=16, kerDim=1, MinDistance=12, notLinear, #C=64, #costeRep=16, 0 in C, d not assigned, 1<=d<=13";
load "./data/CHadq4n16ker1c.m";
load "./data/CHadq4n16ker1c_seqgen.m";
L := CHadq4n16ker1c;
kernel := LinearCode(Matrix(L`Kernel));
C := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);
C`MinimumDistanceUpperBound:=13;

print "Code D: q=4, n=4, kerDim=2, MinDistance= 1, notLinear, #D=32, #costeRep=2, 0 in D, d not assigned, 1<=d<=2";
load "./data/Cq4n4ker2.m";
load "./data/Cq4n4ker2_seqgen.m";
L := Cq4n4ker2;
kernel := LinearCode(Matrix(L`Kernel));
D := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);

n1 := C`Length;
n2 := D`Length;
V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(d)): c in Set(C), d in Set(D)]);
OutputKC := DirectSum([C, D]);
Output2 := DirectSumI([C, D]);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistanceLowerBound eq 1;
assert Output2`MinimumDistanceLowerBound eq 1;
assert OutputKC`MinimumDistanceUpperBound eq 2;
assert Output2`MinimumDistanceUpperBound eq 2;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert not OutputKC`IsLinear;
assert not Output2`IsLinear;

/****************************************************************/
print "test 9: Repetition code of lenght 20 over GF(8) and code with dim kernel 3";
repetitionCode := RepetitionCode(GF(8), 3); 
C := QaryCode(repetitionCode);
d1 := C`MinimumDistance;

print "Code D: q=8, n=20, kerDim=3, dlower=1, dupper=17, notLinear, #D=2560, #costeRep=5, 0 in D, d not assigned, 1<=d<=17";
load "./data/Cq8n20ker3.m";
load "./data/Cq8n20ker3_seqgen.m";
L := Cq8n20ker3;
kernel := LinearCode(Matrix(L`Kernel));
D := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);

n1 := C`Length;
n2 := D`Length;
V := VectorSpace(C`BaseField, n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(d)): c in Set(C), d in Set(D)]);
OutputKC := DirectSum([C, D]);
Output2 := DirectSumI([C, D]);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistanceLowerBound eq 1;
assert Output2`MinimumDistanceLowerBound eq 1;
assert OutputKC`MinimumDistanceUpperBound eq 3;
assert Output2`MinimumDistanceUpperBound eq 3;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert not OutputKC`IsLinear;
assert not Output2`IsLinear;

/****************************************************************/
print "test 10: Linear codes of length 7 over GF(3)";
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
V := VectorSpace(GF(3), n1+n2);
expectedOutput1 := QaryCode([V!(Eltseq(c) cat Eltseq(d)): c in Set(C), d in Set(D)]);
OutputKC := DirectSum([C, D]);
Output2 := DirectSumI([C, D]);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistanceLowerBound eq 1;
assert Output2`MinimumDistanceLowerBound eq 1;
assert OutputKC`MinimumDistanceUpperBound eq 4;
assert Output2`MinimumDistanceUpperBound eq 4;
assert OutputKC`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert Output2`PAutSubgroup eq DirectProduct(C`PAutSubgroup, D`PAutSubgroup);
assert OutputKC`IsLinear;
assert Output2`IsLinear;

/****************************************************************/
print "test 11: Three zero codes over GF(5) of length 3, 5 and 1";
C1 := QaryCode(ZeroCode(GF(5), 3));
C2 := QaryCode(ZeroCode(GF(5), 5));
C3 := QaryCode(ZeroCode(GF(5), 1));
Q := [C1, C2, C3];

V := VectorSpace(C1`BaseField, 3+5+1);
expectedOutput1 := QaryCode([V!(Eltseq(c1) cat Eltseq(c2) cat Eltseq(c3)): 
                                c1 in Set(C1), c2 in Set(C2), c3 in Set(C3)]);
OutputKC := DirectSum(Q);
Output2 := DirectSumI(Q);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistanceLowerBound eq 0;
assert Output2`MinimumDistanceLowerBound eq 0;
subGroup := DirectProduct([C`PAutSubgroup : C in Q]);
assert OutputKC`PAutSubgroup eq subGroup;
assert Output2`PAutSubgroup eq subGroup;
assert OutputKC`IsLinear;
assert Output2`IsLinear;

/****************************************************************/
print "test 12: Several universal codes";
C1 := QaryCode(UniverseCode(GF(3), 3));
C2 := QaryCode(UniverseCode(GF(3), 2));
C3 := QaryCode(UniverseCode(GF(3), 1));
Q := [C1, C2, C3];

V := VectorSpace(C1`BaseField, 3+2+1);
expectedOutput1 := QaryCode([V!(Eltseq(c1) cat Eltseq(c2) cat Eltseq(c3)): 
                            c1 in Set(C1), c2 in Set(C2), c3 in Set(C3)]);
OutputKC := DirectSum(Q);
Output2 := DirectSumI(Q);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistanceLowerBound eq 1;
assert Output2`MinimumDistanceLowerBound eq 1;
subGroup := DirectProduct([C`PAutSubgroup : C in Q]);
assert OutputKC`PAutSubgroup eq subGroup;
assert Output2`PAutSubgroup eq subGroup;
assert OutputKC`IsLinear;
assert Output2`IsLinear;

/****************************************************************/
print "test 13: Three universal codes over GF(3) of length 2, 1, 3 and three zero codes over GF(3) of length 5, 3, 1";
C1 := QaryCode(UniverseCode(GF(3), 2));
C2 := QaryCode(ZeroCode(GF(3), 5));
C3 := QaryCode(UniverseCode(GF(3), 1));
C4 := QaryCode(ZeroCode(GF(3), 3));
C5 := QaryCode(UniverseCode(GF(3), 3));
C6 := QaryCode(ZeroCode(GF(3), 1));
Q := [C1, C2, C3, C4, C5, C6];

V := VectorSpace(C1`BaseField, 3+2+1+3+5+1);
expectedOutput1 := QaryCode([V!(Eltseq(c1) cat Eltseq(c2) cat Eltseq(c3) cat Eltseq(c4) cat Eltseq(c5) cat Eltseq(c6)): 
                            c1 in Set(C1), c2 in Set(C2), c3 in Set(C3), c4 in Set(C4), c5 in Set(C5), c6 in Set(C6)]);
OutputKC := DirectSum(Q);
Output2 := DirectSumI(Q);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistanceLowerBound eq 1;
assert Output2`MinimumDistanceLowerBound eq 1;
subGroup := DirectProduct([C`PAutSubgroup : C in Q]);
assert OutputKC`PAutSubgroup eq subGroup;
assert Output2`PAutSubgroup eq subGroup;
assert OutputKC`IsLinear;
assert Output2`IsLinear;

/****************************************************************/
print "test 14: Several codes";
print "Code C1: q=4, n=16, kerDim=0, MinDistance= 6, notLinear, #C=16, #costeRep=16, 0 in C, d not assigned, 1<=d<=15";
load "./data/Cq4n16ker0.m";
load "./data/Cq4n16ker0_seqgen.m";
L := Cq4n16ker0;
kernel := LinearCode(Matrix(L`Kernel));
C1 := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);

print "Code C2: q=4, n=4, kerDim=2, MinDistance= 1, notLinear, #C=32, #costeRep=2, 0 in C, d not assigned, 1<=d<=2";
load "./data/Cq4n4ker2.m";
load "./data/Cq4n4ker2_seqgen.m";
L := Cq4n4ker2;
kernel := LinearCode(Matrix(L`Kernel));
C2 := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);

print "Code C3: q=4, n=16, kerDim=1, MinDistance= 12, notLinear, #C=64, #costeRep=16, 0 in C, d not assigned, 1<=d<=13";
load "./data/CHadq4n16ker1c.m";
load "./data/CHadq4n16ker1c_seqgen.m";
L := CHadq4n16ker1c;
kernel := LinearCode(Matrix(L`Kernel));
C3 := QaryCode(kernel, L`CosetRepresentatives : IsFinalKernel := true);
C3`MinimumDistanceUpperBound := 13;

V := VectorSpace(C1`BaseField, 16+16+4);
expectedOutput1 := QaryCode([V!(Eltseq(c1) cat Eltseq(c2) cat Eltseq(c3)): 
                            c1 in Set(C1), c2 in Set(C2), c3 in Set(C3)]);
Q := [C1, C2, C3];
OutputKC := DirectSum(Q);
Output2 := DirectSumI(Q);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistanceLowerBound eq 1;
assert Output2`MinimumDistanceLowerBound eq 1;
assert OutputKC`MinimumDistanceUpperBound eq 2;
assert Output2`MinimumDistanceUpperBound eq 2;
subGroup := DirectProduct([C`PAutSubgroup : C in Q]);
assert OutputKC`PAutSubgroup eq subGroup;
assert Output2`PAutSubgroup eq subGroup;
assert not OutputKC`IsLinear;
assert not Output2`IsLinear;

/****************************************************************/
print "test 15: Previous codes C1, C2, C3 and zero code of length 5 over GF(4)";
Z := QaryCode(ZeroCode(GF(4), 5));

V := VectorSpace(C1`BaseField, 16+16+4+5);
expectedOutput1 := QaryCode([V!(Eltseq(c1) cat Eltseq(c2) cat Eltseq(z) cat Eltseq(c3)): 
                            c1 in Set(C1), c2 in Set(C2), z in Set(Z), c3 in Set(C3)]);
Q := [C1, C2, Z, C3];
OutputKC := DirectSum(Q);
Output2 := DirectSumI(Q);

assert expectedOutput1 eq OutputKC;
assert expectedOutput1 eq Output2;
assert OutputKC`MinimumDistanceLowerBound eq 1;
assert Output2`MinimumDistanceLowerBound eq 1;
assert OutputKC`MinimumDistanceUpperBound eq 2;
assert Output2`MinimumDistanceUpperBound eq 2;
subGroup := DirectProduct([C`PAutSubgroup : C in Q]);
assert OutputKC`PAutSubgroup eq subGroup;
assert Output2`PAutSubgroup eq subGroup;
assert not OutputKC`IsLinear;
assert not Output2`IsLinear;