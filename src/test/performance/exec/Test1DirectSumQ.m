/****************************************************************/
/* Project name: Qary nonlinear codes in MAGMA    Test  4sd     */
/****************************************************************/

load "TestPerformanceDirectSumQ.m";
load "TestPerformanceGeneral.m";
methods := ["Iterative", "Direct"];
fileOutput := createFile("TestDirectSumQIterativeVsDirect");

load "./data/Cq4n16ker0.m";//       q=4, n=16, kerDim=0, #C=16
load "./data/CHadq4n4ker2.m";//     q=4, n=4,  kerDim=2, #C=16
load "./data/CHadq4n8ker1.m";//     q=4, n=8,  kerDim=1, #C=32
load "./data/Cq4n4ker2.m";//        q=4, n=16, kerDim=2, #C=32
load "./data/CHadq4n12ker1.m";//    q=4, n=12, kerDim=1, #C=48
load "./data/CHadq4n16ker1c.m";//   q=4, n=16, kerDim=1, #C=64
load "./data/CHadq4n16ker1b.m";//   q=4, n=16, kerDim=1, #C=64
load "./data/CHadq4n16ker2a.m";//   q=4, n=16, kerDim=2, #C=64
//load "./data/CHadq4n16ker3t.m";//   q=4, n=16, kerDim=3, #C=64
load "./data/CHadq4n16ker3.m";//    q=4, n=16, kerDim=3, #C=64

load "./data/Cq2n10ker0.m";//       q=2, n=10, kerDim=0, #C=100
load "./data/Cq2n20ker5.m";//       q=2, n=20, kerDim=5, #C=128
load "./data/CHamq2n15ker11.m";//   q=2, n=15, kerDim=11,#C=2048


C:= RepetitionQaryCode(GF(2), 2);// #C=2, linear
TestAllMethods(fileOutput, methods, C, 30);
print("acabat test 1");

C:= UniverseQaryCode(GF(2), 2);//   #C=4, linear
TestAllMethods(fileOutput, methods, C, 30);
print("acabat test 2");

C :=  createCode(CHadq4n4ker2);//     q=4, n=4,  kerDim=2, #C=16, is linear
TestAllMethods(fileOutput, methods, C, 6);
print("acabat test 3");

C :=  createCode(CHadq4n8ker1);//     q=4, n=8,  kerDim=1, #C=32, #cosetRep=8
TestAllMethods(fileOutput, methods, C, 5);
print("acabat test 4");

C :=  createCode(Cq4n4ker2);//     q=4, n=16, kerDim=2, #C=32, #cosetRep=2
TestAllMethods(fileOutput, methods, C, 5);
print("acabat test 5");

C :=  createCode(CHadq4n12ker1);//     q=4, n=12, kerDim=1, #C=48, #cosetRep=12
TestAllMethods(fileOutput, methods, C, 5);
print("acabat test 6");

C :=  createCode(CHadq4n16ker1c);//  q=4, n=16, kerDim=1, #C=64, #cosetRep=16
TestAllMethods(fileOutput, methods, C, 4);
print("acabat test 7");

C :=  createCode(CHadq4n16ker1b);//   q=4, n=16, kerDim=1, #C=64, #cosetRep=16
TestAllMethods(fileOutput, methods, C, 4);
print("acabat test 8");

C :=  createCode(CHadq4n16ker2a);//   q=4, n=16, kerDim=2, #C=64, #cosetRep=4
TestAllMethods(fileOutput, methods, C, 4);
print("acabat test 9");


C :=  createCode(CHadq4n16ker3);//    q=4, n=16, kerDim=3, #C=64, linear
TestAllMethods(fileOutput, methods, C, 4);
print("acabat test 10");

//------- los que tardan
L := [ VectorSpace(GF(2, 2), 2) |
    [ GF(2, 2) | a, a^2 ] where a := GF(2, 2).1,
    [ GF(2, 2) | 1, 0 ] where a := GF(2, 2).1];
C := QaryCode(L);// #C=2, kerDim=0
TestAllMethods(fileOutput, methods, C, 18);// en r=18 tarda 5 min
print("acabat test 11");

C :=  createCode(Cq4n16ker0);//       q=4, n=16, kerDim=0, #C=16, 
TestAllMethods(fileOutput, methods, C, 5);//en r=5 tarda 10 min.
print("acabat test 12");
