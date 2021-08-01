/****************************************************************/
/* Project name: Qary nonlinear codes in MAGMA    Test  4sd     */
/****************************************************************/

load "TestPerformanceDirectSum.m";
load "TestPerformanceGeneral.m";
methods := ["BruteForce", "KernelRep"];
fileOutput := createFile("TestDirectSumKCvsBF");

load "./data/CHadq4n8ker1.m";//     q=4, n=8, kerDim=1, #D=32
load "./data/Cq4n4ker2.m";//        q=4, n=16, kerDim=2, #D=32
load "./data/CHadq4n16ker1c.m";//   q=4, n=16, kerDim=1, #C=64
load "./data/CHadq4n16ker1b.m";//   q=4, n=16, kerDim=1, #C=64
load "./data/Cq4n16ker0.m";//       q=4, n=16, kerDim=0, #C=16
load "./data/CHadq4n16ker2a.m";//   q=4, n=16, kerDim=2, #C=64
load "./data/CHadq4n16ker3t.m";//   q=4, n=16, kerDim=3, #D=64
load "./data/CHadq4n16ker3.m";//    q=4, n=16, kerDim=3, #C=64

load "./data/Cq2n20ker5.m";//       q=2, n=20, kerDim=5, #C=128
load "./data/Cq2n10ker0.m";//       q=2, n=10, kerDim=0, #D=100
load "./data/CHamq2n15ker11.m";//   q=2, n=15, kerDim=11,#C=2048

C :=  createCode(CHadq4n8ker1);//   q=4, n=8, kerDim=1, #D=32, #cosets=8
TestAllMethods(fileOutput, methods, C, C);

zeroCode := ZeroCode(GF(2), 8);
C := QaryCode(zeroCode);
TestAllMethods(fileOutput, methods, C, C);

universeCode := UniverseCode(GF(3), 2);
C := QaryCode(universeCode);
TestAllMethods(fileOutput, methods, C, C);

C := createCode(Cq4n4ker2);//q=4, n=4, kerDim=2, #D=32, #cosets=2
TestAllMethods(fileOutput, methods, C, C);

D := createCode(CHadq4n16ker1c);//q=4, n=16, kerDim=1, #C=64, #cosets=16
TestAllMethods(fileOutput, methods, D, D);

M:=KMatrixSpace(GF(3),4,7);
G:=M![0,0,0,2,0,1,1,1,2,0,0,0,0,0,1,1,2,0,0,1,1,0,2,1,0,1,0,0];
C:=QaryCode(LinearCode(G));//q=3, n=7, kerDim=4, #C=81
M:=KMatrixSpace(GF(3),3,7);
G:=M![1,1,2,0,3,1,1,0,0,1,0,1,0,2,1,1,1,0,2,0,0];
D:=QaryCode(LinearCode(G));//q=3, n=7, kerDim=3, #D=27
TestAllMethods(fileOutput, methods, C, C);
TestAllMethods(fileOutput, methods, D, D);
TestAllMethods(fileOutput, methods, C,D);

C := createCode(CHadq4n16ker1b);//q=4, n=16, kerDim=1, #C=64, #cosets=16
D := createCode(Cq4n16ker0); //q=4, n=16, kerDim=0, #C=16, #cosets=16
TestAllMethods(fileOutput, methods, C, C);
TestAllMethods(fileOutput, methods, D, D);
TestAllMethods(fileOutput, methods, C,D);

C := createCode(CHadq4n16ker2a);//q=4, n=16, kerDim=2, #C=64, #cosets=4
D:= createCode(CHadq4n16ker3t);//q=4, n=16, kerDim=3, #D=64, is linear
TestAllMethods(fileOutput, methods, C, C);
TestAllMethods(fileOutput, methods, D, D);
TestAllMethods(fileOutput, methods, C,D);
print("empiezan los lentos");

C:=createCode(Cq2n10ker0);//q=2, n=10, kerDim=0, #D=100,
M:=KMatrixSpace(GF(2),4,7);
G:=M![1,0,0,0,0,1,1,1,0,1,0,0,1,0,1,0,1,0,0,1,1,0,1,1,0,1,0,0];
D:=QaryCode(LinearCode(G));//q=2, n=7, #D=16
TestAllMethods(fileOutput, methods, C, C);// tardará 10 min
TestAllMethods(fileOutput,methods, C,D);

C:=createCode(CHamq2n15ker11);//q=2, n=15, kerDim=11, #C=2048
TestAllMethods(fileOutput, methods, C, C);

C:=createCode(Cq2n20ker5);// q=2,n=20,kerDim=5,#C=128, #coset = 4
TestAllMethods(fileOutput, methods, C, C);

load "./data/Cq2n20ker0.m";// q=2,n=20,kerDim=0,MinDistance= 4 ,notLinear, #C=127, #costeRep=127, 0 notin C.
C:=createCode(Cq2n20ker0);// q=2,n=20,kerDim=0,MinDistance= 4 ,notLinear, #C=127, #costeRep=127, 0 notin C.
TestAllMethods(fileOutput, ["KernelRep"], C, C);

load "./data/CHamq2n16ker4t.m";
C:=createCode(CHamq2n16ker4t);// q=2,n=16,kerDim=4,MinDistance=4,notLinear, #C=2048, #costeRep=128, 0 notin C.
TestAllMethods(fileOutput, ["KernelRep"], C, C);

load "./data/COsterq2n24ker0.m";
C:=createCode(COsterq2n24ker0); // q=2,n=24,kerDim=0,MinDistance= 10 ,notLinear, #C=136, #costeRep=136, 0 in C.
TestAllMethods(fileOutput, ["KernelRep"], C, C);

load "./data/Cq3n13ker9c.m";
C:=createCode(Cq3n13ker9c);// q=3,n=13,kerDim=9,MinDistance= 3 , notLinear, #C=59049, 0 in C, #cosets=3
TestAllMethods(fileOutput, ["KernelRep"], C, C);

load "./data/COsterq2n18ker3.m";//  q=2,n=18,kerDim=3, #C=5632
C:=createCode(COsterq2n18ker3);
TestAllMethods(fileOutput, ["KernelRep"], C, C);

load "./data/COsterq2n25ker7.m";//  q=2,n=25,kerDim=7,#C=17920
C:=createCode(COsterq2n25ker7);
TestAllMethods(fileOutput, ["KernelRep"], C, C);

load "./data/Cq3n13ker10.m";//      q=3,n=13,kerDim=10,#C=59049
load "./data/Cq3n13ker8b.m";//      q=3,n=13,kerDim=8,#C=59049
load "./data/CKerq2n256ker9.m";//   q=2,n=25,kerDim=9,#C=65536
load "./data/COsterq2n24ker12.m";//   q=2,n=24,kerDim=12,#C=327680


// print "Codigos no usados";
// load "./data/Cq3n13ker10.m";// q=3,n=13,kerDim=10,MinDistance= 3 , Linear, #C=59049, 0 in C.
// load "./data/Cq3n13ker10_seqgen.m";
// L := Cq3n13ker10_seq;
// kernel := LinearCode(Matrix(Cq3n13ker10`Kernel));
// C := QaryCode(kernel,Cq3n13ker10`CosetRepresentatives : IsFinalKernel := true);// lower=1,upper=4


// load "./data/Cq3n13ker9c.m";// q=3,n=13,kerDim=9,MinDistance= 3 , notLinear, #C=59049, 0 in C.
// load "./data/Cq3n13ker9c_seqgen.m";
// L := Cq3n13ker9c_seq;
// kernel := LinearCode(Matrix(Cq3n13ker9c`Kernel));
// D := QaryCode(kernel,Cq3n13ker9c`CosetRepresentatives : IsFinalKernel := true);//lower=1,upper=4


// load "./data/CHamq2n16ker4t.m";// q=2,n=16,kerDim=4,MinDistance=4,notLinear, #C=2048, #costeRep=128, 0 notin C.
// load "./data/CHamq2n16ker4t_seqgen.m";
// L := CHamq2n16ker4t_seq;
// kernel := LinearCode(Matrix(CHamq2n16ker4t`Kernel));
// print "        d2 not assigned, 1<= d2 <= 6";
// D := QaryCode(kernel,CHamq2n16ker4t`CosetRepresentatives : IsFinalKernel := true);

// load "./data/COsterq2n24ker0.m"; // q=2,n=24,kerDim=0,MinDistance= 10 ,notLinear, #C=136, #costeRep=136, 0 in C.
// load "./data/COsterq2n24ker0_seqgen.m";
// L := COsterq2n24ker0_seq;
// kernel := LinearCode(Matrix(COsterq2n24ker0`Kernel));
// C := QaryCode(kernel,COsterq2n24ker0`CosetRepresentatives : IsFinalKernel := true);//aqui no se conserva min distance
// d1:=MinimumDistance(C);

// load "./data/Cq2n20ker0.m";// q=2,n=20,kerDim=0,MinDistance= 4 ,notLinear, #C=127, #costeRep=127, 0 notin C.
// load "./data/Cq2n20ker0_seqgen.m";
// L := Cq2n20ker0_seq;
// kernel := LinearCode(Matrix(Cq2n20ker0`Kernel));
// D := QaryCode(kernel,Cq2n20ker0`CosetRepresentatives : IsFinalKernel := true);// lower =1, upper=14
