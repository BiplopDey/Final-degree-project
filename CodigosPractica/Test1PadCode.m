/****************************************************************/
/* Project name: Qary nonlinear codes in MAGMA    Test  4sd     */
/****************************************************************/

load "TestPerformancePadCode.m";
load "TestPerformanceGeneral.m";
methods := ["BruteForce", "KernelRep"];
fileOutput := createFile("TestPadCodeKCvsBF");

load "./data/Cq4n16ker0.m";//       q=4, n=16, kerDim=0, #C=32
load "./data/Cq4n4ker2.m";//        q=4, n=16, kerDim=2, #D=32
load "./data/CHadq4n16ker3t.m";//   q=4, n=16, kerDim=3, #D=64

load "./data/Cq2n10ker0.m";//       q=2, n=10, kerDim=0, #D=100
load "./data/Cq2n20ker5.m";//       q=2, n=20, kerDim=5, #C=128
load "./data/CHamq2n15ker11.m";//   q=2, n=15, kerDim=11,#C=2048
load "./data/COsterq2n18ker3.m";//  q=2,n=18,kerDim=3, #C=5632
load "./data/COsterq2n25ker7.m";//  q=2,n=25,kerDim=7,#C=17920
load "./data/Cq3n13ker9c.m";//      q=3,n=13,kerDim=9, #C=59049
load "./data/Cq3n13ker10.m";//      q=3,n=13,kerDim=10,#C=59049
load "./data/Cq3n13ker8b.m";//      q=3,n=13,kerDim=8,#C=59049
load "./data/CKerq2n256ker9.m";//   q=2,n=25,kerDim=9,#C=65536
load "./data/COsterq2n24ker12.m";//   q=2,n=24,kerDim=12,#C=327680


C := createCode(Cq4n16ker0); //q=4, n=16, kerDim=0, #C=32
TestAllMethods(fileOutput, methods, C);
print("acabat  test 1");

C := createCode(Cq4n4ker2);//q=4, n=16, kerDim=2, #D=32
TestAllMethods(fileOutput, methods, C);
print("acabat  test 2");

C := createCode(CHadq4n16ker3t);//q=4, n=16, kerDim=3, #D=64
TestAllMethods(fileOutput, methods, C);
print("acabat  test 3");

C:=createCode(Cq2n10ker0);//q=2, n=10, kerDim=0, #D=100,
TestAllMethods(fileOutput, methods, C);
print("acabat  test 4");

C:=createCode(CHamq2n15ker11);//q=2, n=15, kerDim=11, #C=2048
TestAllMethods(fileOutput, methods, C);
print("acabat  test 5");

C:=createCode(Cq2n20ker5);// q=2,n=20,kerDim=5,#C=128
TestAllMethods(fileOutput, methods, C);
print("acabat  test 6");

C:=createCode(COsterq2n18ker3);// q=2,n=18,kerDim=3, #C=5632
TestAllMethods(fileOutput, methods, C);
print("acabat  test 7");

C:=createCode(COsterq2n25ker7);// q=2,n=25,kerDim=7,#C=17920
TestAllMethods(fileOutput, methods, C);//tarda 7 min
print("acabat  test 8");

C:=createCode(Cq3n13ker9c);// q=3,n=13,kerDim=9, #C=59049
TestAllMethods(fileOutput, methods, C);
print("acabat  test 9");

C:=createCode(Cq3n13ker10);// q=3,n=13,kerDim=10,#C=59049
TestAllMethods(fileOutput, methods, C);
print("acabat  test 10");

C:=createCode(Cq3n13ker8b);// q=3,n=13,kerDim=8,#C=59049
TestAllMethods(fileOutput, methods, C);
print("acabat  test 11");

C:=createCode(CKerq2n256ker9);//   q=2,n=25,kerDim=9,#C=65536
TestAllMethods(fileOutput, methods, C);

C:=createCode(COsterq2n24ker12);//   q=2,n=24,kerDim=12,#C=327680
TestAllMethods(fileOutput, methods, C);
print("acabat  test 13");
