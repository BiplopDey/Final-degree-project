/****************************************************************/
/*                                                              */
/* Project name: Qary nonlinear codes in MAGMA	                */
/* Test file name:   TestPerformancePadCode.m                   */
/*                                                            	*/
/* Comments: Performance test for the function Pad code         */
/*           included in the QaryCodes_ConstructionsBiplop.m    */
/*           file                                               */
/*                                                              */
/* Authors: Biplop Dey and M. Villanueva                        */
/*                                                             	*/
/* Revision version and last date: 1.0,  28/05/2021            	*/
/*                                                              */
/*                                                              */
/****************************************************************/

Attach("QaryCodes_Core.m");
Attach("QaryCodes_Extension.m");
Attach("QaryCodes_Distances.m");
Attach("QaryCodes_ConstructionsBiplop.m");

/****************************************************************/
/*                                                              */
/* Function name: TestPerformancePadCode                        */
/* Parameters: C, D, typeAlgMethod                              */
/* Procedure description: Construct the test performance        */
/*    of the Pad code of C and n, where n is the lengt of C     */
/* Input parameters description:                                */
/*   - C: A q-ary code                                          */
/*   - typeAlgMethod: Method to test the performance            */
/* Output parameters description:                               */
/*   - Average time of 20 reptitions of the Pad code            */
/*     construction                                             */
/*                                                              */
/****************************************************************/
function TestPerformancePadCode(C, typeAlgMethod) 
    n:=C`Length;

    tstart := Cputime();
    for iteration in [1..20] do
        if typeAlgMethod eq "BruteForce" then
            newPadCode :=  PadCodeBF(C, n);
        else
            newPadCode :=  PadCode(C, n);
        end if;
    end for;
    tend := Cputime(tstart);
    tendAvg := tend/20;
      
    return tendAvg;

end function;

/****************************************************************/
/*                                                              */
/* Function name: DrawChart                                     */
/* Parameters: dimK, methods, timeVersions                      */
/* Function description:                                        */
/* Input parameters description:                                */
/*   - dimK: Dimension of the kernel of C                       */  
/*   - methods: Array that has the execution                    */                              
/*   - timeVersions: Array that has the execution time for the  */
/*                  different methods of Pad code               */
/*                                                              */   
/****************************************************************/
procedure DrawChart(dimK, methods, timeVersions)

        printf "k = %o, ", dimK;
        for j:= 1 to #methods do
            printf " %o = %o, ", methods[j], timeVersions[j];
        end for;
        printf " \n";
   
end procedure;

/****************************************************************/
/*                                                              */
/* Function name: TestAllMethods                                */
/* Parameters: fileOutput, methods, C                           */
/* Function description: Given a q-ary code C compute the       */
/* average time spending of constricting Pad code in each       */
/* method.                                                      */     
/* Input parameters description:                                */
/*   - fileOutput: Name of the file where the outputs of the    */ 
/*                 test execution are printed                   */
/*   - methods: Array that has the execution                    */
/*   - C: A q-ary code                                          */
/*                                                              */
/****************************************************************/
procedure TestAllMethods(fileOutput, methods, C)

    SetOutputFile(fileOutput);
    printf "C: %o \n", C;
    UnsetOutputFile();

    timeVersions := [];    
    SetOutputFile(fileOutput);
    for method in methods do
        TF := TestPerformancePadCode(C, method);
        Append(~timeVersions, TF);
    end for;

    dimK := Dimension(C`Kernel);

    printf "---------Results-------\n";
    DrawChart(dimK,  methods, timeVersions);
    printf "-----------------------\n";
    UnsetOutputFile();

    SetOutputFile(fileOutput);
    printf "Test %o acabat \n", fileOutput;
    UnsetOutputFile();
    
end procedure; 

