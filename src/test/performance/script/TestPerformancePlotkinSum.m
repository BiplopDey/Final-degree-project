/****************************************************************/
/*                                                              */
/* Project name: Qary nonlinear codes in MAGMA	                */
/* Test file name:   TestPerformancePlotkinSum.m                */
/*                                                            	*/
/* Comments: Performance test for the function Plotkin sum      */
/*           included in the QaryCodes_Constructions.m    */
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
Attach("QaryCodes_Constructions.m");

/****************************************************************/
/*                                                              */
/* Function name: TestPerformancePlotkinSum                     */
/* Parameters: C, D, typeAlgMethod                              */
/* Procedure description: Construct the test performance        */
/*    of the Plotkin sum                                        */
/* Input parameters description:                                */
/*   - C: A q-ary code                                          */
/*   - D: A q-ary code                                          */
/*   - typeAlgMethod: Method to test the performance            */
/* Output parameters description:                               */
/*   - Average time of 20 reptitions of the Plotkin sum         */
/*     construction                                             */
/*                                                              */
/****************************************************************/
function TestPerformancePlotkinSum(C, D, typeAlgMethod) 

    tstart := Cputime();
    for iteration in [1..20] do
        if typeAlgMethod eq "BruteForce" then
            plotkinSumCode :=  PlotkinSumBF(C, D);
        else
            plotkinSumCode :=  PlotkinSum(C, D);
        end if;
    end for;
    tend := Cputime(tstart);
    tendAvg := tend/20;
  
    return tendAvg;

end function;

/****************************************************************/
/*                                                              */
/* Function name: DrawChart                                     */
/* Parameters: C, D, timeVersions                               */
/* Function description:                                        */
/* Input parameters description:                                */
/*   - C: A q-ary code                                          */
/*   - D: A q-ary code                                          */    
/*   - methods: Array that has the execution                    */                              
/*   - timeVersions: Array that has the execution time for the  */
/*                  different methods of Plotkin sum            */
/*                                                              */   
/****************************************************************/
procedure DrawChart(C, D, methods, timeVersions)
    dimKC:=Dimension(C`Kernel);
    dimKD:=Dimension(D`Kernel);
    NumberCodewords:=(#C)*(#D);
    sumDim:=dimKC+dimKD;

    printf "k_1 = %o, k_2 = %o, k = %o, num. codewords = %o,", dimKC, dimKD, sumDim, NumberCodewords;
    for j:= 1 to #methods do
        printf " %o = %o, ", methods[j], timeVersions[j];
    end for;
    printf " \n";
  
end procedure;

/****************************************************************/
/*                                                              */
/* Function name: TestAllMethods                                */
/* Parameters: fileOutput, methods, C, D                        */
/* Function description: Given C and D q-ary codes compute the  */
/* average time spending of constricting Plotking sum in each   */
/* method.                                                      */     
/* Input parameters description:                                */
/*   - fileOutput: Name of the file where the outputs of the    */ 
/*                 test execution are printed                   */
/*   - methods: Array that has the execution                    */
/*   - C: A q-ary code                                          */
/*   - D: A q-ary code                                          */
/*                                                              */
/****************************************************************/
procedure TestAllMethods(fileOutput, methods, C, D)

    SetOutputFile(fileOutput);
    printf "C: %o \n", C;
    printf "D: %o \n", D;
    UnsetOutputFile();

    timeVersions := [];    
    SetOutputFile(fileOutput);
    for method in methods do
        TF := TestPerformancePlotkinSum(C, D, method);
        Append(~timeVersions, TF);
    end for;

    printf "---------Results-------\n";
    DrawChart(C,D, methods, timeVersions);
    printf "-----------------------\n";
    UnsetOutputFile();

    SetOutputFile(fileOutput);
    printf "Test %o acabat \n \n", fileOutput;
    UnsetOutputFile();
    
end procedure; 

