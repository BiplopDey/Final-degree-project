/****************************************************************/
/*                                                              */
/* Project name: Qary nonlinear codes in MAGMA	                */
/* Test file name:   TestPerformanceDirectSumQ.m                */
/*                                                            	*/
/* Comments: Performance test for the function Direct sum of    */
/*           a list of codes included in the                    */
/*           QaryCodes_ConstructionsBiplop.m file               */
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
/* Function name: TestPerformanceDirectSum                      */
/* Parameters: C, r, typeAlgMethod                              */
/* Procedure description: Construct the test performance        */
/*    of the Direct sum                                         */
/* Input parameters description:                                */
/*   - C: A q-ary code                                          */
/*   - r: Number repetitions codes of C                         */
/*   - typeAlgMethod: Method to test the performance            */
/* Output parameters description:                               */
/*   - Average time of 20 reptitions of the Direct sum          */
/*     construction                                             */
/*                                                              */
/****************************************************************/
function TestPerformanceDirectSum(C,r, typeAlgMethod) 

    tstart := Cputime();
    for iteration in [1..20] do
        if typeAlgMethod eq "Iterative" then
            directSumCode :=  DirectSumI([C^^r]);
        else
            directSumCode :=  DirectSum([C^^r]);
        end if;
    end for;
    tend := Cputime(tstart);
    tendAvg := tend/20;
  
    return tendAvg;

end function;

/****************************************************************/
/*                                                              */
/* Function name: DrawChart                                     */
/* Parameters: r, timeVersions                                  */
/* Function description:                                        */
/* Input parameters description:                                */
/*   - r: Number repetitions codes of C                         */                                   
/*   - methods: Array that has the execution                    */
/*   - typeAlgMethod: Method to test the performance            */                               
/*                time for the different methods of Direct sum  */
/*                                                              */   
/****************************************************************/
procedure DrawChart(r, methods, timeVersions)

    printf "  num. codes = %o,", r;
    for j:= 1 to #methods do
        printf " %o = %o, ", methods[j], timeVersions[j];
    end for;
    printf " \n";
  
end procedure;

/****************************************************************/
/*                                                              */
/* Function name: TestAllMethods                                */
/* Parameters: fileOutput, methods, C, r                        */
/* Function description: Given a q-ary code C compute a  the    */
/* average time spending of constricting Direct sum in each     */
/* method of the vector consisting of r times C.                */     
/* Input parameters description:                                */
/*   - fileOutput: Name of the file where the outputs of the    */ 
/*                 test execution are printed                   */
/*   - methods: Array that has the execution                    */
/*   - C: A q-ary code                                          */
/*   - r: Number repetitions codes of C                         */  
/*                                                              */
/****************************************************************/
procedure TestAllMethods(fileOutput, methods, C, r)

    SetOutputFile(fileOutput);
    printf "C: %o \n", C;
    printf "Dimension of kernel: %o , num. codewords = %o\n", Dimension(C`Kernel), #C;
    printf "---------Results-------\n";
    UnsetOutputFile();

    SetOutputFile(fileOutput);
    for i in {2..r} do
        timeVersions := [];
        for method in methods do
            TF := TestPerformanceDirectSum(C,i, method);
            Append(~timeVersions, TF);
        end for;
        DrawChart(i, methods, timeVersions);
    end for;

    printf "-----------------------\n";
    UnsetOutputFile();

    SetOutputFile(fileOutput);
    printf "Test %o acabat \n \n", fileOutput;
    UnsetOutputFile();
    
end procedure; 

