/****************************************************************/
/*                                                              */
/* Project name: Qary nonlinear codes in MAGMA	                */
/* Test file name:   TestPerformanceFix.m                       */
/*                                                            	*/
/* Comments: Performance test for the function Fix              */
/*           included in the QaryCodes_ActionGroup.m file       */
/*                                                              */
/* Authors: A. Figuerola and M. Villanueva                      */
/*                                                            	*/
/* Revision version and last date: 1.1,  2020/08/8            	*/
/*                                                              */
/*                                                              */
/****************************************************************/

Attach("QaryCodes_Core.m");
Attach("QaryCodes_Extension.m");
Attach("QaryCodes_Distances.m");
Attach("QaryCodes_Constructions.m");

/****************************************************************/
/*                                                              */
/* Function name: TestPerformanceC_G                            */
/* Parameters: C, G, typeAlgMethod                              */
/* Procedure description: Construct the test performance        */
/*   of the function C^G considering artifical codes with       */
/*   diferent dimension of the kernel from the same code.       */                                     
/* Input parameters description:                                */
/*   - C: A q-ary code                                          */
/*   - G: A group of permutations                               */
/*   - typeAlgMethod: Method to test the performance            */
/* Output parameters description:                               */
/*   - A sequence with the execution times of C^G by using the  */
/*     given method and considering pseudo-codes constructed    */
/*     from C with different dimensions of the kernel           */
/*                                                              */
/****************************************************************/
function TestPerformancePlotkinSum(C, D, typeAlgMethod) 

    plotkinSumCD := PlotkinSum(C, D); 

    kernel := C`Kernel;
    cosetRep := C`CosetRepresentatives;

    maxDimPartialK := Dimension(kernel);
    minDimPartialK := 0; // ponerlo igual al maxDimPartialK si se usa PlotkinSum, 
                         // ya que si a un codigo le da el partial kernel como final,
                         // El brutforce usa QaryCode que este calcula el kernel total y por tanto es mas costoso.

    timePerformance := [];

    for dimPartialK := maxDimPartialK to minDimPartialK by -1 do
        //partial kernel and coset representatives
        partialK := Subcode(kernel, dimPartialK);
        compK := CodeComplement(kernel, partialK);
        partialRep := [k + c : k in Set(compK), c in cosetRep ];
        newCode := QaryCode(partialK, partialRep : IsFinalKernel := true);

        tstart := Cputime();
        for iteration in [1..20] do
            if typeAlgMethod eq "BruteForce" then
                plotkinSumCode :=  PlotkinSumBF(newCode, D);
            else
                plotkinSumCode :=  PlotkinSum(newCode, D);
            end if;
        end for;
        tend := Cputime(tstart);
        tendAvg := tend/20;
        Append(~timePerformance, tendAvg);

        i := maxDimPartialK - dimPartialK + 1;
        //printf "--------------k = %o---------------\n ", dimPartialK;
        //printf "#PAutC: %o \n", #PAutC;
        //printf "#PAutCNew: %o \n", #PAutCNew;

        assert Set(plotkinSumCD) eq Set(plotkinSumCode);// para comprovar que sean iguales

    end for;

    return timePerformance;

end function;

/****************************************************************/
/*                                                              */
/* Function name: DrawChart                                     */
/* Parameters: maxDimPartialK, timeVersions                     */
/* Function description:                                        */
/* Input parameters description:                                */
/*   - maxDimPartialK: Maximum dimension of the partial kernel  */                                   
/*   - timeVersions: Bidimensional array that has the execution */
/*                   time for the different methods of C^G      */
/*        timeVersions[1] time sequence for the brute force     */
/*        timeVersions[2] time sequecen for the CosetRepV1      */
/*        timeVersions[3] time sequence for the CosetRepV2      */      
/*                                                              */
/****************************************************************/
procedure DrawChart(maxDimPartialK, methods, timeVersions)

    for i := 0 to maxDimPartialK do 
        printf "k = %o, ", maxDimPartialK-i;
        for j:= 1 to #methods do
            printf " %o = %o, ", methods[j], timeVersions[j][i+1];
        end for;
        printf " \n";
    end for;

end procedure;

/****************************************************************/
/*                                                              */
/* Function name: TestAllMethods                                */
/* Parameters: fileOutput, C, seqG, PAutC                       */
/* Function description: Given a q-ary code C and a sequence    */
/*   of groups [G1, G2, ...], return a sequence with the        */
/*   orbit of the code for each group [C^G1, C^G2, ...].        */
/* Input parameters description:                                */
/*   - fileOutput: Name of the file where the outputs of the    */ 
/*                 test execution are printed                   */
/*   - C: A q-ary code                                          */
/*   - seqG: A sequence of permutation groups                   */
/*   - PAutC: The permutation automorphism group of C           */
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

    maxDimPartialK := Dimension(C`Kernel);

    printf "---------Results-------\n";
    DrawChart(maxDimPartialK,  methods, timeVersions);
    printf "-----------------------\n";
    UnsetOutputFile();

    SetOutputFile(fileOutput);
    printf "Test %o acabat \n", fileOutput;
    UnsetOutputFile();
    
end procedure; 

