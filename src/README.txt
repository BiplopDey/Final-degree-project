             The "Qary Codes" Package
             --------------------------

"Qary Codes" is  a package  that supports  the basic  facili-
ties to  work with  q-ary codes (see [1]) in  Magma (see [2]).
With this package, this kind of codes can be created  and mani-
pulated and information about these codes can be calculated.

"Qary Codes" consists of files written in the Magma language.
Please send your bug reports to Combinatorics, Coding and Secu-
rity Group (CCSG) at  support-ccsg@deic.uab.cat  or if  it is a
Magma problem  to  magma-trouble (magma@maths.usyd.edu.au). See
the section below.

"Qary Codes" was originally written in Magma by several assis-
tant developers supervised by Mercè Villanueva and  Jaume Pujol
as  a support for a research project on q-ary codes  developed
by  the  Combinatorics, Coding and Security Group (CCSG) within
the  Department of  Information  and Communications Engineering
(dEIC) at the Autonomous University of Barcelona (UAB).

This version of the package has been developped in Magma version
2.23-4.


                    Composition of the package
                        ------------------

The   "Qary  Codes"   package  is  composed  of  four directo-
ries:

/src: The files to attach to Magma "QaryCodes_Core.m", 
      "QaryCodes_Extension.m", "QaryCodes_Distances.m".
/license: The license of the package.
/doc: The manual to use the package in pdf format.
/examples: Some examples  from  the  manual.  They can  be loaded
           in Magma as soon as the package is attached.
/test: Some test files that can be used to check the package.

            Using/Installing "Qary Codes"
            -------------------------------

To use  "Qary Codes"  temporally  (as a Magma Package)
unpack  the  archive  file in a directory.   Enter to the ./src
directory. Call Magma and then write:
   Attach("QaryCodes_Core.m");
   Attach("QaryCodes_Extension.m");
   Attach("QaryCodes_Distances.m");
in the Magma command line.

To install "Qary Codes" permanent (as a Magma Package):

1. Unpack the archive file in a directory.

2. Enter  to  the  directory  where  Magma  is installed, go to
   package directory    $PATHMAGMA/package/    and create a new
   directory.

     mkdir QaryCodes

3. Copy $PATH/src/QaryCodes_Core.m $PATH/src/QaryCodes_Extension.m
    $PATH/src/QaryCodes_Distances.m
    to the Magma directory $PATHMAGMA/package/QaryCodes/ 

   cp $PATH/src/QaryCodes_Core.m
              $PATHMAGMA/package/QaryCodes/
            
   cp $PATH/src/QaryCodes_Extension.m
              $PATHMAGMA/package/QaryCodes/
            
   cp $PATH/src/QaryCodes_Distances.m
              $PATHMAGMA/package/QaryCodes/            

4. Edit the file      $PATHMAGMA/package/spec     and write the
   following lines at the end:

     QaryCodes
     {
        QaryCodes_Core.m
        QaryCodes_Extension.m
        QaryCodes_Distances.m
     }


                             Bug reports
                             -----------

When  sending a  bug  report to support-ccsg@deic.uab.cat or to
magma@maths.usyd.edu.au,    remember we will need to be able to
reproduce the problem; so please include:

 * The  version  of  Magma  you  are  using; either look at the
   header when you start up Magma.
 * The  operating  system you are using e.g. Linux, SunOS 5.8 =
   Solaris 2.8, IRIX 6.5, Windows, ...
 * A script that demonstrates the bug, along with a description
   of why it's a bug (e.g.  by  adding comments to  the  script
   - recall  comments  in Magma  begin  with  a  //  or between
   /*  */).


                             Bibliography
                             ------------

[1] M. Villanueva, F. Zeng, and J. Pujol,
    "Efficient representation of binary nonlinear codes: constructions
    and minimum distance computation", Designs, Codes and Cryptography, 
    vol. 76(1), pp. 3-21, 2015.

[2] J.J. Cannon and W. Bosma (Eds.) "Handbook of Magma Functions",
   Edition 2.13, 4350 pages, 2006.


July 31, 2017
