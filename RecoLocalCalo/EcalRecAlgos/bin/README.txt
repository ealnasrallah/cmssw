// 
*** Intstructions for CMS Ecal Rechit Validation Codes ***
*** Eissa Alnasrallah, August 23, 2019 ***
//
The instructions are for files in the following directory: RecoLocalCalo/EcalRecAlgos/bin/
/ Contents:
-BuildFile.xml
-makeEcalRechitValidationPlots.cpp
-makeEcalMultifitResultsGpuValidationPlots.cpp
-makeEcalRechitEtaValidationPlots.cpp
-rechitValTest.cpp
-uncalibRechitValTest.cpp
//
//
-- BuildFile.xml:
Any file that needs to be compiled should be included in the buildfile with the following format:

<bin name= "NAME_OF_EXECUTABLE" file="FILE_NAME.cpp">
    <use name="root"/>
    <use name="rootgraphics"/>
    <use name="CUDADataFormats/EcalRecHitSoA"/>
    <use name="DataFormats/Common"/>
    <use name="DataFormats/EcalRecHit"/>
</bin>

where "FILE_NAME.cpp" is the file to be compiled; "NAME OF EXECUTABLE" is the command to run in the shell.
//
//
-- makeEcalRechitValidationPlots.cpp (*********************************** The main validation code ***********************************)
This is the main validation code for Ecal Rechits. It currently tests for: event size(No. of Rechits), did, energy, chi2, flags, and extras. It produces six plots, 3 for barrel (EB) and 3 for end-cap (EE):
a) The values in GPU and CPU
b) 2D histogram of GPU vs CPU. so that CPU value is the x-axies and GPU value is the y-axis.
c) 1D histogram of (GPU/CPU) ratio.
+) Plots for delta vs CPU are included but currently commented since it is not needed currently. It is useful for more accurate comparisons.
The code retrieves variables of each rechit (did, energy, flag,...) on GPU, then allocates the CPU rechit with the same "did" and compares the variables. and the plots the results for EB and EE individually.
//
** Commented parts of the code:
1- There is an assurance part in the code to ensure that GPU and CPU Events are of the same Rechit sizes. This is commented since thus far not all GPU rechits are resolved.
2- For Both EB and EE, in case a GPU rechit had a DetId that code fails to find an equivalent for on CPU a message will appear. This is commented for now.
//
** Notes
To call for "flags" and "extras" changes were done to DataFormats/EcalRecHit/interface/EcalRecHit.h to include them. Please refer to the file on my branch if needed.
**
To run the File:
1- cmsenv
2- Compile: scramv1 b -j 10
3- run the code: makeEcalRechitValidationPlots <input> <output>
      Example: makeEcalRechitValidationPlots test.root result.root
4- The output should appear in the directory along with .root files for each variable of the rechits. 
//
//
-- makeEcalMultifitResultsGpuValidationPlots.cpp
This is the validation code for the Uncalibrated rechits (for amplitudes) which was the work of Viktor and Andrea. It was the basis for the validation code for the rechits. I modified this version to present the results in a similar fashion to the Rechits
** To Run:
1- cmsenv
2- Compile: scramv1 b -j 10
3- run the code: makeEcalMultifitResultsGpuValidationPlots <input> <output>
      Example: makeEcalMultifitResultsGpuValidationPlots test.root result.root
4- The output should appear in the directory along with .root files for each variable of the uncalibrated rechits.
//
//
-- makeEcalRechitEtaValidationPlots.cpp
This code runs only the Energies and produces two .root files:
1-The energies of EB and EE for the rechits (Just like the "makeEcalRechitValidationPlots" plot)
2- plots of the GPU/CPU energy ration for three different regions of eta.
The advantage is that it is a much simpler code to go through and only has the energy parts. eta is plotted for endcaps (EE).
** To Run:
1- cmsenv
2- Compile: scramv1 b -j 10
3- run the code: makeEcalRechitEtaValidationPlots <input> <output>
      Example: makeEcalRechitEtaValidationPlots test.root result.root
4- The output should appear in the directory along with .root files for the energies of the rechits and energies of different eta regions.
//
//
-- rechitValTest.cpp
This is a compact validation code that produces any variable (currently set to did) in a text file. It is helpful for debugging in case plots can't be produced in the main validation code.
** The code only produces results for endcap (EE) GPU part.
** can be adjusted to produce EB results by swtiching all 'ee' to 'eb
It makes sure we have values from GPU for the variables and allows for checking the values.
** To Run:
1- cmsenv
2- Compile: scramv1 b -j 10
3- run the code: valRechitTest <input> <output>
      Example: valRechitTest test.root result.root
4- The output should appear in the directory along with 'test_ee.txt'. 
//
//
-- uncalibRechitValTest.cpp
This is a compact validation code that produces any variable (currently set to amplitude) of the uncalibrated rechits in a text file. It is helpful for debugging and to compare to the rechits values.
** The code only produces results for endcap (EE) GPU part.
** can be adjusted to produce EB results by swtiching all 'ee' to 'eb
It allows for comparison between calibrated and uncalibrated results if needed.
** To Run:
1- cmsenv
2- Compile: scramv1 b -j 10
3- run the code: valUncalibRechitTest <input> <output>
      Example: valUncalibRechitTest test.root result.root
4- The output should appear in the directory along with 'test_ee_unc.txt'. 
******************
 test files are produced by running: 
 cmsRun <input.py> <output.root>