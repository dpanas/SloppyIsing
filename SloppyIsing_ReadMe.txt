  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++


	This file describes the workflow of the SloppyIsing package. 
	Following the steps below will allow the user to reproduce the main results of the paper:

	Sloppiness in Spontaneously Active Neuronal Networks
	D. Panas, H. Amin, A. Maccione, O. Muthmann, M. van Rossum, L. Berdondini, M.H. Hennig
	The Journal of Neuroscience 2015


	All the scripts and functions provided here were written by me, DAP (D. Panas).
	Any mistakes, bugs, etc, are solely my responsibility. If you have any comments or criticism,
	please e-mail me: 
					dagger _at_ autograf.pl

	The code is free to use, share and modify. 
	Should it inspire your research or grow into something bigger and better, I will be glad to know.
	If you feel my code has helped you in your research in any appreciable degree, feel free to cite our work :-)


	The scripts are all written in Matlab and some of the scripts use functions from Statistical Toolbox.


  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++


	BASIC INFORMATION NEEDED TO REPRODUCE RESULTS / SLOPPY ISING DEMO


These steps generate the maximum entropy model fits and Fisher Information Matrix decompositions for 160 groups of randomly sampled neurons:

1. script_SampleRandomGroups
- run the script, it is set up to sample 8 neuron groups from the population of common channels between 6 recordings from a hippocampal culture
(8 neuron groups, instead of 10, to save time and to show that our results hold also for different group size)

2. script_Ising_Main
- run the script, it is set up to use the result of previous script
- inspect whether all groups converged (information will be displayed on the screen and plots will be generated)
- few groups might not converge, but likely it will be minor discrepancies from the data that can be ignored
(this can be verified in the plots)
! this script takes a long time - it took about 15 hours on a desktop machine (4GB RAM, Intel Core i3-2100 3.1GHz, 2 cores, 4 threads) 

2a. script_Ising_Additional 
- if there are unconverged groups this script can be run to correct the convergence
- should there still remain unconverged groups, learn_rate_ising2 needs to be changed (e.g. 0.99), and max_iter2 can be changed as well

3. script_FIM_Decomposition
- run the script, it is set up to use results of Ising analysis described above
- it automatically generates results that can then be plotted


To view / compare results, the following scripts are provided, reproducing some of the figures from the paper. All are currently set up
in demo mode to use the results of above analysis:

4. script_Compare_RatesCorrelations
- run this script to inspect the distributions of firing rates and correlations, and compute some statistical tests

5. script_Compare_Results
- run this script to inspect sloppiness and sparsity, and see how groups of neurons differ between recordings in various measures

6. script_Compare_Projections
- run this script to see how projections of parameters onto eigenparameters of FIM change over time

7. script_Compare_ChangesVsSensitivity
- run this script to see if there is any relationship between changes between recordings and presumed sensitivity of corresponding parameters


  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++


	EXTENDED INFORMATION FOR INTERESTED USERS

The purpose of this package is to be able to analyse appropriately formatted files containing spikes recorded from neural preparations.
This analysis is described in detail in our paper. The general idea behind this work was to identify units active throughout an experiment, 
then subsample the population extensively, so that many random groups of neurons can be fitted with a Pairwise Maximum Entropy model. 
Once the model is satisfactorily fit to the data, it can then be further analysed with Fisher Information Matrix Decomposition, which provides
insight into the parameter sensitivity of each model.

As a sample to work on, 6 spike files are provided in the ./results_spikes/ folder, the same data as in the paper (Culture 1).
The functions necessary to run all the scripts are in the ./functions/ folder. The remaining folders are empty, prepared for the results
of the DEMO analysis.

The DEMO analysis is set up for almost automatic generation of results and plots, but for anyone interested in further digging into the data, 
or analysis of their own data (how spike files are formatted is described at the end), all scripts come with a short description that
can be viewed in the source file, or by calling 'help script_OfInterest' in Matlab command line. This should be sufficient for you to be able
to change any adjustable variables in the scripts (which are all grouped together) and generate your own new results.

General tip about convergence:
The more asychronous the data, the smaller the difference between Independent and Ising - for some groups it might happen that Ising does not appear
to offer an improvement over the Independent (which will be apparent as Multiinformatio ratio out of range of [0 1]). Also, although the algorithm
is theoretically bound to converge, in practice it might sometimes not, since this is discrete data, and this can sometimes result in very big differences
in spike counts. All these cases will be caught by Ising_Main and flagged as uncoverged groups. From experience, changing the learning rate to rather low
values, such as 0.2 or 0.3, and increasing the number of maximum iterations helps. 

Extra scripts:

1. script_Ising_Main_ResampleRefit
- run this script to create artificial spike trains from the modelled distributions, and re-fit them again with Ising;
this can be used later on as a check against overfitting; currently set up for 20 re-samples
! this script takes a long time - it took about 44 hours on a desktop machine (4GB RAM, Intel Core i3-2100 3.1GHz, 2 cores, 4 threads) 

2. script_FIM_Decomposition_ResampleRefit
- run this script after the previous one, to get FIM decompositions for re-samples;
the output of this can then be included in the 'script_Compare_ChangesVsSensitivity' to plot alongside the main results; 

3. script_Get_Spikes
- this is a script meant as an illustration, to outline how spike files are formatted;
- if you want to use it for your own spikes, I would advise copying the GetSpikesText function and re-writing bits of it to suit the particular spike
data you have, e.g. a differnt array size, or input format not in .txt but .mat or .hdf5, etc. 

  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++  ++

									THE END

