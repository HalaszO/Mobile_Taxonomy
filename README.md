# Mobile_Taxonomy
MATLAB-based software kit for electrophysiological data analysis, co-developed by MTA KOKI (Institute of Experimental Medicine Budapest) and University of Szeged.

Contains multiple MATLAB function files developed outside of the project. All credits of those belong to the original creators.

Required input: text files containing the measured membrane voltage values for each sweep (current injection) in a columnal format. (1st column: 1st sweep, 2nd column: 2nd sweep...) The file must not contain any metadata, only the voltage values. To edit the time  and currents of the experiment protocol, you must do that in the script called "convertIVfromascii".

To use it, first run the script called "convertIVfromascii". This is a script for converting raw ASCII-format txt files into .mat metadata files. It also filters out any corrupted text files and selects the usable in-vitro measurement data files.
The second one ("NGtax_MAIN") is the main analysis script. It extracts certain features (see the attached document called "Field descriptions") solely based on the measurement data and saves them into .mat files. It also creates a summary for each cell and plots some of the features (e.g. number of action potentials relative to the injected current). 
The datasum2table is an conversion script, designed to write the contents of the 'datasum.mat' summary into a .dat-format file.
