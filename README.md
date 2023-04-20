# NetworkPlastic

To recover plots go to folder scriptFigs. 

Fig. 5 Bicrystal rotatio
Run file "Bicrystal_Rot_SS_script.m" 
The output will be stored in "Bicrystal_Rot_SS.txt" in CSV format.
The columns of the data are organized as stress (GPa) and strain for some rotation (degrees) of a bicrystal boundary: 
strain, theta=0, theta=5, theta=15, theta=25, theta=45, theta=55, theta=65, theta=70, theta=90

Fig. 6 Dissipation energy
Run file "DispEnergy.m"
The output will be stored in two files "DispEnergyinterp_data.txt" and "DispEnergyfit_data.txt" in CSV format
The columns of the data are organized as STGB angle (degrees) and dissipation energy (GPa) for both files.

Fig. 7 Twin Boundaries
Run file "TwinBoundaries.m"
The output will be stored in the file "TwinBoundariesSS.txt" along with intermediate files "GrainData_GB_NONeng_twins.mat" and "GrainData_GB_NONeng_random.mat". The intermediate files contain the microstructure information of both samples.
The columns of the data are organized as positive strain correlated to stress (GPa) for random boundaries and stress (GPa) for Twin boundaries, and negative strain correlated to stress (GPa) twin boundaries.
strain, stress (random), stress (twin), strain, stress (twin).

Fig. 8 EBSD 3 distributions
Run file "oriDist.m"
The output will be stored in file "oriDist.txt" along with intermediate files "GrainData_EBSD_GBE.mat", "GrainData_EBSD_GBESubet.mat", and "GrainData_EBSD_nonGBE.mat". The intermediate files contain the microstructure information of the samples.
The columns are organized as:
GBE EBSD data- Phi, Cumulative Frequency, phi1, Cumulative Frequency, phi2, Cumulative Frequency, GBE distribution data- Phi, Cumulative Frequency, phi1, Cumulative Frequency, phi2, Cumulative Frequency, non-GBE EBSD data- Phi, Cumulative Frequency, phi1, Cumulative Frequency, phi2, Cumulative Frequency, non-GBE distribution data- Phi, Cumulative Frequency, phi1, Cumulative Frequency, phi2, Cumulative Frequency,GBE subset EBSD data- Phi, Cumulative Frequency, phi1, Cumulative Frequency, phi2, Cumulative Frequency, GBE subset distribution data- Phi, Cumulative Frequency, phi1, Cumulative Frequency, phi2, Cumulative Frequency

Fig. 9 Volume flux


Fig. 10 EBSD stress strain
Run file "EBSDss.m"
The output will be stored in files "EBSD_stats.txt" and "EBSD_samples.txt" in CSV format. 
