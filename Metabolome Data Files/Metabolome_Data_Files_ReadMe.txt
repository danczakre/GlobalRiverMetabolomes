# ReadMe for the metabolome data files
'S19S_Sed_Water_08.12_Report.csv' is the raw, aligned FTICR-MS report generated using
Formularity (Tolic et al., 2017). Raw data can be obtained via ESS-DIVE. This file was
compressed to save space and will need to be unzipped before use.

'S19S_Sed-Water_08.12_Poorly_Calibrated_Samples.csv' contains all of the samples that did
not calibrate well enough in Formularity to be included in our analysis.

'Processed_S19S_Sed-Water_08.12_Data.csv' and 'Processed_S19S_Sed-Water_08.12_Mol.csv' are files
obtained by using the R package ftmsRanalysis (https://github.com/EMSL-Computing/ftmsRanalysis)
and are the processed versions of the raw FTICR-MS report. These files have been processed
to include non-isotopic metabolites with a mass between 200-900 m/z. These files were
compressed to save space and will need to be unzipped before use.

'S19S_Sed-Water_08.12_Trans_Profiles.csv' contains the counts of transformations identified
in each of our samples. This file was obtained by running the 'Data' file through the
transformation analysis found at https://github.com/danczakre/Meta-Metabolome_Ecology.

'S19S_Metadata.csv' contains meta-data information (i.e., latitude, longitude, site, etc.)