# Carfilzomib induced cardiotoxicity

All data from this project is accesible under the following link: https://figshare.com/collections/Carfilzomib_mediated_cardiotoxicity/6716034. 
1) Cell-barcode-ID-to-clusters from all scRNAseq datasets are available as .tsv files.
2) Seurat files from scRNAseq analysis are accesible as .rds files. 
3) The R code for the analysis of the murine ECG as well as example files are provided as the .zip file "CardiotoxicityCFZ". A Shinny example can be found in: https://maxmendezl.shinyapps.io/ShinnyApp/


# April 2026
Major Updates:
- Refactored to improve robustness and maintainability
- Peak detection now uses QRS-anchored workflo instead of local maxima as before. This improves P-wave identification across experimental conditions.
- Data loading is safer, added validation for columns and fallback handling.
- App structure was simplified to reduce code duplication.
- Now the App calculates the ECG intervasl directly from the input.
- New helper added for signal smoothing, polarity aware QRS/R detection, beginn and end of wave, live calculation of intervals.

- Remaining issues (for future updates)
  -- QRS width is still simplified: in the future update a QRS width calculator will be added.
  -- Condition-specific tuning still needs refinement, specially for difficult traces.
  -- Debug message will be added
  -- complete migration of all modules to ecg_data and merge all plots and overlapping modules. 
