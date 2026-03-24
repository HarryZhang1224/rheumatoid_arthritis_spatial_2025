Quantification of IF signal

Inputs:
- pSTAT3 panel CSV: /Users/ianmantel/Documents/PhD_Postdoc/Research/Experimental Log/IM104/final_images/viewer_parameter_exports/RA1_opn_cd14_fibrin_pSTAT3.figure_parameters_031426.csv
- ECM panel CSV: /Users/ianmantel/Documents/PhD_Postdoc/Research/Experimental Log/IM104/final_images/viewer_parameter_exports/RA1_opn_cd14_fibrin_CollagenI.figure_parameters_031426.csv
- OME search root: /Users/ianmantel/Documents/PhD_Postdoc/Research/Experimental Log/IM104/final_images/raw_ome_data

Outputs:
- if_signal_quantification_summary.csv
- pstat3_seed_calls.csv
- ecm_seed_calls.csv
- qc/*.mask_qc.png

Notes:
- DAPI is used to define nucleus-anchored events.
- CD14 and OPN are assigned in perinuclear neighborhoods.
- pSTAT3 is assigned in the nuclear neighborhood.
- Fibrin and CollagenI are assigned as local ECM-region signals.
- Optional coarse-mask thresholds can be provided per panel to suppress background artifact.
