# Genomic features and GIE frequency (and simulated GIE frequency as control)

1. Run the 1.preprocess_variables.ipynb to extract most of genomic features information to later paerform the regression with the GIE (and simulated GIE) frequency
2. Run the 2.drivers_and_gie.ipynb to extract the driver information per sample and perform the analysis of association between drivers and GIE frequency (simulated GIE respectively)
3. Run the 3.perform_analysis_discrete_variables.ipynb notebook, to analyze the associated between discrete genomic features and GIE events (as well as control for simulated GIE)
4. Run the 3.perform_regression_continue_variables.ipynb notebook, to analyze the associated between continuos genomic features and GIE events (as well as control for simulated GIE)

5. Visualize the results displayed in Figure 7 and Supp. Fig. 4 using the visualize_enrichments.ipynb notebook 