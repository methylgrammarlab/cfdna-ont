# cfdna-ont

## Preparation 

To deconvolute the samples to cell types, the needed file foramt is bedGraph, with the columns: chr, start, end, methylation_level

The bedGraph file needs to be mapped to the 450k illumina human methylation array, using a relevant manifest file, that map coordinate to the relevant probes.
The file for the analysis should contain the columns:
1. chr
2. start
3. end
4. illumina's probe name
5. chr
6. start
7. end
8. methylation value

## Run Moss deconvolution

To run moss deconvolution on the samples:
First, you should run the [merge_all_samples_for_deconvolusion.R](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/deconvolution_moss/merge_all_samples_for_deconvolusion.R),to merge all samples to be on one table.

Second, you should take the a reference atlas (you can use the original atlas that was used in [Moss paper](https://www.nature.com/articles/s41467-018-07466-6#Sec13) [](https://github.com/nloyfer/meth_atlas/blob/master/reference_atlas.csv), and run deconvolution as explained in the [usage page](https://github.com/nloyfer/meth_atlas#usage)
This will also create a plot of the deconvolution.

To plot the deconvolution results with the editing, you can use the script [plot_deconv.py](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/deconvolution_moss/plot_deconv.py)

To plot the deconvolution results with redued cell types, you can use the script [deconvolution_plot.R](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/deconvolution_moss/deconvolution_plot.R), with or without a group file. An example for a group file is [here](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/deconvolution_moss/group_file_for_plot_green_epithilial.csv)

## Run 2-component deconvolution

To run the deconvolution, separating the samples only to lung component and healthy component, I used pure lung epitelial probes, and TCGA LUAD probes.

For the pure lung epithelial cells, I used the script: [](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/cell_type_probes/plot_score_of_methylation_according_to_nanopore_sampeles_create_table.R) to create a table with the nnls results

For the TCGA LUAD probes, I first corrected the methylation values using the script [](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/TCGA_probes/correct_methylation_values_by_tumor_purity_create_table.R)

Then, I run the script [](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/TCGA_probes/plot_score_of_methylation_according_to_nanopore_sampeles_create_table.R) to create a table with the nnls results


To plot the results, I used the script: [plot_tumor_fractions_vs_score.R](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/script_for_plot_probes/plot_tumor_fractions_vs_score.R)
