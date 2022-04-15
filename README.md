# cfdna-ont

## Preparation 

To deconvolute the samples to cell types, the data must be in the bedGraph format containing  columns: chr, start, end, methylation_level.
The bedGraph file must then be mapped to the 450k Illumina human methylation array, using a relevant manifest file that maps coordinates to the relevant probes. The file for the analysis should contain the columns:
1. chr
2. start
3. end
4. illumina's probe name
5. chr
6. start
7. end
8. methylation value


## Run Moss Deconvolution

To run Moss deconvolution on the samples: First, run [merge_all_samples_for_deconvolusion.R](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/deconvolution_moss/merge_all_samples_for_deconvolusion.R) in order to merge all samples into one table.

Next, provide a reference atlas (this can be the original atlas that was used in the [Moss paper](https://www.nature.com/articles/s41467-018-07466-6#Sec13) [](https://github.com/nloyfer/meth_atlas/blob/master/reference_atlas.csv), and run the deconvolution as explained in the [usage page](https://github.com/nloyfer/meth_atlas#usage). This will also create a plot of the deconvolution.

The script [plot_deconv.py](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/deconvolution_moss/plot_deconv.py) will plot the deconvolution results with some editing options different from the original. To plot the deconvolution results with fewer cell types, the script [deconvolution_plot.R](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/deconvolution_moss/deconvolution_plot.R) can be used, with or without a file with sample groups. An example group file is found [here](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/deconvolution_moss/group_file_for_plot_green_epithilial.csv).


## Run 2-component Deconvolution

To deconvolute the samples into the lung component and the healthy component, only pure lung epithelial and TCGA LUAD probes were used.

For the pure lung epithelial cells, the [](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/cell_type_probes/plot_score_of_methylation_according_to_nanopore_sampeles_create_table.R) was used to create a table with the NNLS results

For the TCGA LUAD probes, the [](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/TCGA_probes/correct_methylation_values_by_tumor_purity_create_table.R) was used to adjust the methylation values according to the tumor purity. The [](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/TCGA_probes/plot_score_of_methylation_according_to_nanopore_sampeles_create_table.R) script was then used to create a table with the NNLS results.
To plot the results, use the script: [plot_tumor_fractions_vs_score.R](https://github.com/methylgrammarlab/cfdna-ont/blob/main/deconvolution_code/script_for_plot_probes/plot_tumor_fractions_vs_score.R)

## CTCF nucleosome positioning Heatmaps

To produce the CTCF heatmaps we used the 9,780 evolutionarily conserved CTCF motifs occurring in distal ChIP-seq peaks. 
To get the coverage for the regions for each sample [deepTools](https://deeptools.readthedocs.io/en/develop/) (Version 3.5.0).

[bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) was used with the parameters `--ignoreDuplicates --binSize -bl ENCODE_blacklist -of bedgraph --effectiveGenomeSize 2913022398 --normalizeUsing RPGC`. For Illumina WGS the additional parameter `--extendedReads 145` was used. The bedgraph was converted to a bigwig file using bigWigToBedGraph downloaded from UCSC Genome Browser. This bigwig file was passed to deepTools computeMatrix with the command line parameters `reference-point --referencePoint center -out table.out`, and the table.out file was imported into R to create fragment coverage heatmap using the [](https://github.com/methylgrammarlab/cfdna-ont/blob/main/heatmapOfCTCFrerun.R) script.
For the tumor fraction annotation on the side of the plot, a table with the tumor fraction is needed (in this case ichorCNA was used to obtain the tumor fraction. 
