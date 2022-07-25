# Deseq2 Bulk RNA-seq Report

The aim of this Rmarkdown website is to automatically run Clustering, Deseq2 and Enrichment analysis using bulk RNA-seq data.

Documentation for this code can be accessed [here](https://acribbs.github.io/deseq2_report/)

Ths report requires either:
1. a counts table (samples x genes) generated from running running mapping (hisat2) followed by featurecounts (using this [pipeline](https://github.com/cgat-developers/cgat-flow/blob/master/cgatpipelines/tools/pipeline_rnaseqdiffexpression.py)) or
2. a directory containing  kallisto quant output (using this [pipeline](https://github.com/Acribbs/cribbslab/blob/master/cribbslab/pipeline_pseudobulk.py)).

To run the analysis workflow:

* Clone the repository:

`git clone https://github.com/Acribbs/deseq2_report.git`


## Running the Rmarkdown using featurecounts output

In order to run the report you require the following input files for the report to generate a report correctly:

1. A meta data file following the naming convention design_<test>_<control>_<test>_<column>.csv
2. A counts table called featurecounts.tsv.gz
3. A config.yml for modifying the output of the deaseq_report
4. A full.csv file listing all the input files

Example files are currently within the repo so if you are unsure on how each of these files should be layed out please refer to the repo.

* Navigate to the directory and rename the deseq2.Rproj file to something of your choosing

* modify the config.yml file

* Remove the design_* files and make your own based on the following naming convention:

  `design_<test>_<control>_<test>_<column>.csv`
1. `<test>` - refers to the test that you plan to run. There are two options "ltr" or "wald".
2. `<control>` - This is the name of your control condition i.e. the samples you want to test against. This should match on of the samples in the <column> of the file.
3. `<test>` - This is the name of the test condition. This should match with one of the samples in the <column> of the file
4. `<column>` - This is a column that you want to use for your Deseq model in the design_* file

You can have multiple design_* files in the folder.

* Double click the Rproj folder and the project should open in Rstudio then click "Build Website"

![Location of Build Website in Rstudio](https://raw.githubusercontent.com/Acribbs/deseq2_report/master/img/build_img.png)

Make sure rmarkdown is installed in your library and hit the build tab on the environment window and then click "Build Website". When the website has finished building a window will pop up with the rendered site. The final report is in the directory "Final_report" and can be accessed by opening the index.html file in a web browser.

## Running pipeline using kallisto pipeline output

In order to run the report you require the following input files for the report to generate a report correctly:

1. A meta data file following the naming convention design_<test>_<control>_<test>_<column>.csv
2. A kallisto quant directory containing the output from kallisto
3. A config.yml for modifying the output of the deaseq_report
4. A kallisto_input.csv file listing all the files


Example files are currently within the repo so if you are unsure on how each of these files should be layed out please refer to the repo.

* Navigate to the directory and rename the deseq2.Rproj file to something of your choosing

* modify the config.yml file

* Remove the design_* files and make your own based on the following naming convention:

  `design_<test>_<control>_<test>_<column>.csv`
1. `<test>` - refers to the test that you plan to run. There are two options "ltr" or "wald".
2. `<control>` - This is the name of your control condition i.e. the samples you want to test against. This should match on of the samples in the <column> of the file.
3. `<test>` - This is the name of the test condition. This should match with one of the samples in the <column> of the file
4. `<column>` - This is a column that you want to use for your Deseq model in the design_* file

You can have multiple design_* files in the folder.

* Double click the Rproj folder and the project should open in Rstudio then click "Build Website"

![Location of Build Website in Rstudio](https://raw.githubusercontent.com/Acribbs/deseq2_report/master/img/build_img.png)

Make sure rmarkdown is installed in your library and hit the build tab on the environment window and then click "Build Website". When the website has finished building a window will pop up with the rendered site. The final report is in the directory "Final_report" and can be accessed by opening the index.html file in a web browser.
