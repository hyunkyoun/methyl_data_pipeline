# methyl_data_pipeline


This is done for each combat normalization separately. 
parse sampleTable txt files -> run R script -> check for Sentrix Barcode + Sample Section in the samplesTable -> if same one is found, replace the header with the sampleID (if not, remove the column) 

After this, 
run combat processing for each respective runs. 