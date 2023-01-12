# Metadata processing

After downloading raw metadata from NCBI for contigs, we process it to get a clean metadata file. 

Process involves:

* run `01_process-metadata.R`: `metadata-raw.csv` > `metadata-processed-v1.csv`
* Manually edit spreadsheet: `metadata-processed-v1.csv` > `metadata-processed-v2.csv`
* run `02_process-metadata.R`: `metadata-processed-v2.csv` > `metadata-processed-v3.csv`

That final `v3` file becomes `data/metadata.csv` 

