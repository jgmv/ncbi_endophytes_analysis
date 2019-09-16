# ncbi_endophytes_analysis
Analysis of NCBI records of fungal endophyte ITS sequences.

Retrieves and analyzes ITS sequence data and metadata of fungal endophytes from NCBI GenBank. Checks how records represent different fungal taxa or plant hosts, and how they are distributed geographically and across plant organs.

## Steps
### Fecth data from NCBI GenBank
Run script `fetch_data.sh`. This fetches records in GenBank format and extracts sequence data and metadata. Then, it identifies *de novo* the sequences using the Naïve Bayesian Classifier of `mothur` (system-wide installation required), using UNITE fungal ITS reference sequences (downloaded automatically). 
Fungi other than endophytes can be targeted by changing the search query in line 35 of `fetch_data.sh`. 

The files resulting from this step are provided in folder `files`, so that this part can be skipped to proceed with the data analysis.

### Analyze records data
Run R script `R_script.R`. Requires installation of the following packages: `ggplot2`, `ggthemes`, `maps`, `maptools`, `plotrix`, `RColorBrewer`, `rgeos`, `rnaturalearth`, `sf`, `sp`, and `Taxonstand`.
Output files are stored in folder `output`.
