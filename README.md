# Aedes aegypti microbiota Read Classification From WGS

These scripts are for my third rotation project where I classified the Aedes aegypti microbiota using whole genome sequences. The premise behind this project was to use the unaligned reads for read classification. To run the analyses, clone the the scripts located in the scripts folder. Run the 'create folders' script. This will create three folders - "Working" "Results" "Data". From the working folder, run all the required scripts in the desired order to re-run the analysis pipeline.

To re-perform the analyses the following steps must be taken
1. Load the scripts into your working directory
2. run the create_folders command to create the neccessary folders for analyses. This command will create three folders; data, results and working. The working directory is where you want to be for all the subsequent analysis steps. 
3. Download the desired paried WGS for Aedes aegypti from the ENA browser website for which you want to classify the microbiota
4. Download the desired reference to map the reads against. I would recommend the recent Aedes aegypti reference genome from Vectorbase located here: https://vectorbase.org/vectorbase/app/record/dataset/TMPTX_aaegLVP_AGWG
5. Make sure BOTH of these files are situated in the data wd with the aedes reference in the aedes_ref folder and fastq files in fastq folder (I would recommend using enaBrowser tools to pull the sequences directly)
6. Change to the working folder and type create_samples.txt. This command will create a .txt file with your sample accession numbers.
7. Now run all all anlaysis steps using the given scripts. You will be able to run the whole analysis from the working directory folder and .txt file! 

(I understand this isn't the most efficent way to recapitulate the analysis...still a learning process). Any questions please email 248064@lstmed.ac.uk

