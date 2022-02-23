# wdl_workflows

WDL workflows that automates runing GFAse, yak, dipcall, and whatshap on a shasta assembly small bubble gfa. 

Install cromwell tools to run wdl worflows locally:
```
cd bin
wget https://github.com/broadinstitute/cromwell/releases/download/71/cromwell-71.jar
wget https://github.com/broadinstitute/cromwell/releases/download/71/womtool-71.jar
```
Clone this git repo to obtain the wdl workflows:
```
git clone https://github.com/meredith705/wdl_workflows.git
```
Download the yak parental count files, phased GIAB hg002 v hg38 truthset vcf, hg38 assembly, and wdl inputs.json. These are downloaded outside of the WDL because they are used over and over for each WDL execution. 
The inputs.json can be chaned as needed it is set up so that the input data is a subdirectory of the working directory where the wdl will be executed.
```
cd shasta_eval
mkdir input_data
cd input_data
wget http://public.gi.ucsc.edu/~memeredith/hg002_shasta_eval_wdl_files/hg03.hiseq.k31.pe.yak && \
wget http://public.gi.ucsc.edu/~memeredith/hg002_shasta_eval_wdl_files/hg04.hiseq.k31.pe.yak && \
wget http://public.gi.ucsc.edu/~memeredith/hg002_shasta_eval_wdl_files/GIAB_HG002_GRCh38_1_22_v4.2.1_phased.vcf && \
# check into adding the medically relevant genes vcf
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz && \
wget http://public.gi.ucsc.edu/~memeredith/hg002_shasta_eval_wdl_files/inputs.json
cd ../
```

Run the wdl workflow:
```
java -jar bin/cromwell-71.jar run -i inputs.json gfase_yak_dipcall_whatshap.wdl
```

