# Shasta Evaluation WDL Workflows

WDL workflows that automates runing GFAse, yak, dipcall, and whatshap on shasta assembly small bubble gfa files. 

# WDL Cromwell Dependencies 
Download cromwell tools to run wdl worflows locally:
```
wget https://github.com/broadinstitute/cromwell/releases/download/71/cromwell-71.jar
wget https://github.com/broadinstitute/cromwell/releases/download/71/womtool-71.jar
```
Install jdk
```
sudo apt-get update
sudo apt install openjdk-17-jre-headless
```
Install docker ( directions from here https://docs.docker.com/engine/install/ubuntu/ )
```
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
```

# Installation
Clone this git repo to obtain the wdl workflows:
```
git clone https://github.com/meredith705/wdl_workflows.git
```

# WDL Input Files
Download the yak parental count files, phased GIAB_hg002_v_hg38_truthset.vcf, and hg38 assembly. These are downloaded outside of the WDL and are re-used in each WDL execution. 
The files and paths in inputs.json can be updated as needed. The inputs.json file assumes the input data is a subdirectory of the shasta_eval working directory where the wdl will be executed.
```
cd wdl_workflows/shasta_eval
mkdir input_data
cd input_data
wget http://public.gi.ucsc.edu/~memeredith/hg002_shasta_eval_wdl_files/hg03.ilmn.k31.pe.yak && \
wget http://public.gi.ucsc.edu/~memeredith/hg002_shasta_eval_wdl_files/hg04.ilmn.k31.pe.yak && \
wget http://public.gi.ucsc.edu/~memeredith/hg002_shasta_eval_wdl_files/hg02.ilmn250.k31.pe.yak && \
wget http://public.gi.ucsc.edu/~memeredith/hg002_shasta_eval_wdl_files/hg03.48_55.unique.k31.fa && \
wget http://public.gi.ucsc.edu/~memeredith/hg002_shasta_eval_wdl_files/hg04.54_61.unique.k31.fa && \
wget http://public.gi.ucsc.edu/~memeredith/hg002_shasta_eval_wdl_files/GIAB_HG002_GRCh38_1_22_v4.2.1_phased.vcf && \
wget http://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
# check into adding the medically relevant genes vcf
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz && \
cd ../
```
# Usage
Edit the `inputs.json` file:
1) specify the path to your shasta assembly file: `Assembly-Detailed.gfa`
2) Name the `'assemblyRunID'` as you see fit
3) Change any other input files or add other options.
To see all the possible inputs for the wdl run womtool
```
sudo java -jar bin/womtool-71.jar inputs gfase_yak_dipcall_whatshap.wdl
```

Run the wdl workflow:
```
sudo java -jar bin/cromwell-71.jar run -i inputs.json gfase_yak_dipcall_whatshap.wdl
```

# Output
wdl output is found in `cromwell-executions` followed by the apporpriate run id (###).
If `'assemblyRunID'` is 'chr11':
```
cat cromwell-executions/GFAseYakDipcallWhatshap/###/call-makeSummaryStatFile/execution/chr11.summary.txt
```

To look at yak and dipcall/whatshap output individually:
note: this output is for the shasta chr11 assembly (run 15)
```
# yak summary
cat cromwell-executions/GFAseYakDipcallWhatshap/###/call-yakAssemblyStats/execution/paternal.summary.txt 
# mat qv
FR      0.000994        0.00149
ER      132200147       826381.441
CV      0.049
QV      37.085  36.941
# pat qv
FR      0.00104 0.00154
ER      133187440       857064.852
CV      0.049
QV      36.975  36.814
# mat switch
W       10713   65621   0.163256    # switch error rate
H       8674    65636   0.132153    # hamming rate
N       8734    56904   0.133063
# pat switch
W       10695   64596   0.165568    # switch error rate
H       9018    64617   0.139561    # hamming rate
N       55370   9249    0.143131

# dipcall/whatshap summary
awk '{print $2,$33,$37}' cromwell-executions/GFAseYakDipcallWhatshap/###/call-coalesceResults/execution/sample.full.tsv
chromosome all_switch_rate blockwise_hamming_rate
chr11 0.00042188159189987344 0.00021093338020741783
```
