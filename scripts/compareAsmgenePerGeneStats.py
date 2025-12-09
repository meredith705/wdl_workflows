#!/usr/bin/env python3

import os
import argparse
import pandas as pd

"""

python3 ~/Desktop/paten/wdl_git_workflows/wdl_workflows/scripts/compareAsmgenePerGeneStats.py 
    -a ../../AsmgeneStats/perGeneStats/AsmgenePerGeneStats/RUSH_027_FTX_hapdup_dual_1.0.97.per_gene_stats.txt 
    -b hifiasmAsmgenePerGeneStats/RUSH_027_FTX_hifiasm.ont.bp.hap1.p_ctg.0.97.per_gene_stats.txt 
    -o asmgenePerStats_RUSH_027_FTX.tsv

"""

def parse_asmgene_gene_stats(report_path):
    """
    Parse a asmgene gene_stats.txt file into a dictionary of metrics.
    """
    metrics = {}
    with open(report_path) as f:
        for line in f:
            # skip first line, though it has headers
            if line.strip() == "":
                continue

            headers = ['H', 'Metric', 'genesToRef', 'hapdup_dual_1']
            if line.startswith("H"):
                headers = line.strip().split()

            # Each line looks like: "X  full_sgl    34173   31413"
            # split on whitespace
            parts = line.strip().split()


            # asmgene has 4 columns  
            # key = " ".join(parts[:-2]) if len(parts) > 2 else parts[0]
            key = parts[1]

            # ref value
            ref_value = parts[2]
            # asm value
            asm_value = parts[3]
                   
            # ref metrics
            metrics['ref_'+key] = ref_value
            metrics['asm_'+key] = asm_value


    return metrics

def compare_asmgene_problems(singleCopy, multiCopy, output_file):
    """
    Of genes that had problematic output for each assembler which genes had similar problems
    """
    # open a file to store outputs and write header
    countsFile = open(f"comparison_counts_{output_file}", 'w')
    countsFile.write(f"same_sc_problem\tperctSame_sc_problem\ttotal_sc_problems\tdiff_sc_problem\tShasta_not_hifiasm_scGenes\thifiasm_not_shasta_scGenes\tsame_mc_problem\n")

    # single copy genes where gene problem matches. 
    same_sc_problem = singleCopy.loc[ (singleCopy['asmgeneVal_hapdup_dual_1']==singleCopy['asmgeneVal_hifiasm.ont.bp.hap1.p_ctg']) & (singleCopy['asmgeneVal_hapdup_dual_1']!="X") & (singleCopy['asmgeneVal_hapdup_dual_1']!="H")]
    print('same single copy problem:', same_sc_problem.shape[0], f'thats { round(100*(same_sc_problem.shape[0]/(singleCopy.shape[0]-28)),2)}% of {singleCopy.shape[0]-28}',
     '\n', same_sc_problem['asmgeneVal_hapdup_dual_1'].value_counts() )
    # write to counts file
    countsFile.write(f"{same_sc_problem.shape[0]}\t{ round(100*(same_sc_problem.shape[0]/(singleCopy.shape[0]-28)),2)}\t{singleCopy.shape[0]-28}\t")

    # this doesn't account for many of Shasta only problems.. I wonder if this is not including the rows where there is no value for hifiasm TODO
    dif_sc_problem = singleCopy.loc[ (singleCopy['asmgeneVal_hapdup_dual_1']!=singleCopy['asmgeneVal_hifiasm.ont.bp.hap1.p_ctg']) & (singleCopy['asmgeneVal_hapdup_dual_1']!="X") & (singleCopy['asmgeneVal_hapdup_dual_1']!="H")]
    dif_sc_problem = dif_sc_problem.dropna()
    print( 'ShastaHD problem counts:','\n', dif_sc_problem['asmgeneVal_hapdup_dual_1'].value_counts() )
    print('genes shasta gets that hifiasm does not (0/M)', dif_sc_problem.loc[ dif_sc_problem['asmgeneVal_hifiasm.ont.bp.hap1.p_ctg'].isin(["0", 'M'])][['asmgeneVal_hapdup_dual_1','asmgeneVal_hifiasm.ont.bp.hap1.p_ctg','chrom_hapdup_dual_1', 'start_hapdup_dual_1', 'end_hapdup_dual_1']])
    print( 'hifiasm problem counts:','\n', dif_sc_problem['asmgeneVal_hifiasm.ont.bp.hap1.p_ctg'].value_counts() )
    print('genes hifiasm gets that shasta does not (0/M)', dif_sc_problem.loc[ dif_sc_problem['asmgeneVal_hapdup_dual_1'].isin(["0", 'M'])][['asmgeneVal_hapdup_dual_1','asmgeneVal_hifiasm.ont.bp.hap1.p_ctg','chrom_hapdup_dual_1', 'start_hapdup_dual_1', 'end_hapdup_dual_1']])
    #write to counts
    countsFile.write(f"{dif_sc_problem['asmgeneVal_hapdup_dual_1'].shape[0]}\t{dif_sc_problem.loc[ dif_sc_problem['asmgeneVal_hifiasm.ont.bp.hap1.p_ctg'].isin(['0', 'M'])].shape[0]}\t{dif_sc_problem.loc[ dif_sc_problem['asmgeneVal_hapdup_dual_1'].isin(['0', 'M'])].shape[0]}\t")

    same_mc_problem = multiCopy.loc[ multiCopy['asmgeneVal_hapdup_dual_1']==multiCopy['asmgeneVal_hifiasm.ont.bp.hap1.p_ctg']]
    print('same multi copy problem:', same_mc_problem.shape )
    countsFile.write(f"{same_mc_problem.shape[0]}\n")
    countsFile.close()

def main(input_asmgene_a, input_asmgene_b, output_file):
    """
    load the 2 files and merge on gene id.
    """

    columnNames = ['asmgeneVal', 'assembly', 'geneID', 'length', 'chrom', 'start', 'end', 'dup_end']
    asmgene_a = pd.read_csv(input_asmgene_a, sep="\t", header=None, names=columnNames, engine="python")
    asmgene_b = pd.read_csv(input_asmgene_b, sep="\t", header=None, names=columnNames, engine="python")

    suffixA = asmgene_a.iloc[0,1]
    suffixB = asmgene_b.iloc[0,1]
    print('asm A', suffixA, 'asm B', suffixB)

    # print('whole df\n', asmgene_a.head())
    # print(suffixA, 'shape', asmgene_a.shape, suffixB, 'shape', asmgene_b.shape)

    # isolate rows describing duplicated genes
    asmgene_a_dup = asmgene_a.loc[asmgene_a['asmgeneVal']=='d']
    asmgene_b_dup = asmgene_b.loc[asmgene_b['asmgeneVal']=='d']
    # print('duplicated\n', asmgene_a_dup.head())
    # print(suffixA, ' dup shape', asmgene_a_dup.shape, suffixB, ' dup shape', asmgene_b_dup.shape)

    columnNamesDup = ['asmgeneVal', 'assembly', 'copyNumber', 'geneID', 'length', 'chrom', 'start', 'end']
    asmgene_a_dup.columns = columnNamesDup
    asmgene_b_dup.columns = columnNamesDup

    # isolate rows for single copy genes
    columnNames = ['asmgeneVal', 'assembly', 'geneID', 'length', 'chrom', 'start', 'end']
    asmgene_a_singleCopy = asmgene_a.loc[asmgene_a['asmgeneVal']!="d"]
    asmgene_b_singleCopy = asmgene_b.loc[asmgene_b['asmgeneVal']!="d"]
    
    # drop the extra column
    asmgene_a_singleCopy = asmgene_a_singleCopy.drop(columns=["dup_end"])
    asmgene_b_singleCopy = asmgene_b_singleCopy.drop(columns=["dup_end"])

    # print('single copy\n', asmgene_a_singleCopy.head())
    # print(suffixA, 'shape', asmgene_a_singleCopy.shape, suffixB, 'shape', asmgene_b_singleCopy.shape)

    
    asmgene_a_singleCopy.columns = columnNames
    asmgene_b_singleCopy.columns = columnNames

    

    df_single = pd.merge(asmgene_a_singleCopy, asmgene_b_singleCopy, on="geneID", how="outer", suffixes=("_"+suffixA,"_"+suffixB))
    df_double = pd.merge(asmgene_a_dup, asmgene_b_dup, on="geneID", how="outer", suffixes=("_"+suffixA,"_"+suffixB))



    # Save to file
    df_single.to_csv("singleCopy_"+output_file, sep="\t")
    df_double.to_csv("multiCopy_"+output_file, sep="\t")
    print(f"Compiled report written to singleCopy_{output_file} and multiCopy_{output_file}")

    # which genes have shared problems
    compare_asmgene_problems(df_single, df_double, output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare genes assembled pairwise from asmgene per gene stats files")
    parser.add_argument(
        "-a", "--input_asmgene_a",
        required=True,
        help="file containing asmgene per-stats assembler A"
    )

    parser.add_argument(
        "-b", "--input_asmgene_b",
        required=True,
        help="file containing asmgene per-stats assembler B"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file (TSV format)"
    )
    args = parser.parse_args()

    main(args.input_asmgene_a, args.input_asmgene_b, args.output)


'''
Asmgene output per gene stats
https://github.com/lh3/minimap2/blob/6d49eb690f3c32ae2b95a796951397bf598396f0/misc/paftools.js#L954

output here: https://github.com/lh3/minimap2/blob/6d49eb690f3c32ae2b95a796951397bf598396f0/misc/paftools.js#L968

single copy genes are not in the pergene stats output (gene[g][i][0] == 1)

type.   Input           Gene (QueryName)    QuerySeqLen?    targetSeqName targetStart   targetEnd
                                            ExpectedGeneLen 
D       hapdup_dual_1   ENST00000620569.1       409     chr7    143917290       143917799                gene[g][i][0] > 1 Duplicated? 
M       hapdup_dual_1   ENST00000553483.3       931     chr14   88593831        88596213                 unmapped (M) // missing completely
F       hapdup_dual_1   ENST00000402588.6       937     chr22   24387464        24391117                 gene[g][i][1] >= opt.min_cov     


| Code    | Meaning                           | Trigger condition              |
| ------- | --------------------------------- | ------------------------------ |
| M       | Missing gene                      | `gene[g][i] == null`           |
| D       | Duplicated gene                   | `gene[g][i][0] > 1`            |
| F       | Fragmented gene (meets min_cov)   | `gene[g][i][1] >= opt.min_cov` |
| 5       | Partial (≥50% coverage)           | `gene[g][i][1] >= 0.5`         |
| 1       | Partial (≥10% but <50%)           | `gene[g][i][1] >= 0.1`         |
| 0       | Absent (<10% coverage)            | otherwise                      |



# convert paf to sortened storage 
while (file.readline(buf) >= 0) {
            var t = buf.toString().split("\t");
            var ql = parseInt(t[1]), 
                qs = parseInt(t[2]), 
                qe = parseInt(t[3]), 
                mlen = parseInt(t[9]), 
                blen = parseInt(t[10]), 
                mapq = parseInt(t[11]);
            if (i == getopt.ind) refpos[t[0]] = [t[0], t[1], t[5], t[7], t[8]];
            if (gene[t[0]] == null) gene[t[0]] = [];
            if (a.length && t[0] != a[0][0]) {
                gene[a[0][0]][i - getopt.ind] = process_query(opt, a);
                a = [];
            }
            a.push([t[0], ql, qs, qe, mlen, blen]);
        }

paf stored as:
a[j][1] → query length (QLEN)
a[j][2] → query start (QS)
a[j][3] → query end (QE)
a[j][4] → number of matching bases
a[j][5] → alignment block length

function process_query(opt, a) {
        # b stores alignments
        var b = [], cnt = [0, 0, 0];

        for (var j = 0; j < a.length; ++j) {
            # identity thresholding cutoff
            # matches / aligned_length >= min_iden
            #if #_matching bases < (alignment lenghth * min_iden) skip alignment
            if (a[j][4] < a[j][5] * opt.min_iden)
                continue;
            # store alignments that pass iden threshold
            b.push(a[j].slice(0));
        }
        
        # if the matches is less than min_iden return 0's
        if (b.length == 0) return cnt;

        // count full length transcripts
        # innitialize full count
        var n_full = 0;
        # loop through all the alignments and count up those that cover min_cov
        for (var j = 0; j < b.length; ++j)
            #if the query_end - query_start >= (query_length * min_cov) 
            if (b[j][3] - b[j][2] >= b[j][1] * opt.min_cov)
                # count it as full
                ++n_full;
        # store counts in the count array 0ith position
        cnt[0] = n_full;

        // compute coverage
        b = b.sort(function(x, y) { return x[2] - y[2] });
        var l_cov = 0, st = b[0][2], en = b[0][3];
        for (var j = 1; j < b.length; ++j) {
            if (b[j][2] <= en)
                en = b[j][3] > en? b[j][3] : en;
            else l_cov += en - st;
        }
        l_cov += en - st;
        # coverage fraction
        cnt[1] = l_cov / b[0][1];
        # number of alignments
        cnt[2] = b.length;
        return cnt;
        # cnt = [ n_full, cov_frac, n_aln ]
    }


    // count and print
    var col1 = ["full_sgl", "full_dup", "frag", "part50+", "part10+", "part10-", "dup_cnt", "dup_sum"];
    # innitialize rst
    var rst = [];
    for (var k = 0; k < col1.length; ++k) {
        rst[k] = [];
        for (var i = 0; i < n_fn; ++i)
            rst[k][i] = 0;
    }

    for (var g in gene) { // count single-copy genes
        if (gene[g][0] == null || gene[g][0][0] != 1) continue;          // must be single-copy gene in ref
        if (gene_nr[g] == null) continue;                                // skip if no reference info
        if (auto_only && /^(chr)?[XY]$/.test(refpos[g][2])) continue;    // skip sex chromosomes if requested

        for (var i = 0; i < n_fn; ++i) {
            if (gene[g][i] == null) {
                rst[5][i]++;                                                   // missing completely
                if (print_err) print('M', header[i], refpos[g].join("\t"));
            } else if (gene[g][i][0] == 1) {
                rst[0][i]++;                                                   // full single-copy
            } else if (gene[g][i][0] > 1) {
                rst[1][i]++;
                if (print_err) print('D', header[i], refpos[g].join("\t"));    // duplicated gene
            } else if (gene[g][i][1] >= opt.min_cov) {
                rst[2][i]++;
                if (print_err) print('F', header[i], refpos[g].join("\t"));    // fragmented (meets min_cov threshold)
            } else if (gene[g][i][1] >= 0.5) {
                rst[3][i]++;
                if (print_err) print('5', header[i], refpos[g].join("\t"));    // partially covered (≥50%)
            } else if (gene[g][i][1] >= 0.1) {
                rst[4][i]++;
                if (print_err) print('1', header[i], refpos[g].join("\t"));    // weakly covered (≥10%)
            } else {
                rst[5][i]++;                                                    
                if (print_err) print('0', header[i], refpos[g].join("\t"));    // not found (<10%)
            }
        }
    }

    for (var g in gene) { // count multi-copy genes
        if (gene[g][0] == null || gene[g][0][0] <= 1) continue; // not in ref or single copy, skip
        if (gene_nr[g] == null) continue;                       // skip genes with no ref or autosomal optional
        if (auto_only && /^(chr)?[XY]$/.test(refpos[g][2])) continue;
        

        for (var i = 0; i < n_fn; ++i) {
            if (gene[g][i] != null) rst[7][i] += gene[g][i][0];    // if gene in asm add it to num copies found dup_sum
            if (gene[g][i] != null && gene[g][i][0] > 1) {
                rst[6][i]++;                                       // if there ismore than one copy increment the number of multicopy genes dup_count
            } else if (print_err) {
                print('d', header[i], gene[g][0][0], refpos[g].join("\t"));  // print out the missing or collapsed multigenes, and expected ref copy number
            }
        }
    }
'''