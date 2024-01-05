#!/usr/bin/env python3
from multiprocessing import Pool
from matplotlib.lines import Line2D
from matplotlib import pyplot
from random import randint
import matplotlib
import subprocess
import argparse
import numpy
import sys
import os


matplotlib.use('Agg')


# Using https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
# and
HG38_ARM_LENGTHS = [
    123400000,
    125556422,
    39800000,
    93997422,
    53400000,
    81686622,
    35500000,
    97775309,
    17700000,
    96664328,
    17200000,
    89843718,
    19000000,
    82991189,
    36800000,
    53538345,
    25100000,
    58157441,
    18500000,
    61873285,
    26200000,
    32417616,
    93900000,
    148293529,
    28100000,
    36344167,
    12000000,
    34709983,
    15000000,
    35818468,
    90900000,
    107395559,
    50000000,
    140214555,
    48800000,
    132738259,
    59800000,
    111005979,
    60100000,
    99245973,
    45200000,
    99938636,
    43000000,
    95394717,
    61000000,
    95040895,
    10400000,
    46827415,
]


CHM13_ARM_LENGTHS = [
    124048267,
    124339061,
    40649191,
    94108943,
    52743313,
    82384456,
    35911664,
    97412884,
    16522942,
    97043744,
    11400261,
    89761231,
    17186630,
    82566565,
    36838903,
    59491471,
    25689679,
    58587218,
    18449624,
    62092914,
    27792923,
    33914441,
    93503283,
    149193469,
    28012753,
    38197502,
    11134529,
    33956153,
    14249622,
    37075304,
    94076514,
    107029434,
    52452474,
    141122471,
    48317879,
    133727560,
    59672548,
    112454080,
    62064435,
    98502993,
    45270456,
    100988875,
    46267185,
    104350062,
    59373565,
    94886001,
    10724418,
    51735611
]


PRESET_LENGTHS = {
    "hg38_chromosome_arms":HG38_ARM_LENGTHS,
    "chm13_chromosome_arms":CHM13_ARM_LENGTHS,
}

VALID_FILE_TYPES = {".fa",".fasta",".fna",".fai"}

def get_name_of_unit(n):
    if n == 1_000_000_000:
        return "Gbp"
    if n == 1_000_000:
        return "Mbp"
    if n == 1_000:
        return "Kbp"
    if n == 0:
        return "bp"


def plot_ngx(fig, axes, name, color, lengths, genome_size, output_dir):
    fig.set_size_inches(12,10)

    total = float(sum(lengths))
    genome_size = float(genome_size)

    if total > genome_size:
        sys.stderr.write("WARNING: observed total sequence length is greater than genome size\n")

    unit = float(10)**float(max(0,round((numpy.log10(genome_size) - 3) / 3)*3))

    ng50 = 0.0
    n50 = 0.0

    x = list()
    y = list()

    x_prev = 0.0
    s_prev = 0.0
    for i,s in enumerate(sorted(lengths, reverse=True)):
        s = float(s)/float(unit)

        x0 = x_prev
        x1 = x_prev + s

        if x0*unit <= float(genome_size)/2.0 < x1*unit:
            ng50 = s*unit
            print("ng50\t%.0f" % ng50)

        if x0*unit <= float(total)/2.0 < x1*unit:
            n50 = s*unit
            print("n50\t%.0f" % n50)

        if i > 0:
            # pyplot.plot([x0,x0],[s_prev,s], color=color)
            x.extend([x0,x0])
            y.extend([s_prev,s])

        x.extend([x0,x1])
        y.extend([s,s])

        x_prev = x1
        s_prev = s

    # Bring down to 0 for easier comparison
    x.extend([x_prev,x_prev])
    y.extend([s_prev,0])

    pyplot.plot(x,y,color=color)

    axes.axvline(genome_size/2/unit, linestyle='--', linewidth=0.6, color='gray')
    axes.ticklabel_format(style='plain')


    axes.set_xlim([0,float(genome_size)/float(unit)])

    axes.set_xlabel("Cumulative coverage (%s)" % get_name_of_unit(unit))
    axes.set_ylabel("Length (%s)" % get_name_of_unit(unit))
    axes.set_title("NGx")

    all_lengths_txt_path = os.path.join(output_dir, name + "_all_lengths.txt")
    ng50_txt_path = os.path.join(output_dir, "ng50.txt")
    n50_txt_path = os.path.join(output_dir, "n50.txt")

    with open(all_lengths_txt_path,'a') as file:
        for l in lengths:
            file.write(str(l))
            file.write('\n')

    with open(ng50_txt_path,'a') as file:
        file.write("%s,%.0f\n" % (name,ng50))

    with open(n50_txt_path,'a') as file:
        file.write("%s,%.0f\n" % (name,n50))


def index_fasta(path):
    index_path = path + ".fai"

    if not os.path.exists(index_path):
        command = ["samtools", "faidx", path]

        print("Running: " + ' '.join(command))

        result = subprocess.run(command, check=True, stderr=sys.stderr, universal_newlines=True)
    else:
        sys.stderr.write("Found existing fai: " + index_path + "\n")
        sys.stderr.flush()

    return True


def get_lengths_from_fasta(path):
    index_path = path + ".fai"

    lengths = list()

    # This won't do anything if the files are indexed already
    index_fasta(path)

    with open(index_path, 'r') as file:
        for l,line in enumerate(file):
            data = line.strip().split()
            lengths.append(int(data[1]))

    return lengths


def main(input_paths, output_dir, genome_size, color_indexes, n_threads, legends):
    if os.path.exists(output_dir):
        exit("ERROR: output directory already exists")
    else:
        os.makedirs(output_dir)

    colors = [
        (175/256.0,   48/256.0,   51/256.0),        #0 red
        (224/256.0,   99/256.0,   58/256.0),        #1 orange
        (215/256.0,   219/256.0,  84/256.0),        #2 yellow
        (110/256.0,   170/256.0,  100/256.0),       #3 light green
        (80/256.0,    180/256.0,  150/256.0),       #4 green
        (100/256.0,   189/256.0,  197/256.0),       #5 green-blue
        (0/256.0,     170/256.0,  231/256.0),       #6 turquoise
        (51/256.0,    87/256.0,   182/256.0),       #7 blue
        (37/256.0,    36/256.0,   93/256.0),        #8 indigo
        (95/256.0,    51/256.0,   139/256.0),       #9 purple
        (200/256.0,   53/256.0,   93/256.0)         #10 pink
    ]

    fig = pyplot.figure()

    axes = pyplot.axes()

    names = list()
    label_colors = list()
    max_len = 0

    sys.stderr.write("Indexing fastas...\n")

    fasta_paths = [p for p in input_paths if os.path.splitext(p)[-1] in VALID_FILE_TYPES]

    # if len(fasta_paths)>len(colors):
    #     colors = []
    #
    #     for i in range(len(fasta_paths)):
    #         colors.append('#%06X' % randint(0, 0xFFFFFF))

    pool = Pool(n_threads)
    success = pool.imap(index_fasta, fasta_paths, 1)
    pool.close()
    pool.join()

    sys.stderr.write("Plotting...\n")
    print(colors)

    for p,path in enumerate(input_paths):
        print(path)

        lengths = list()

        _, file_type = os.path.splitext(path)

        if file_type not in VALID_FILE_TYPES:

            if path not in PRESET_LENGTHS:
                exit("ERROR: '%s' file type or input preset does not match any accepted value, see help for details" % path)
            else:
                lengths = PRESET_LENGTHS[path]

        else:
            lengths = get_lengths_from_fasta(path)

        name = os.path.basename(path).split('.fa')[0]
        names.append(name)

        color_index = p
        if color_indexes is not None:
            color_index = color_indexes[p]

        color = colors[color_index]
        label_colors.append(color)

        plot_ngx(
            name=name,
            color=color,
            fig=fig,
            axes=axes,
            lengths=lengths,
            genome_size=genome_size,
            output_dir=output_dir
        )

        if max(lengths) > max_len:
            max_len = max(lengths) + 2500000

    if legends:
        # Generate custom legend/key lines
        print('making legend')
        custom_lines = list()
        for c in list(set(label_colors)):
            custom_lines.append(Line2D([0], [0], color=c, lw=4))

        pyplot.legend(custom_lines, ["HBCC","NABEC"])



    fig_path = os.path.join(output_dir, "ngx.png")
    pyplot.savefig(fig_path, dpi=200)

    fig_path = os.path.join(output_dir, "ngx.pdf")
    pyplot.savefig(fig_path, dpi=200)

    pyplot.close()


def parse_comma_separated_string(s):
    return s.strip().split(',')


def parse_comma_separated_int_string(s):
    if s is None:
        return []
    else:
        return list(map(int,s.strip().split(',')))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="Comma separated list of: input fasta file paths and/or a preset name: " + ','.join(PRESET_LENGTHS.keys())
    )

    parser.add_argument(
        "-o",
        required=True,
        type=str,
        help="Output directory"
    )

    parser.add_argument(
        "-c","--colors",
        required=False,
        default=None,
        type=str,
        help="Comma separated list of color indexes to use for each item in the plot: "
             "0=red, "
             "1=orange, "
             "2=yellow, "
             "3=light green, "
             "4=green, "
             "5=green-blue, "
             "6=turquoise, "
             "7=blue, "
             "8=indigo, "
             "9=purple, "
             "10=pink, "
    )

    parser.add_argument(
        "-g",
        required=True,
        type=int,
        help="Size of genome"
    )

    parser.add_argument(
        "-t",
        required=False,
        default=1,
        type=int,
        help="Maximum number of threads to use"
    )

    parser.add_argument(
        "-l",
        required=False,
        default=True,
        type=bool,
        help="legend Boolean"
    )



    args = parser.parse_args()

    args.i = parse_comma_separated_string(args.i)

    if args.colors is not None:
        args.colors = parse_comma_separated_int_string(args.colors)

        if len(args.colors) != len(args.i):
            exit("ERROR: color list length doesn't match list of input files")

        # for c in args.colors:
        #     if 0 > c > 10:
        #         exit("ERROR: color argument must be between 0 and 10 (inclusive), see help for details")

    main(input_paths=args.i, output_dir=args.o, genome_size=args.g, color_indexes=args.colors, n_threads=args.t,
         legends=args.l)
