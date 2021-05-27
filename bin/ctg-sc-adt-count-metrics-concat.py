#!/usr/bin/python3.4

import pandas as pd
import numpy as np
import sys, getopt
import os

def main(argv):
    projdir = ''
    outputdir = ''

    usage='> Usage: ctg-sc-count-mterics-concat.py -i PROJECT-OUTDIR -o SUMMARY-OUTDIR'

    try:
        opts, args = getopt.getopt(argv,"hi:o:",["projdir=", "outdir="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    if len(sys.argv) <= 2:
        print("> Error: No project dir / output dir entered:")
        print(usage)
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-i", "--projdir"):
            projdir = arg
        elif opt in ("-o", "--outdir"):
            outputdir = arg

    out1 = outputdir + "/ctg-cellranger-adt-count-summary_metrics.csv"
    out2 = outputdir + "/cellranger-adt-count_summary_metrics_mqc.csv"

    # list all metricsfiles
    samples = os.listdir(projdir + '/count/')

    # get id of metricfile
    sname = samples[0]
    fname = projdir + "/qc/cellranger/" + sname + ".metrics_summary.csv"
    final = pd.read_csv(fname)
    final.index = [sname]

    # concatenate all sample tables to one
    for sname in samples[1:]:
        fname = projdir + "/qc/cellranger/" + sname + ".metrics_summary.csv"
        data = pd.read_csv(fname)
        data.index = [sname]
        final = pd.concat([final,data],axis=0)
    
    # Write csv file        
    cols = final.columns.tolist()
    cols = list( cols[i] for i in [19,20,26,21,32,22,23,24,25,27,28,29,30,31] )
    final = final[cols]
    final.replace(",","",regex=True)
    final.index.name = "Sample"
    final.to_csv(out1,sep=",")

    # Parse csv file to mqc
    # parse % input
    def p2f(x):
        return float(x.strip('%'))
    # parse integers with comma
    def s2i(x):
        return int(x.replace(",",""))

    mqdf = pd.read_csv(out1,
                       converters={'Antibody: Number of Reads':s2i,
                                   'Antibody: Mean Reads per Cell':s2i,
                                   'Antibody: Fraction Antibody Reads':p2f,
                                   'Antibody: Valid Barcodes':p2f,
                                   'Antibody: Median UMIs per Cell (summed over all recognized antibody barcodes)':s2i,
                                   'Antibody: Sequencing Saturation':p2f,
                                   'Antibody: Q30 Bases in Barcode':p2f,
                                   'Antibody: Q30 Bases in Antibody Read':p2f,
                                   'Antibody: Q30 Bases in UMI':p2f,
                                   'Antibody: Fraction Antibody Reads Usable':p2f,
                                   'Antibody: Antibody Reads Usable per Cell':s2i,
                                   'Antibody: Fraction Antibody Reads in Aggregate Barcodes':p2f,
                                   'Antibody: Fraction Unrecognized Antibody':p2f,
                                   'Antibody: Antibody Reads in Cells':p2f
                               })
    
    orig_cols = mqdf.columns
    mqdf.columns = ['SampleID','col2','col3','col4','col5','col6','col7','col8','col9','col10','col11','col12','col13','col14','col15']
    
    f = open(out2,'a')
    f.write("# plot_type: 'table'" + "\n")
    f.write("# section_name: 'Cellranger Metrics (Antibodies)'\n")
    f.write("# description: 'Cellranger 10x-sc-adt-rna count metrics for antibody libraries'\n")
    f.write("# pconfig:\n")
    f.write("#     namespace: 'CTG'\n") 
    f.write("# headers:\n")
    f.write("#     col1:\n")
    f.write("#         title: 'Sample'\n")
    f.write("#         description: 'Sample ID'\n")
    f.write("#     col2:\n")
    f.write("#         title: 'Antibody Reads: Total Number'\n")
    f.write("#         description: 'Number of Antibody reads'\n")
    f.write("#         format: '{:,.0f}'\n")
    f.write("#     col3:\n")
    f.write("#         title: 'Average Antibody Reads per Cell'\n")
    f.write("#         description: 'Mean Antibody Reads per Cell'\n")
    f.write("#         format: '{:,.0f}'\n")
    f.write("#     col4:\n")
    f.write("#         title: 'Fraction antibody reads'\n")
    f.write("#         description: 'Fraction antibody reads'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#     col5:\n")
    f.write("#         title: 'Valid Barcodes'\n")
    f.write("#         description: 'Valid barcodes'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#     col6:\n")
    f.write("#         title: 'Antibody UMI: Median UMIs per Cell'\n")
    f.write("#         description: 'Median UMI per Cell (summed over all recognized antibody barcodes)'\n")
    f.write("#         format: '{:,.0f}'\n")
    f.write("#     col7:\n")
    f.write("#         title: 'Sequencing Saturation'\n")
    f.write("#         description: 'Sequencing Saturation of Antibody library'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#     col8:\n")
    f.write("#         title: 'Q30 bases in Barcodes'\n")
    f.write("#         description: 'Fraction of bases with >Q30 in barcodes'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#     col9:\n")
    f.write("#         title: 'Q30 bases in Antibody Reads'\n")
    f.write("#         description: 'Fraction of bases with >Q30 in Antibody Reads'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#     col10:\n")
    f.write("#         title: 'Q30 bases in UMI'\n")
    f.write("#         description: 'Fraction of bases with >Q30 in UMI'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#     col11:\n")
    f.write("#         title: 'Fraction Antibody Reads Usable'\n")
    f.write("#         description: 'Fraction Antibody Reads Usable'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col12:\n")
    f.write("#         title: 'Usable Antibody Reads Per Cell'\n")
    f.write("#         description: 'Antibody Reads Usable Per Cell'\n")
    f.write("#         format: '{:,.0f}'\n")
    f.write("#     col13:\n")
    f.write("#         title: 'Fraction Antibody Reads in Aggregate Barcodes'\n")
    f.write("#         description: 'Fraction Antibody Reads in Aggregate Barcodes'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col14:\n")
    f.write("#         title: 'Fraction Unrecognized Antibody'\n")
    f.write("#         description: 'Fraction Unrecognized Antibody'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    f.write("#     col15:\n")
    f.write("#         title: 'Antibody Reads in Cells'\n")
    f.write("#         description: 'Antibody Reads in Cells'\n")
    f.write("#         format: '{:.1f}'\n")
    f.write("#         suffix: '%'\n")
    f.write("#         min: 0 \n")
    f.write("#         max: 100 \n")
    mqdf.to_csv(f,sep="\t",index=False)
    f.close()
    
if __name__ == "__main__":
    main(sys.argv[1:])
