#!/usr/bin/python3.4
import pandas as pd
import numpy as np
import sys, getopt
import os
import re

def main(argv):
    projdir = ''
    outputdir = ''
    pipeline = ''
    usage='> Usage: ctg-sc-count-mterics-concat.py -i PROJECT-OUTDIR -o SUMMARY-OUTDIR -t PIPELINE'

    try:
        opts, args = getopt.getopt(argv,"hi:o:p:",["projdir=", "outdir=", "pipeline="])
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
        elif opt in ("-p", "--pipeline"):
            pipeline = arg

    out1 = outputdir + "/ctg-cellranger-count-summary_metrics.csv"
    out2 = outputdir + "/cellranger-count_summary_metrics_mqc.csv"

    projid=os.path.basename(os.path.normpath(projdir))

    # list all metricsfiles
    samples = os.listdir(projdir + '/count-cr/')

    # get id of metricfile
    sname = samples[0]
    fname = projdir + "/qc/cellranger/" + sname + ".metrics_summary.csv"
    final = pd.read_csv(fname)
    final.index = [projid + "-" + sname]

    # concatenate all sample tables to one
    for sname in samples[1:]:
        fname = projdir + "/qc/cellranger/" + sname + ".metrics_summary.csv"
        data = pd.read_csv(fname)
        data.index = [projid + "-" + sname]
        final = pd.concat([final,data],axis=0)

    # Write csv file
    final.replace(",","",regex=True)
    final.index.name = "Sample"
    final.to_csv(out1,sep=",")

    mqdf = pd.read_csv(out1)
    # Write header
    f = open(out2,'a')
    f.write("# plot_type: 'table'" + "\n")
    f.write("# section_name: '10x Cellranger Metrics'\n")
    f.write("# description: '10x %s count metrics'\n" % pipeline)
    f.write("# pconfig:\n")
    f.write("#     namespace: 'CTG'\n")
    f.write("# headers:\n")
    # iterate over column names
    colidx=1
    for line in mqdf:
        f.write("#     col"+str(colidx)+":\n")
        f.write("#         title: '"+line+"'\n")
        f.write("#         description: '"+line+"'\n")
        #check datatype for formating
        form="na"
        if "," in mqdf[line][1]:
            f.write("#         format: '{:.0f}' \n")
        if "%" in mqdf[line][1]:
            f.write("#         min: 0\n")
            f.write("#         max: 100\n")
            f.write("#         format: '{:.1f}'\n")
            f.write("#         suffix: '%'\n")
        colidx = colidx+1
    f.close()
    # fix colnames
    columns = mqdf.columns
    newcols = ['Sample-ID']
    colidx=2
    for col in range(1,len(columns)):
        newcols.append("col"+str(colidx))
        colidx = colidx+1
    newcols
    mqdf.columns = newcols
    outdf=mqdf.copy(deep=True)
    # Iterate over columns to get format of each
    for rowIndex, row in mqdf.iterrows(): #iterate over rows
        for columnIndex, value in row.items():
            newval=""
            if "%" in value:
                newval=value[:-1]
                outdf.at[rowIndex,columnIndex] = newval
            if "," in value:
                newval=re.sub(",","",str(value))
                outdf.at[rowIndex,columnIndex] = newval
    # Write header with format
    outdf.to_csv(out2,mode="a",sep="\t",index=False)

if __name__ == "__main__":
    main(sys.argv[1:])
