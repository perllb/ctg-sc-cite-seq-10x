#!/usr/bin/env nextFlow

// Base params
exp = params.experiment
basedir = params.basedir
metaid = params.metaid

// Output dirs
outdir = params.outdir
fqdir = params.fqdir
qcdir = params.qcdir
countdir = params.countdir
aggdir = params.aggdir
ctgqc = params.ctgqc
metadir = params.metadir

// Input feature reference
features = params.features

// Demux args
b2farg_rna = params.bcl2fastqarg_rna
b2farg_adt = params.bcl2fastqarg_adt
index_rna = params.index_rna
index_adt = params.index_adt
demux = params.demux

// Read and process CTG samplesheet 
sheet = file(params.sheet)

// create new samplesheet in cellranger mkfastq IEM (--samplesheet) format. This will be used only for demultiplexing
newsheet_rna = "$metadir/samplesheet.nf.cite-seq-adt-10x.rna.csv"
newsheet_adt = "$metadir/samplesheet.nf.cite-seq-adt-10x.adt.csv"

println "============================="
println ">>> cite-seq-adt-10x pipeline >>>"
println ""
println "> INPUT: "
println ""
println "> run-meta-id		: $metaid "
println "> basedir		: $basedir "
println "> runfolder		: $exp "
println "> sample-sheet		: $sheet "
println "> feature-ref          : $features "
println ""
println " - demux settings " 
println "> bcl2fastq-arg-rna    : '${b2farg_rna}' "
println "> bcl2fastq-arg-adt    : '${b2farg_adt}' "
println "> demux                : $demux " 
println "> index-rna            : $index_rna "
println "> index-adt            : $index_adt "
println "> metadir              : $metadir "
println ""
println " - output directory structure "
println "> outdir               : $outdir "
println "> fastq                : $fqdir "
println "> qc                   : $qcdir "
println "> count                : $countdir " 
println "> aggregated           : $aggdir "
println "> ctg-qc               : $ctgqc "
println "> metadata             : $metadir "
println ""
println "============================="


// extract RNA samplesheet info
Channel
    .fromPath(sheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Project, row.Sample_Species, row.Sample_Lib, row.Sample_Pair ) }
    .tap{infoall}
    .into { crlib_ch; cragg_ch; fqc_ch }

// RNA Projects
Channel
    .fromPath(sheet)
    .splitCsv(header:true)
    .map { row -> row.Sample_Project }
    .unique()
    .tap{infoProject}
    .set { count_summarize  }


// Channel to start count if demux == 'n'
if ( demux == 'n' ) {
   Channel
	 .fromPath(sheet)
	 .splitCsv(header:true)
    	 .map { row -> tuple( row.Sample_ID, row.Sample_Project, row.Sample_Species, row.Sample_Lib, row.Sample_Pair ) }
    	 .set { count_lib_csv }

   Channel
	 .fromPath(sheet)
    	 .into { count_rna; count_adt }

}

println " > Samples to process: "
println "[Sample_ID,Sample_Project,Sample_Species,Sample_Lib,pair]"
infoall.subscribe { println "Info: $it" }

println " > RNA Projects to process : "
println "[Sample_Project]"
infoProject.subscribe { println "Info Projects: $it" }

// Parse RNA samplesheet 
process parsesheet_rna {

	tag "$metaid"

	input:
	val sheet
	val index_rna

	output:
	val newsheet_rna into demux_sheet_rna

	when:
	demux == 'y'

	"""
mkdir -p $metadir
cat $sheet | grep \',rna,\\|Lane,Sample_ID\' > tmp.sheet.RNA.csv
python $basedir/bin/ctg-parse-samplesheet.10x.py -s tmp.sheet.RNA.csv -o $newsheet_rna -i $index_rna
rm tmp.sheet.RNA.csv
	"""
}

// Parse ADT samplesheet 
process parsesheet_adt {

	tag "$metaid"

	input:
	val sheet
	val index_adt

	output:
	val newsheet_adt into demux_sheet_adt

	when:
	demux == 'y'

	"""
mkdir -p $metadir
cat $sheet | grep \',adt\\|Lane,Sample_ID\' > tmp.sheet.ADT.csv
python $basedir/bin/ctg-parse-samplesheet.10x.py -s tmp.sheet.ADT.csv -o $newsheet_adt -i $index_rna
rm tmp.sheet.ADT.csv
	"""
}

// aggregation
process gen_libraries_csv {

    	tag "${sid}_${projid}"

	input:
	val sheet
	set sid, projid, ref, lib, pair from crlib_ch

 	output:
	set sid, projid, ref, lib, pair into count_lib_csv

	when:
	lib == 'rna'

	"""
mkdir -p $metadir

libcsv=$metadir/${projid}_${sid}_libraries.csv

# Print header
echo 'fastqs,sample,library_type' > \$libcsv
# Print RNA entry
echo '${fqdir}/rna/${projid},$sid,Gene Expression' >> \$libcsv
# Get paired ADT sample
adtid=\$(grep ',adt,$pair' $sheet | cut -f2 -d ',')
echo "${fqdir}/adt/${projid},\$adtid,Antibody Capture" >> \$libcsv

        """
}

// Run RNA mkFastq
process mkfastq_rna {

	tag "${metaid}-rna"

	input:
	val rnasheet from demux_sheet_rna
	
	output:
	val "count" into count_rna
	val "gorna" into fastqc_go_rna

	when:
	demux == 'y'

	"""
if [ '$index_rna' == 'dual' ]; then
   indexarg='--filter-dual-index'
else
   indexarg='--filter-single-index'
fi

cellranger mkfastq \\
	   --id=${metaid}_adt \\
	   --run=$exp \\
	   --samplesheet=$rnasheet \\
	   --jobmode=local \\
	   --localmem=150 \\
	   --output-dir ${fqdir}/rna \\
	   \${indexarg} \\
	   $b2farg_rna
"""
}

// Run ADT mkFastq
process mkfastq_adt {

	tag "${metaid}-adt"

	input:
	val adtsheet from demux_sheet_adt

	output:
	val "count" into count_adt
	val "goadt" into fastqc_go_adt

	when:
	demux == 'y'

	"""
if [ '$index_rna' == 'dual' ]; then
   indexarg='--filter-dual-index'
else
   indexarg='--filter-single-index'
fi

cellranger mkfastq \\
	   --id=${metaid}_adt \\
	   --run=$exp \\
	   --samplesheet=$adtsheet \\
	   --jobmode=local \\
	   --localmem=150 \\
	   --output-dir ${fqdir}/adt \\
	   \${indexarg} \\
	   $b2farg_adt
"""
}

// count RNA + ADT
process count {

	tag "${sid}-${projid}"
	publishDir "${countdir}/", mode: "copy", overwrite: true

	input: 
	val rnaready from count_rna
	val adtready from count_adt
	set sid, projid, ref, lib, pair from count_lib_csv

	output:
	file "${sid}/outs/" into samplename
	val "${qcdir}/cellranger/${sid}.metrics_summary.csv" into count_metrics
	val "${aggdir}/${sid}.molecule_info.h5" into count_agg

	when:
	lib == 'rna'

	"""
if [ $ref == "Human" ] || [ $ref == "human" ]
then
	genome="/projects/fs1/shared/references/hg38/cellranger/refdata-gex-GRCh38-2020-A"
elif [ $ref == "mouse" ] || [ $ref == "Mouse" ]
then
	genome="/projects/fs1/shared/references/mm10/cellranger/refdata-gex-mm10-2020-A"
elif [ $ref == "custom"  ] || [ $ref == "Custom" ] 
then
	genome=${params.custom_genome}
else
	echo ">SPECIES NOT RECOGNIZED!"
	genome="ERR"
fi

mkdir -p ${countdir}

libcsv=$metadir/${projid}_${sid}_libraries.csv

cellranger count \\
	--id=$sid \\
	--libraries=\$libcsv \\
	--transcriptome=\$genome \\
	--feature-ref=$features \\
	--localmem=150 \\
	--jobmode=local \\
	--localcores=${task.cpus} 

## Copy h5 file for aggregation
mkdir -p $aggdir
cp ${sid}/outs/molecule_info.h5 ${aggdir}/${sid}.molecule_info.h5

## Copy metrics file for qc
# Remove if it exists
if [ -f ${qcdir}/cellranger/${sid}.metrics_summary.csv ]; then
	rm -r ${qcdir}/cellranger/${sid}.metrics_summary.csv
fi
mkdir -p ${qcdir}
mkdir -p ${qcdir}/cellranger/
cp ${sid}/outs/metrics_summary.csv ${qcdir}/cellranger/${sid}.metrics_summary.csv

## Copy to delivery folder 
mkdir -p ${outdir}/summaries
mkdir -p ${outdir}/summaries/cloupe
mkdir -p ${outdir}/summaries/web-summaries
cp ${sid}/outs/web_summary.html ${outdir}/summaries/web-summaries/${sid}.web_summary.html
cp ${sid}/outs/cloupe.cloupe ${outdir}/summaries/cloupe/${sid}_cloupe.cloupe

## Copy to CTG-QC dir 
mkdir -p ${ctgqc}
mkdir -p ${ctgqc}/web-summaries
cp ${sid}/outs/web_summary.html ${ctgqc}/web-summaries/${sid}.web_summary.html

	"""

}

process summarize_count {

	tag "${projid}"

	input:
	val metrics from count_metrics.collect()

	output:
	val "y" into mqc_count 	
	val "x" into run_summarize

	"""
cd $outdir

mkdir -p ${qcdir}
mkdir -p ${qcdir}/cellranger

# RNA summaries
python $basedir/bin/ctg-sc-rna-count-metrics-concat.py -i ${outdir} -o ${qcdir}/cellranger
# ADT summaries
python $basedir/bin/ctg-sc-adt-count-metrics-concat.py -i ${outdir} -o ${qcdir}/cellranger
	"""
}

// aggregation
process gen_aggCSV {

	tag "${sid}_${projid}"

	input:
	set sid, projid, ref, lib, pair from cragg_ch

	output:
	val projid into craggregate

	when:
	lib == 'rna'

	"""
mkdir -p ${aggdir}

aggcsv=${aggdir}/${projid}_libraries.csv

if [ -f \${aggcsv} ]
then
	if grep -q $sid \$aggcsv
	then
		echo ""
	else
		echo "${sid},${aggdir}/${sid}.molecule_info.h5" >> \$aggcsv
	fi
else
	echo "sample_id,molecule_h5" > \$aggcsv
	echo "${sid},${aggdir}/${sid}.molecule_info.h5" >> \$aggcsv
fi


	"""
}

process aggregate {

	publishDir "${outdir}/aggregate/", mode: 'move', overwrite: true
	tag "$projid"
  
	input:
	val projid from craggregate.unique()
	val moleculeinfo from count_agg.collect()

	output:
	file "${projid}_agg/outs" into doneagg
	val projid into md5_agg_go

	"""
cellranger aggr \
   --id=${projid}_agg \
   --csv=${aggdir}/${projid}_libraries.csv \
   --normalize=mapped

## Copy to delivery folder 
cp ${projid}_agg/outs/web_summary.html ${outdir}/summaries/web-summaries/${projid}_agg.web_summary.html
cp ${projid}_agg/outs/count/cloupe.cloupe ${outdir}/summaries/cloupe/${projid}_agg_cloupe.cloupe

## Copy to CTG QC dir 
cp ${outdir}/summaries/web-summaries/${projid}_agg.web_summary.html ${ctgqc}/web-summaries/
cp ${outdir}/summaries/cloupe/${projid}_agg_cloupe.cloupe ${ctgqc}/web-summaries/

## Remove the molecule_info.h5 files that are stored in the aggregate folder (the original files are still in count-cr/../outs 
rm ${aggdir}/*h5

	"""
}

process fastqc {

	tag "${sid}-${projid}"

	input:
	val gorna from fastqc_go_rna
	val goadt from fastqc_go_adt
	set sid, projid, ref, lib, pair from fqc_ch	
		
	output:
	val projid into mqc_cha
	val projid into md5_fastqc_go

	"""
mkdir -p ${qcdir}
mkdir -p ${qcdir}/fastqc

for file in ${fqdir}/${lib}/${projid}/${sid}*fastq.gz
	do fastqc -t ${task.cpus} \$file --outdir=${qcdir}/fastqc
done
	"""
}

process multiqc_count_run {

	tag "${metaid}"

	input:
	val x from run_summarize.collect()
	val projid from mqc_cha.collect()
		
	output:
	val "x" into summarized

	"""
# make copy of adt Stats.json
cp ${fqdir}/adt/Stats/Stats.json ${fqdir}/adt/Stats/Stats.adt.json

cd ${outdir}
mkdir -p ${qcdir}/multiqc
multiqc -f ${fqdir} ${qcdir}/fastqc/ ${qcdir}/cellranger/ ${fqdir}/adt/Stats/Stats.adt.json ${fqdir}/rna/Stats/Stats.json --outdir ${qcdir}/multiqc -n ${metaid}_sc-aft-rna-10x_summary_multiqc_report.html

cp -r ${qcdir} ${ctgqc}

	"""

}

// Final process, when all is done: md5 recursively from output root folder
process md5sum {

	input:
	val v from summarized.collect()
	val projid from md5_fastqc_go.unique()
	val x from md5_agg_go.collect()
	
	"""
# Remove Undetermined files!
rm ${fqdir}/adt/Undetermined*
rm ${fqdir}/rna/Undetermined*

cd ${outdir}
find -type f -exec md5sum '{}' \\; > ctg-md5.${projid}.txt
        """ 

}
