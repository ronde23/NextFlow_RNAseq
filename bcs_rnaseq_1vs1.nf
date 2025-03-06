#!/usr/bin/env nextflow

import java.nio.file.Path
import java.nio.file.Paths

// You need to change the input and output filenames of the float2intNF.R

// NFCoreRNASEQ % docker buildx build .
// Run: nextflow bcs_rnaseq.nf -params-file bcs_rnaseq.yaml
// Or Run: $ nextflow bcs_rnaseq.nf --basename 'treatment_vs_control' --outdir 'Output' --replicates yes

nextflow.enable.dsl=2

params.str = 'BCS RNA-Seq Pipeline...\nDeveloped by Ronika De.'
println(params.str)

// Inputs
params.yamlfile = "bcs_rnaseq.yaml" // input yaml file
params.organism = 'Mouse' // input organism 
params.sample_input = 'bcs_samplesheet.csv' // input samplesheet file
params.refdir = "/database/nextflow/pipelines/downstream_expression/"
params.basename = 'comparison_vs_control' // comparison groups for e.g., treatment_vs_baseline
params.folder = '' // The upstream nf-core/rnaesq output folder. By default it is 'star_salmon/'. If you change the default setting in the upstream nf-core/rnaesq, please change this input as well.  // Not for yaml file
params.infile = ''; // Gene count file: 'salmon.merged.gene_counts.tsv'
params.tpmfile = ''; // Normalized count (TPM) file 'salmon.merged.gene_tpm.tsv'

// If the raw and normalized gene count files are not provided, it will use the default output files from Salmon.
if ((params.infile != '') && (params.tpmfile != '')) {
	rawf = params.infile
	normf = params.tpmfile
}

if ((params.infile == '') && (params.tpmfile == '')) {
	//Path rinput = Paths.get(params.folder+'/'+'salmon.merged.gene_counts.tsv'); rawf = rinput.getFileName().toString();
	rawf = params.folder+'/'+'salmon.merged.gene_counts.tsv'
	

	//Path ninput = Paths.get(params.folder+'/'+'salmon.merged.gene_tpm.tsv'); normf = ninput.getFileName().toString();
	normf = params.folder+'/'+'salmon.merged.gene_tpm.tsv'
	
}
println(rawf); println(normf); 
//println(params.infile); println(params.tpmfile); 

// Gene count files with rounded (integer) gene count values
params.matrix_feat = 'salmon.merged.gene_roundReadCounts.tsv' // Not for yaml file
params.pca_input = 'salmon.gene_roundReadCounts.tsv' // Not for yaml file
params.tpmout = 'salmon.gene_tpm.tsv' // Not for yaml file

params.outdir = 'Results/' // Output directory
params.replicates = 'yes' // "yes" if replicates are present otherwise "no"

params.foldchange = 0.58 // Expression fold-change threshold to be used
params.significance = 0.05 // Significance (P-value or FDR) threshold to be used

// PCA parameters
params.pca_groups = '' // User chosen groups for performing PCA and Differential gene expression analysis. Recommended: "treatment" & "baseline". You can use other group names as well.
params.pca_topGenes = 500
params.pca_outfname = 'PCA_RNAseq'
params.width = 650
params.height = 600
params.pca_title = 'PCA'
params.scale = 1 // Ranges from 0-1

// GMT files to use for GSEA analysis
params.p1gmtfile = '' // '/data/BCS/Databank/GMT_files/m2.all.v2023.2.Mm.symbols.gmt'
params.p2gmtfile = '' // '/data/BCS/Databank/GMT_files/m5.all.v2023.2.Mm.symbols.gmt'
params.p3gmtfile = '' // '/data/BCS/Databank/GMT_files/mh.all.v2023.2.Mm.symbols.gmt'

if ((params.p1gmtfile != '') || (params.p2gmtfile != '') || (params.p3gmtfile != '')) {
	p1gmtf = params.p1gmtfile
	p2gmtf = params.p2gmtfile
	p3gmtf = params.p3gmtfile
}

params.plotTopX = 50
params.gsea_output = "GSEA_Output" // GSEA output folder name


/* Generting inputs for PCA and Downstream DESeq2 Analysis */
process data_floor {
	publishDir path: "${params.outdir}", mode: 'copy'

	input:
	path rawfile
	path normfile
	
	output:
	path "${params.matrix_feat}", emit: feat_matrix
	path "${params.pca_input}", emit: pca_matrix
	path "${params.tpmout}", emit: tpm_matrix

	script:

	"""
	#!/bin/bash
	echo "${params.sample_input}"
	module load R/4.4.0
	# Rounding the raw and normalized gene counts to the nearest integer
	${params.refdir}float2intNF.R --r ${rawfile} --n ${normfile}
	# Copying the input files to the output directory
	cp ${rawfile} ../../../${params.outdir}
	cp ${normfile} ../../../${params.outdir}
	# Copying the output files to the current working directory
	cp ${params.matrix_feat} ../../../
	cp ${params.pca_input} ../../../
	cp ${params.tpmout} ../../../
	
	"""
}

/* Performing PCA and DEG analysis. Generate PCA plot, volcano plot and heatmap. */
process pca_deseq2 {
	publishDir path: "${params.outdir}", mode: 'copy'

	input:
	path params.pca_input
	path params.matrix_feat
	path normalized_genecount_file
	
	output:
	path "${params.pca_outfname}.pdf", emit: pcafl
	path "*.dsq", emit: dsq

	shell:

	if (params.replicates == 'yes') { 
		rlangs = 'R/4.4.0'
	}
	else { 
		rlangs = 'R/3.4.3' 
	}
	//println(rlangs)
	
	def grpfile = new File(params.pca_input).text.readLines(); //println(grpfile);
	grpheader = grpfile[0]
	grpheader_list = grpheader.split('\t')
	grpheader_list = grpheader_list; println(grpheader_list);

	// Getting the "treatment" and "baseline" groups for performing differential gene expression analysis
	def samplefile = new File(params.sample_input).text.readLines(); //println(samplefile);
	sampleheader = samplefile[0]
	sampleheader_list = sampleheader.split(','); //println(sampleheader_list);
	treatment_comparison_list = ''; baseline_comparison_list = ''; non_comparison_list = '';
	for (int x = 1; x < samplefile.size(); x++) {
		//print(samplefile[x]) 
		split_sample = samplefile[x].split(',')
		sample = split_sample[0]
		grp = split_sample[1]
		comparison = split_sample[2]; //println(split_sample);
		
		if ((sample in grpheader_list) && ((comparison == 'treatment') || (comparison == 'Treatment') || (comparison == 'TREATMENT'))) {
			for (j=0; j<grpheader_list.size(); j++) {
				if (grpheader_list[j] == sample) {
					r = j+1
					treatment_comparison_list += r.toString()+','
				}
			}
		}
		if ((sample in grpheader_list) && ((comparison == 'baseline') || (comparison == 'Baseline') || (comparison == 'BASELINE'))) {
			for (j=0; j<grpheader_list.size(); j++) {
				if (grpheader_list[j] == sample) {
					r = j+1
					baseline_comparison_list += r.toString()+','
				}
			}
		}
		if ((sample in grpheader_list) && (comparison == '')) {
			for (j=0; j<grpheader_list.size(); j++) {
				if (grpheader_list[j] == sample) {
					r = j+1
					non_comparison_list += r.toString()+','
				}
			}		
		}
	}
	
	// Extracting sample groups information for performing PCA analysis
	group_list = ''
	for (int g = 1; g < grpheader_list.size(); g++) {
		for (int k = 1; k < samplefile.size(); k++) {
			splith = samplefile[k].split(','); //println(splith);
			smpl = splith[0]; //println(grouping); println(grpheader_list[g]);
			grouping = splith[1]
			if (smpl == grpheader_list[g]) {
				group_list += grouping+','
			}
		}
	}
	
	
	group_list = group_list.substring(0, group_list.length()-1)
	baseline_comparison_list = baseline_comparison_list.substring(0, baseline_comparison_list.length()-1)
	treatment_comparison_list = treatment_comparison_list.substring(0, treatment_comparison_list.length()-1)
	
	if (non_comparison_list.length() > 0) {
		non_comparison_list = non_comparison_list.substring(0, non_comparison_list.length()-1)
	}
	
	println(treatment_comparison_list); println(baseline_comparison_list); println(group_list);
	
	// If comparison groups are not provided by user, use default groups. This part is not in use.
	if (params.pca_groups != '') { pcagroups = params.pca_groups }
	
	if (params.pca_groups == '') {
		def filecont = new File(params.pca_input).text.readLines(); //println(filecont);
		header = filecont[0]
		header_list = header.split('\t'); //println(header_list);
		header_list = header_list - header_list[0]
		
		pcagroups = ''
		
		for (int i = 0; i < header_list.size(); i++) {
			split_header = header_list[i].split('_')
			samplename = split_header[0]; //println(samplename);
			if (i < header_list.size()-1) {
				pcagroups += samplename+','
				}
			else {
				pcagroups += samplename
			}
		}
	//println(header_list)
	}
	
	//runPCA -m../../../${params.pca_input} -g$pcagroups -v${params.pca_topGenes} -o${params.pca_outfname} -x${params.width} -y${params.height} -t${params.pca_title} -s${params.scale}
	
	"""
	
	#!/bin/bash 
	module purge
	module load $rlangs
	module unload perl/5.36.0
	module load perl/5.20.1
	# Generating matrix file to be used for PCA analysis if you want to use only the comparison groups in the samplesheet.csv file
	#cut -f 1,${treatment_comparison_list},${baseline_comparison_list} ../../../${params.pca_input} > ${params.basename}_matrix.tsv
	# Generating matrix file to be used for PCA analysis if you want to use all the groups in the samplesheet.csv file
	if [ -z "${non_comparison_list}" ]; then
	  cut -f 1,${treatment_comparison_list},${baseline_comparison_list} ../../../${params.pca_input} > ${params.basename}_matrix.tsv
	else
	  cut -f 1,${treatment_comparison_list},${baseline_comparison_list},${non_comparison_list} ../../../${params.pca_input} > ${params.basename}_matrix.tsv
	fi
	
	# Performing differential gene expression analysis by calling the script "runDESeq2_plots"
	${params.refdir}runDESeq2_plots -m../../../${params.pca_input} -t../../../${normalized_genecount_file}  -f${baseline_comparison_list}:${treatment_comparison_list} -o${params.basename}
	# Performing PCA analysis by calling the script "runPCA"
	${params.refdir}runPCA -m${params.basename}_matrix.tsv -g$group_list -v${params.pca_topGenes} -o${params.pca_outfname} -x${params.width} -y${params.height} -t${params.pca_title} -s${params.scale}
	# Copying the output of PCA and differential gene expression analysis to the output directory
	cp -r ${params.pca_outfname}* ../../../${params.outdir}
	cp -r *.dsq ../../../${params.outdir}
	cp -r *.dsq ../../../
	cp -r *.pdf ../../../${params.outdir}
	cp -r *.png ../../../${params.outdir}

	"""
}

/* Filtering the list of DEGs based on expression fold-change and significance threshold  */
process volcano_run {
	publishDir path: "${params.outdir}", mode: 'copy'

	input:
	path dsqf
	float fc = 0.58
	float fdr = 0.05

	output:
	path("*_DFs.txt"), emit: dfs

	script:
	// Selecting "p-value" or "FDR" based on presence or absence of replicates
	if (params.replicates == 'yes') { 
	sig = 'c'
	}
	else { 
	sig = 'p' 
	}
	
	"""

	#!/bin/bash
	module load perl/5.20.1
	# Filtering the genes based on expression fold-change and significance
	${params.refdir}getVolcanoData -d${dsqf} -f${params.foldchange} -$sig${params.significance}
	# Copying the filtered list of DEGs to the output directory
	cp -r *_DFs.txt ../../../${params.outdir}
	# Copying the CCHMC logo, yaml file and RMarkdown script to the output directory
	cp ${params.refdir}/childrens_logo.png ../../../${params.outdir}
	cp -r ../../../${params.yamlfile} ../../../${params.outdir}
	cp ${params.refdir}/RNAseqNF_Markdown.Rmd ../../../${params.outdir}
	
	"""

}

/* Perform GSEA analysis */
/* output:
	path "${params.gsea_output}_p1_${dsqf}", emit: p1
	path "${params.gsea_output}_p2_${dsqf}", emit: p2
	path "${params.gsea_output}_p3_${dsqf}", emit: p3
*/
process run_gsea {
	publishDir path: "${params.outdir}", mode: 'copy'

	input:
	path dsqf
	path p1gmtfl
	path p2gmtfl
	path p3gmtfl
	
	output:
	path "*"

	script:
	
	"""

	#!/bin/bash
	# Generating the rnk file
	cut -f1,3 ${dsqf} > "tmp.rnk"
	tail -n +2 "tmp.rnk" > ${dsqf}.rnk
	rm tmp.rnk
	
	module load jdk/1.8.0_321
	# If File path is not empty (i,e., it contains a value) and the file exists 
	if [ -n "$p1gmtfl" ] && [ -f "$p1gmtfl" ]; then
		# Performing GSEA analysis using hallmark, gene ontology and curated genesets
		${params.refdir}GSEA_Linux_4.3.3/gsea-cli.sh GSEAPreranked -gmx ${p1gmtfl} -norm meandiv -nperm 10000 -rnk ${dsqf}.rnk -scoring_scheme weighted -rpt_label ${dsqf} -create_svgs false -make_sets true -plot_top_x ${params.plotTopX} -rnd_seed timestamp -set_max 500 -set_min 5 -zip_report false -out ${params.gsea_output}_p1_${dsqf.baseName.replaceAll(/\\.dsq$/, '')}
		# Copying the GSEA results to output directory
		cp -r ${params.gsea_output}_p1_${dsqf.baseName.replaceAll(/\\.dsq$/, '')} ../../../${params.outdir}
	fi
	
	if [ -n "$p2gmtfl" ] && [ -f "$p2gmtfl" ]; then
		${params.refdir}/GSEA_Linux_4.3.3/gsea-cli.sh GSEAPreranked -gmx ${p2gmtfl} -norm meandiv -nperm 10000 -rnk ${dsqf}.rnk -scoring_scheme weighted -rpt_label ${dsqf} -create_svgs false -make_sets true -plot_top_x ${params.plotTopX} -rnd_seed timestamp -set_max 500 -set_min 5 -zip_report false -out ${params.gsea_output}_p2_${dsqf.baseName.replaceAll(/\\.dsq$/, '')}
		cp -r ${params.gsea_output}_p2_${dsqf.baseName.replaceAll(/\\.dsq$/, '')} ../../../${params.outdir}
	fi
	
	if [ -n "$p3gmtfl" ] && [ -f "$p3gmtfl" ]; then
		${params.refdir}/GSEA_Linux_4.3.3/gsea-cli.sh GSEAPreranked -gmx ${p3gmtfl} -norm meandiv -nperm 10000 -rnk ${dsqf}.rnk -scoring_scheme weighted -rpt_label ${dsqf} -create_svgs false -make_sets true -plot_top_x ${params.plotTopX} -rnd_seed timestamp -set_max 500 -set_min 5 -zip_report false -out ${params.gsea_output}_p3_${dsqf.baseName.replaceAll(/\\.dsq$/, '')}
		cp -r ${params.gsea_output}_p3_${dsqf.baseName.replaceAll(/\\.dsq$/, '')} ../../../${params.outdir}
	fi
	
	# Removing the raw and normalized gene count files, and other undesired files from the working directory and output directory
	rm -r ../../../${params.basename}.dsq
	rm -r ../../../${params.basename}.pdf
	rm -r ../../../${params.basename}.png
	rm -f ../../../*roundReadCounts.tsv
	rm -f ../../../${params.tpmout}
	rmdir ../../../${params.outdir}/*
	rm -f ../../../${params.outdir}/${params.tpmout}
	rm -f ../../../${params.outdir}/Rplots.pdf
	rm -f ../../../${params.outdir}/*roundReadCounts.tsv
	
	module load perl/5.36.0
	module load R/4.4.0
	module load pandoc
	cd ../../../${params.outdir}
	
	# Executing the Rmarkdown script to generate summary of the downstream RNA-Seq analysis
	Rscript -e "rmarkdown::render('RNAseqNF_Markdown.Rmd', 'html_document')"

	"""
}

workflow {
	def raw_genecount_matrix = Channel.fromPath(rawf)
	def normalized_genecount_matrix = Channel.fromPath(normf)
	def p1gmt = Channel.fromPath(p1gmtf)
	def p2gmt = Channel.fromPath(p2gmtf)
	def p3gmt = Channel.fromPath(p3gmtf)
	data_floor(raw_genecount_matrix, normalized_genecount_matrix)
	pca_deseq2(data_floor.out.pca_matrix, data_floor.out.feat_matrix, data_floor.out.tpm_matrix)
	volcano_run(pca_deseq2.out.dsq)
	run_gsea(pca_deseq2.out.dsq, p1gmt, p2gmt, p3gmt)
}




	