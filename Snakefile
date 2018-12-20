###### Config file of parameters #####
configfile: "config.yaml"


rule all:
    input:
        config["OutputPrefix"] + "FinalPrimerList.tab"

### Prepare database and sorted input files ###

rule MakeBlastDatabase:
    input:
        config["InputFiles"]["GenomeFasta"]
    output:
        config["InputFiles"]["GenomeFasta"] + ".nhr",
        config["InputFiles"]["GenomeFasta"] + ".nin",
        config["InputFiles"]["GenomeFasta"] + ".nsq"
    shell:
        "makeblastdb -nucl -in {input}"

rule FastaIndex:
    input:
        config["InputFiles"]["GenomeFasta"]
    output:
        config["InputFiles"]["GenomeFasta"] + ".fai",
    shell:
        "samtools index {input}"

rule SortTargetSpliceSites:
    input:
        config["InputFiles"]["TargetSpliceSites"]
    output:
        config["OutputPrefix"] + "TargetSpliceSites.sorted.bed"
    params:
        GenomeIndex=config["InputFiles"]["GenomeFasta"] + ".fai",
    shell:
        "bedtools sort -g {params.GenomeIndex} -i {input} > {output}"

if config["InputFiles"]["RepeatRegions"]:
    rule SortRepeatRegions:
        input:
            config["InputFiles"]["RepeatRegions"]
        output:
            config["OutputPrefix"] + "RepeatRegions.sorted.bed"
        params:
            GenomeIndex=config["InputFiles"]["GenomeFasta"] + ".fai",
        shell:
            "bedtools sort -g {params.GenomeIndex} -i {input} > {output}"
else:
    rule MakePseudoRepeatRegions:
        output:
            config["OutputPrefix"] + "RepeatRegions.sorted.bed"
        params:
            GenomeIndex=config["InputFiles"]["GenomeFasta"] + ".fai",
        shell:
            "touch {output}"

if config["InputFiles"]["AnnotatedExons"]:
    rule SortAnnotatedExons:
        input:
            config["InputFiles"]["AnnotatedExons"]
        output:
            config["OutputPrefix"] + "AnnotatedExons.sorted.bed"
        params:
            GenomeIndex=config["InputFiles"]["GenomeFasta"] + ".fai",
        shell:
            "bedtools sort -g {params.GenomeIndex} -i {input} > {output}"
else:
    rule MakePseudoAnnotatedExons:
        output:
            config["OutputPrefix"] + "AnnotatedExons.sorted.bed"
        params:
            GenomeIndex=config["InputFiles"]["GenomeFasta"] + ".fai",
        shell:
            "bedtools makewindows -n 1 -g {params.GenomeIndex} > {output}"

if config["InputFiles"]["RNASeq"]:
    rule SortRNASeq:
        input:
            RNASeqBam = config["InputFiles"]["RNASeq"],
            GenomeIndex = config["InputFiles"]["GenomeFasta"] + ".fai",
        output:
            config["OutputPrefix"] + "RNASeq.sorted.bam"
        params:
            GenomeFasta=config["InputFiles"]["GenomeFasta"]
        shell:
            """
            samtools reheader <(cat <(samtools view -H {input.RNASeqBam} | awk '$0 ~ /^@HD/') <(samtools view -H {input.RNASeqBam} | awk '$0 ~ /^@SQ/ {{split($2,a,":"); print a[2],$0}}' | awk 'FNR == NR {{ lineno[$1] = NR; next}} {{print lineno[$1], $0;}}' {input.GenomeIndex} - | sort -k 1,1n | awk -v OFS='\t' '{{print $3,$4,$5}}' ) <(samtools view -H {input.RNASeqBam} | awk '$0 !~ /^@SQ/ && $0 !~ /^@HD/')) {input.RNASeqBam} | samtools sort > {output}
            """
else:
    rule MakePseudoRNASeq:
        output:
            config["OutputPrefix"] + "RNASeq.sorted.bam"
        params:
            GenomeIndex=config["InputFiles"]["GenomeFasta"] + ".fai",
        shell:
            "bedtools makewindows -n 1 -g {params.GenomeIndex} > {output}"


### Design primers ###

rule CreateWindowsForCandidatePrimersFromTargetList:
    #Write primer design windows to new bedfile. If the region overlaps an exon, constrain the region to the exon edges
    input:
        GenomeIndex = config["InputFiles"]["GenomeFasta"] + ".fai",
        AnnotatedExons = config["OutputPrefix"] + "AnnotatedExons.sorted.bed",
        TargetSpliceSites = config["OutputPrefix"] + "TargetSpliceSites.sorted.bed"
    output:
        config["OutputPrefix"] + "PrimerWindows.bed"
    params:
        WindowSize = config["Parameters"]["WindowSize"],
        PathToHelperScript = "scripts/ConstrainTargetsToOverlapSingleExons_helper1.pl"
    shell:
        """
        cat <(bedtools intersect -sorted -a <(bedtools flank -g {input.GenomeIndex} -l 0 -r {params.WindowSize} -s -i {input.TargetSpliceSites}) -b {input.AnnotatedExons} -g {input.GenomeIndex} -s -wo | perl {params.PathToHelperScript} ) <(bedtools intersect -sorted -a <(bedtools flank -g {input.GenomeIndex} -l 0 -r {params.WindowSize} -s -i {input.TargetSpliceSites}) -b {input.AnnotatedExons} -g {input.GenomeIndex} -s -wa -v) | awk -F '\t' -v OFS='\t' '$6=="+" {{print $1,$2,$3,$4,$5,"-"}} $6=="-" {{print $1,$2,$3,$4,$5,"+"}}' | bedtools sort -i - -faidx {input.GenomeIndex} > {output}
        """

rule CreateAllCandidatePrimerSequencesFromWindows:
    #Get sequence of each 24-32mer in each window
    input:
        config["OutputPrefix"] + "PrimerWindows.bed"
    output:
        config["OutputPrefix"] + "CandidatePrimers.AllKmerToLmers.tab"
    params:
        GenomeFasta = config["InputFiles"]["GenomeFasta"],
        MinPrimerLen = config["Parameters"]["MinPrimerLength"],
        MaxPrimerLen = config["Parameters"]["MaxPrimerLength"],
        PathToHelperScript = "scripts/GetOverlapping_Kmer-Lmers.py"
    shell:
        """
        bedtools getfasta -fi {params.GenomeFasta} -bed {input} -s -tab -name+ | python {params.PathToHelperScript} {params.MinPrimerLen} {params.MaxPrimerLen} > {output}
        """

rule FilterOutPrimersByTm:
    #filter for Tm in x-y window
    input:
        config["OutputPrefix"] + "CandidatePrimers.AllKmerToLmers.tab"
    output:
        config["OutputPrefix"] + "CandidatePrimers.tmfiltered.tab"
    params:
        MinTm = config["Parameters"]["TmFilter"]["MinTm"],
        MaxTm = config["Parameters"]["TmFilter"]["MaxTm"],
    shell:
        """
        paste {input} <(awk '{{print $2}}' {input} | xargs -L1 oligotm -n 1 -mv 75 -dv 3 -d 10 -tp 1 -sc 1) | awk '$3<{params.MaxTm} && $3>{params.MinTm}' > {output}
        """

rule FilterOutPrimersThatOverlapRepetitiveRegions:
    #filter out primers that >50% overlap repeat region
    input:
        GenomeIndex = config["InputFiles"]["GenomeFasta"] + ".fai",
        CandidatePrimers = config["OutputPrefix"] + "CandidatePrimers.tmfiltered.tab",
        RepeatRegions = config["OutputPrefix"] + "RepeatRegions.sorted.bed"
    output:
        config["OutputPrefix"] + "CandidatePrimers.RepeatFiltered.bed"
    params:
        config["Parameters"]["RepeatRegionMaxPercentOverlap"]
    shell:
        """
        awk -v OFS='\t' '{{split($1, namelist, "::"); split(namelist[3], coord, "_"); print coord[1], coord[2], coord[3], $1, ".", coord[4], $3, $2}}' {input.CandidatePrimers} | bedtools sort -i - -faidx {input.GenomeIndex} | bedtools intersect -a - -b {input.RepeatRegions} -sorted -g {input.GenomeIndex} -f {params} -v > {output}
        """

rule BlastCandidatePrimers:
    input:
        CandidatePrimers = config["OutputPrefix"] + "CandidatePrimers.RepeatFiltered.bed",
        fasta = config["InputFiles"]["GenomeFasta"],

        #Pre-made Blast database should exist, therefore, specified as input in this rule
        nhr = config["InputFiles"]["GenomeFasta"] + ".nhr",
        nin = config["InputFiles"]["GenomeFasta"] + ".nin",
        nsq = config["InputFiles"]["GenomeFasta"] + ".nsq"
    output:
        config["OutputPrefix"] + "CandidatePrimers.blastresults.tab"
    shell:
        """
        awk -F'\t' '{{print ">"$4"::"$7"\n"$8}}' {input.CandidatePrimers} | blastn -task blastn-short -db {input.fasta} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send btop nident qlen evalue" > {output}
        """

rule CalculateExpressionOfBlastHits:
    input:
        BlastResults = config["OutputPrefix"] + "CandidatePrimers.blastresults.tab",
        RNASeq = config["OutputPrefix"] + "RNASeq.sorted.bam",
        GenomeIndex=config["InputFiles"]["GenomeFasta"] + ".fai"
    output:
        config["OutputPrefix"] + "CandidatePrimers.blastresults.counted.tab"
    shell:
        """
        bedtools sort -faidx {input.GenomeIndex} -i <(awk -F'\t' -v OFS='\t' '$10>$9 {{print $2, $9,$10, ".", ".", "+", $0}} $9>$10 {{print $2, $10,$9, ".", ".", "-", $0}}'  {input.BlastResults}) | bedtools coverage -split -a - -sorted -g {input.GenomeIndex} -b {input.RNASeq} -mean | awk -F'\t' -v OFS='\t' '{{$1=$2=$3=$4=$5=$6=""; print $0}}' > {output}
        """

rule InSilicoPredictionOfOnAndOffTargetCounts:
    input:
        BlastResultsExpressionCounts = config["OutputPrefix"] + "CandidatePrimers.blastresults.counted.tab",
        RNASeq = config["OutputPrefix"] + "RNASeq.sorted.bam",
        GenomeIndex=config["InputFiles"]["GenomeFasta"] + ".fai"
    output:
        config["OutputPrefix"] + "CandidatePrimers.blastresults.OnTargetRatioAdded.tab"
    params:
        PathToHelperScript = "scripts/EstimateOffTarget_helper.pl"
    shell:
        """
        cat {input.BlastResultsExpressionCounts} | awk -F'\t' -v OFS='\t' '{{split($7,a,"::"); print a[1]":"a[2]":"a[3]":"a[5]":"a[4],$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}}' |  awk -F'\t' -v OFS='\t' '$8==$13' | sort -k1,1 -k14g,14 | perl {params.PathToHelperScript} | awk -F '\t' -v OFS='\t' '{{print $1, ($15*$19)}}' | awk -F'\t' -v OFS='\t' '{{arr[$1]+=$2}} END {{for (i in arr) {{print i,arr[i]}} }}' | awk -F'\t' -v OFS='\t' '{{split($1,a,":"); split(a[3],b,"_"); print b[1],b[2],b[3],$1,".",b[4], $2}}' | bedtools sort -faidx {input.GenomeIndex} -i - | bedtools coverage -split -a - -sorted -g {input.GenomeIndex} -b  {input.RNASeq} -mean
        """

rule ScorePrimersAndPickTheBestForEachTarget:
    input:
        PythonPrimerScorerInput = config["OutputPrefix"] + "CandidatePrimers.blastresults.OnTargetRatioAdded.tab",
    output:
        AllPrimersScore = config["OutputPrefix"] + "CandidatePrimers.scores.tab",
        FinalPrimerList = config["OutputPrefix"] + "FinalPrimerList.tab"
    params:
        MinPrimerLen = config["Parameters"]["MinPrimerLength"],
        MaxPrimerLen = config["Parameters"]["MaxPrimerLength"],
        MinTm = config["Parameters"]["TmFilter"]["MinTm"],
        OptimalTm = config["Parameters"]["TmFilter"]["MinTm"],
        MaxTm = config["Parameters"]["TmFilter"]["MaxTm"],
        Tm_weight = config["Parameters"]["Tm_weight"],
        PercentOnTargetEstimate_weight = config["Parameters"]["PercentOnTargetEstimate_weight"],
        GCcontent_weight = config["Parameters"]["GCcontent_weight"],
    shell:
        """
        cat {input.PythonPrimerScorerInput} | python PickBestScoringPrimers.py {params} | tee {output.AllPrimersScore} | sort -nrk13,13 | sort -u -k14,14 | awk -F '\t' -v OFS='\t' '{{split($4,a,":"); print $0, a[5]}}' > {output.FinalPrimerList}
        """
