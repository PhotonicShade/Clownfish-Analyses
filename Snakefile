
data = '/data/clownfish/'
gaboriau = ['SRR295227{}'.format(x) for x in range(58, 77)]
anothernine = ['SRR84426{}'.format(x) for x in range(10, 21)]
srrs = gaboriau + anothernine + ['SRR26235419', 'SRR6685843']

left = ['SRR84426{}'.format(x) for x in range(10, 14)]
right = ['SRR84426{}'.format(x) for x in range(15, 21)]
srrss = gaboriau + left + right + ['SRR26235419', 'SRR6685843'] # omit SRR8442614

trim = 'trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar'
libs = ['TruSeq3-PE', 'TruSeq3-SE']

reference = data + 'A_frenatus/afrenatus.transcripts.uniprot.fa'
rdict = data + 'A_frenatus/afrenatus.transcripts.uniprot.dict' # kludge: could use lambda instead

#----------------------------------------------------------------------

rule master :
    input :
        # table for each library
        expand(data + 'SRRs.{lib}.pruned.csv', lib = libs),

        # read qualities
        expand(data + '{s}_fastqc.html', s = srrs),
        expand(data + '{s}.trim:{lib}_fastqc.html', s = srrs, lib = libs),

# call variants, select and filter
#----------------------------------------------------------------------

# prune down to universal rows (no '-'s) with at least one variation
rule prune :
    input : '{path}/SRRs.{lib}.csv'
    output : '{path}/SRRs.{lib}.pruned.csv'
    shell : 'grep -v "-" {input} | python3 scripts/prune.py - > {output}'

# build table with SNPs on the rows and SRRs on the columns for a given library
rule tabulate :
    input :
        expand(data + '{s}.trim:{{lib}}.bwa.filt.mark.rgs.calls.snps.filt.txt', s = srrss)

    output : '{path}/SRRs.{lib}.csv'
    shell : 'python3 scripts/tabulate.py {input} > {output}'

# just keep the uniqe SNP positions
rule snps :
    input : '{path}/SRR{num}.trim:{lib}.bwa.filt.mark.rgs.calls.snps.filt.vcf.gz'
    output : '{path}/SRR{num,[0-9]+}.trim:{lib}.bwa.filt.mark.rgs.calls.snps.filt.txt'
    shell : '''

  zcat {input} \
    | awk '/^#CHROM/ {{s=1}} s {{if($7=="PASS" && length($5)==1) {{print $1" "$2" "$5}}}}' \
      > {output} '''

# filter the SNPs based on various criteria
rule filtration :
    input :
        ref = reference,
        vcf = '{path}/SRR{num}.trim:{lib}.bwa.filt.mark.rgs.calls.snps.vcf.gz'

    output : '{path}/SRR{num,[0-9]+}.trim:{lib}.bwa.filt.mark.rgs.calls.snps.filt.vcf.gz'
    log : '{path}/SRR{num}.trim:{lib}.bwa.filt.mark.rgs.calls.snps.filt.log'

    shell : '''

  gatk VariantFiltration -R {input.ref} -V {input.vcf} \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 200.0" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -20.0" \
    -filter-name "SOR_filter" -filter "SOR > 10.0" \
      -O {output} > {log} 2>&1 '''

# select only SNPs from the calls
rule select :
    input :
        ref = reference,
        vcf = '{path}/SRR{num}.trim:{lib}.bwa.filt.mark.rgs.calls.vcf.gz'

    output : '{path}/SRR{num,[0-9]+}.trim:{lib}.bwa.filt.mark.rgs.calls.snps.vcf.gz'
    log : '{path}/SRR{num}.trim:{lib}.bwa.filt.mark.rgs.calls.snps.log'

    shell : '''

  gatk SelectVariants -R {input.ref} -V {input.vcf} \
    -select-type SNP -O {output} > {log} 2>&1 '''

# call variants with GATK HaplotypeCaller
rule call :
    input :
        ref = reference,
        fai = reference + '.fai',
        dic = rdict,
        bam = '{path}/SRR{num}.trim:{lib}.bwa.filt.mark.rgs.bam'

    output : '{path}/SRR{num,[0-9]+}.trim:{lib}.bwa.filt.mark.rgs.calls.vcf.gz'
    log : '{path}/SRR{num}.trim:{lib}.bwa.filt.mark.rgs.calls.log'

    shell : '''

  gatk --java-options "-Xmx4g" HaplotypeCaller \
    -R {input.ref} -I {input.bam} -O {output} > {log} 2>&1 '''

rule faidx :
    input : reference
    output : reference + '.fai'
    shell : 'samtools faidx {input} && touch {output}'

rule fadict :
    input : reference
    output : rdict
    shell : 'samtools dict {input} > {output}'        

# prepare bam file for GATK: add readgroups
rule add_readgroups :
    input :
        prgm = 'scripts/picard.jar',
        bam = '{path}/SRR{num}.trim:{lib}.bwa.filt.mark.bam',
    output : '{path}/SRR{num,[0-9]+}.trim:{lib}.bwa.filt.mark.rgs.bam',
    log : '{path}/SRR{num}.trim:{lib}.bwa.filt.mark.rgs.log',

    shell : '''

  java -jar scripts/picard.jar AddOrReplaceReadGroups \
    I={input.bam} O={output} \
    SORT_ORDER=coordinate \
    RGID=id1 RGLB=lib1 RGPL=illumina RGSM=sample1 RGPU=unit1 \
    CREATE_INDEX=True > {log} 2>&1 '''

# obtain picard
picardlink = 'https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar'
rule get_picard :
    output : 'scripts/picard.jar'
    shell : 'wget -O {output} {picardlink}'

# align sequences, filter and remove duplicates
#----------------------------------------------------------------------
            
# remove duplicate reads with markdup
rule markdup :
    input : '{path}/SRR{num}.trim:{lib}.bwa.filt.bam'
    output : '{path}/SRR{num,[0-9]+}.trim:{lib}.bwa.filt.mark.bam'
    log : '{path}/SRR{num}.trim:{lib}.bwa.filt.mark.log'
    shell : 'samtools markdup -r {input} {output}'

# filter unmapped and secondary alignments
rule filter :
    input : '{path}/SRR{num}.trim:{lib}.bwa.bam'
    output : '{path}/SRR{num,[0-9]+}.trim:{lib}.bwa.filt.bam'
    log : '{path}/SRR{num}.trim:{lib}.bwa.filt.log'
    shell : '''

  samtools view -h {input} -F 260 2> {log} \
    | samtools sort -o {output} - >> {log} 2>&1 '''

# align sequences with bwa and store directly as bam
rule align :
    input :
        ref = reference,
        ind = reference + '.bwt',
        reads = '{path}/SRR{num}.trim:{lib}.fastq.gz'

    output : '{path}/SRR{num,[0-9]+}.trim:{lib}.bwa.bam'

    log :
        log = '{path}/SRR{num}.trim:{lib}.bwa.log',
        time = '{path}/SRR{num}.trim:{lib}.bwa.time',

    shell : '''

  /usr/bin/time -vo {log.time} \
    bwa mem {input.ref} {input.reads} 2> {log.log} \
      | samtools sort -o {output} - >> {log.log} 2>&1 '''

rule index :
    input : reference
    output : reference + '.bwt'
    shell : 'bwa index {input} && touch {output}'

# obtain fasta from SRA, compress and trim
#----------------------------------------------------------------------
            
# trim adaptors from a gzipped fastq file using trimmomatic
rule trim :
    input :
        prg = trim,
        lib = '{lib}.fa',
        fa = '{path}/SRR{num}.fastq.gz'

    output : '{path}/SRR{num,[0-9]+}.trim:{lib}.fastq.gz'

    log :
        log = '{path}/SRR{num}.trim:{lib}.fastq.log',
        err = '{path}/SRR{num}.trim:{lib}.fastq.err'

    shell : '''

  java -jar {input.prg} SE -phred33 {input.fa} {output} -trimlog {log.log} \
    ILLUMINACLIP:{wildcards.lib}.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > {log.err} 2>&1 '''

# gzip a fastq file
rule gzip :
    input : '{path}/SRR{num}.fastq'
    output : '{path}/SRR{num,[0-9]+}.fastq.gz'

    shell : 'gzip {input} && touch {output}'

# download a fastq file with fastqdump
rule fastqdump :
    output : '{path}/SRR{num,[0-9]+}.fastq'
    log : '{path}/SRR{num}.fastqdump.log'

    shell : '''

  fastq-dump SRR{wildcards.num} --outdir {data} > {log} 2>&1
  touch {output} '''

# trimmomatic
#----------------------------------------------------------------------

# setup library
rule setup_library :
    input : trim
    output : '{lib,.+-(SE|PE)}.fa'
    shell : 'cp trimmomatic/Trimmomatic-0.36/adapters/{output} {output}'

# unpack trimmomatic
rule unpack_trimmomatic :
    input : 'trimmomatic/Trimmomatic-0.36.zip'
    output : trim
    shell : '''

  cd trimmomatic && unzip Trimmomatic-0.36.zip && cd ..
  touch {output} '''

# obtain trimmomatic
rule get_trimmomatic :
    output : 'trimmomatic/Trimmomatic-0.36.zip'
    shell : '''

  cd trimmomatic && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && cd ..
  touch {output} '''

# collect some statistics
#----------------------------------------------------------------------

# quality check with fastqc
rule fastqc :
    input : '{path}.fastq.gz'
    output : '{path}_fastqc.html'
    log : '{path}_fastqc.log'
             
    shell : 'fastqc {input} > {log} 2>&1 && touch {output}'
