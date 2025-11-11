
data = '/data/clownfish/'
gaboriau = ['SRR295227{}'.format(x) for x in range(58, 77)]
anothernine = ['SRR84426{}'.format(x) for x in range(10, 21)]
srrs = gaboriau + anothernine + ['SRR26235419', 'SRR6685843']

trim = 'trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar'
libs = ['TruSeq3-PE', 'TruSeq3-SE']

reference = data + 'A_frenatus/afrenatus.transcripts.uniprot.fa'

#----------------------------------------------------------------------

rule master :
    input :
        # trimmed readsets
        expand(data + '{s}.trim:{lib}.bwa.filt.mark.bam', s = srrs, lib = libs),

        # read qualities
        expand(data + '{s}_fastqc.html', s = srrs),
        expand(data + '{s}.trim:{lib}_fastqc.html', s = srrs, lib = libs),

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
