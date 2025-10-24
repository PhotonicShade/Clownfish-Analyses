
data = '/data/clownfish/'
gaboriau = ['SRR295227{}'.format(x) for x in range(58, 77)]
anothernine = ['SRR84426{}'.format(x) for x in range(10, 21)]
srrs = gaboriau + anothernine + ['SRR26235419', 'SRR6685843']

trim = 'trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar'
libs = ['TruSeq3-PE', 'TruSeq3-SE']

#----------------------------------------------------------------------

rule master :
    input :
        # trimmed readsets
        expand(data + '{s}.trim:{lib}.fastq.gz', s = srrs, lib = libs),

        # read qualities
        expand(data + '{s}_fastqc.html', s = srrs),
        expand(data + '{s}.trim:{lib}_fastqc.html', s = srrs, lib = libs),

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
    output : '{lib}.fa'
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
