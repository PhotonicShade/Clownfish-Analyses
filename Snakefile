
data = '/data/clownfish/'
gaboriau = ['SRR295227{}'.format(x) for x in range(58, 77)]
anothernine = ['SRR84426{}'.format(x) for x in range(10, 21)]
srrs = gaboriau + anothernine + ['SRR26235419', 'SRR6685843']

trim = 'trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar'
alib = 'TruSeq3-PE.fa'

#----------------------------------------------------------------------

rule master :
    input :
        expand(data + '{s}.fastq.trim.gz', s = srrs)

#----------------------------------------------------------------------

# trim adaptors from a gzipped fastq file using trimmomatic
rule trim :
    input :
        prg = trim,
        lib = alib,
        fa = '{path}.fastq.gz'

    output : '{path}.fastq.trim.gz'

    log :
        log = '{path}.fastq.trim.log',
        err = '{path}.fastq.trim.err'

    shell : '''

  java -jar {input.prg} SE -phred33 {input.fa} {output} -trimlog {log.log} \
    ILLUMINACLIP:{input.lib}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > {log.err} 2>&1 '''

# gzip a fastq file
rule gzip :
    input : '{path}.fastq'
    output : '{path}.fastq.gz'

    shell : 'gzip {input} && touch {output}'

# download a fastq file with fastqdump
rule fastqdump :
    params :
        lambda wildcards : wildcards.path.rsplit('/', 1)[1]

    output : '{path}.fastq'
    log : '{path}.fastqdump.log'

    shell : '''

  fastq-dump {params} --outdir {data} > {log} 2>&1
  touch {output} '''

# trimmomatic
#----------------------------------------------------------------------

# setup library
rule setup_library :
    input : trim
    output : alib
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
