
data = '/data/clownfish/'
gaboriau = ['SRR295227{}'.format(x) for x in range(58, 77)]
anothernine = ['SRR84426{}'.format(x) for x in range(10, 21)]
srrs = gaboriau + anothernine + ['SRR26235419', 'SRR6685843']

#----------------------------------------------------------------------

rule master :
    input :
        expand(data + '{s}.fastq', s = srrs)

#----------------------------------------------------------------------

# download a fastq file with fastqdump
rule fastqdump :
    params :
        lambda wildcards : wildcards.path.rsplit('/', 1)[1]

    output : '{path}.fastq'
    log : '{path}.fastqdump.log'

    shell : '''

  fastq-dump {params} --outdir {data} > {log} 2>&1
  touch {output} '''

