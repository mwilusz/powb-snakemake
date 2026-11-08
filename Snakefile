# Indicate the config file
configfile:
    'config.yml'#--configfile nazwa pliku

# Extract query sequence names (sequence IDs) from file paths
# to use them as IDs for query files
seqids, = glob_wildcards('input/queries/{seqid}.fna')

# Extract genome assembly accession numbers (assembly IDs) from file paths
# to use them as IDs for genome files
asmids, = glob_wildcards(config['genomesdir'] + '/{asmid}.fna.gz')

# Rule all: collect the final output (multiple files)
# WARNING: initially runs the pipeline up to the rule merge_flt
rule all:
    input:
        out1='output/tree.newick',
        out2='output/tree.png'

# Rule blastn: for each genome perform a nucletide BLAST search using:
# - a single-sequence (a fragment of gmk gene) FASTA query file
# - a multiple-sequence FASTA file (one staphylococcal genome assembly contigs)
# Uncompress the gzipped genome assembly file and redirect the output
# to blastn's stdin while asking blastn to read the subject FASTA sequence
# from that stream ('-' instead of a subject file path)
rule blastn:
    conda:
        'envs/powb-blast.yml'
    params:
        evalue = 0.01,
        cols   = 'qseqid qlen sseqid length qcovs pident evalue sseq'
    input:
        query = 'input/queries/{seqid}.fna',
        subj  = config['genomesdir'] + '/{asmid}.fna.gz'
    output:
        'output/blastn/blastn_{asmid}_{seqid}.tsv'
    log:
        'log/blastn/blastn_{asmid}_{seqid}.log'
    shell:
        '''gunzip -c {input.subj} 2> {log} |      \
           blastn -evalue         {params.evalue} \
                  -outfmt      "6 {params.cols}"  \
                  -query          {input.query}   \
                  -subject        -               \
                  -out            {output}        \
                   >> {log} 2>&1
        '''

# Rule filter: for each search retain 100% query coverage hits (alignments)
# of the same length as the query sequence (gapless alignment)
rule filter:
    conda:
        'envs/powb-python.yml'
    params:
        qcovs = 100,
        cols  = rules.blastn.params.cols
    input:
        rules.blastn.output
    output:
        'output/filter/filtered_{asmid}_{seqid}.tsv'
    log:
        'log/filter/filter_{asmid}_{seqid}.log'
    shell:
        '''scripts/filter.py --eqlen                    \
                             --asmid  {wildcards.asmid} \
                             --qcovs  {params.qcovs}    \
                             --cols  "{params.cols}"    \
                             --input  {input}           \
                             --output {output}          \
                               > {log} 2>&1
        '''

# Rule merge_flt: join all filtered BLAST result files into one using
# a short Bash script (a list of file path on input)
rule merge_flt:
    params:
        mask = rules.filter.output[0].replace('{seqid}', '*')
    input:
        expand(rules.filter.output, seqid=seqids, allow_missing=True)
    output:
        'output/merge_flt/merged_{asmid}.tsv'
    log:
        'log/merge_flt/merge_{asmid}.log'
    shell:
        '''files=({params.mask})
           head -1 ${{files[0]}} > {output} 2> {log}
           for file in "${{files[@]}}"; do
               tail -n +2 ${{file}}
           done >> {output} 2>> {log}
        '''

# Rule mlst: perform MLST analysis, one genome at a time
rule mlst:
    conda:
        'envs/powb-python.yml'
    params:
        alleles='input/alleles.fna',
        bigsdb='input/bigsdb.tsv',
        asmid='{asmid}'
    input:
        rules.merge_flt.output
    output:
        'output/mlst/mlst_{asmid}.tsv'
    log:
        'log/mlst/{asmid}.log'
    shell:
        '''scripts/mlst.py --asmid {params.asmid} --alleles {input} --allelesdb {params.alleles} --mlstdb {params.bigsdb} --output {output} 2> {log}'''
    
    

# Rule merge_mlst: join all MLST results files into one using a short Bash
# script (a list of file path on input)
rule merge_mlst:
    params:
        mask = rules.mlst.output[0].replace('{asmid}', '*')
    input:
        expand(rules.mlst.output, asmid=asmids, allow_missing=True)
    output:
        'output/merged.tsv'
    log:
        'log/merged.log'
    shell:
        '''files=({params.mask})
           head -1 ${{files[0]}} > {output} 2> {log}
           for file in "${{files[@]}}"; do
               tail +2 $file >> {output} 2>> {log}
           done
        '''
# Rule unique: from merged BLAST results extract unique merged MLST loci
# sequences (i.e., one sequence per ST) to a FASTA file
rule unique:
    conda:
        'envs/powb-python.yml'
    input:
        rules.merge_mlst.output
    output:
        'output/unique.fna'
    log:
        'log/unique.log'
    shell:
        '''scripts/unique.py --input {input}\
                           --output {output} 2> {log}
        '''
    

# Rule tree: take the merged MLST results and extracted merged MLST loci
# sequences, and then compute and visualise a phylogenetic tree (NJ tree)
rule tree:
    conda:
        'envs/powb-python.yml'
    input:
        in1=rules.merge_mlst.output,
        in2=rules.unique.output
    output:
        out1='output/tree.newick',
        out2='output/tree.png'
    log:
        'log/tree.log'
    shell:
        '''scripts/tree.py --mlst     {input.in1}\
                           --seqs     {input.in2}\
                           --outtree  {output.out1}\
                           --outvis   {output.out2} 2> {log}
        '''

