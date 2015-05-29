shell.prefix("set -euo pipefail;")  
# NOTE: seqclean is a pain in the ass when working with it. The easiest solution
# is to download it, use ldd to know which are the missing libraries to execute
# all the components, and add all the binaries and scripts to the path

# Other requirements are
#   gzip - to allow parallel decompression and specially parallel compression
#   gnu parallel - to remove bottlenecks


# List variables
SAMPLES = "test1 test2".split()
ENDS = "1 2 U".split()

# Folder variables TODO
raw      = "data/fastq_raw"
trimmed  = "data/fastq_trimmed"
assembly = "data/assembly"


# ADAPTORS
ADAPTORS = "./src/trimmomatic-0.33/adapters/TruSeq3-PE-2.fa"

# Path to programs (or element on path)
trinity     = "./src/trinityrnaseq-2.0.6/Trinity"
trimmomatic = "java -jar ./src/Trimmomatic-0.33/trimmomatic-0.33.jar "

rule all:
    """
    Generate all individual assemblies and the one of all the combined reads
    """
    input:
        expand("data/assembly/{sample}.fasta", sample= SAMPLES),
        expand("data/assembly/{sample}_normalized.fasta", sample= SAMPLES),
        "data/assembly/all.fasta",
        "data/assembly/all_normalized.fasta"



rule clean:
    """
    Delete everything
    """
    shell:
        """
        rm -rf data/fastq_trimmed
        rm -rf data/assembly
        rm -rf data/univec
        """



rule create_univec:
    """
    If not present, download the latest version of UniVec
    """
    output:
        fasta = "data/univec/univec.fasta",
        nhr   = "data/univec/univec.fasta.nhr",
        nin   = "data/univec/univec.fasta.nin",
        nsq   = "data/univec/univec.fasta.nsq",
    shell:
        """
        mkdir -p data/univec
        wget -O {output.fasta} ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
        makeblastdb                 \
            -in     {output.fasta}  \
            -dbtype nucl
        """



rule until_trimmomatic1:
    """
    Run until the first trimmomatic step (1st on trimming)
    """
    input:
        expand("data/fastq_trimmed/{sample}_{end}.trimmed1.fastq.gz",
            sample= SAMPLES,
            end= ENDS)



rule until_seqclean:
    """
    Run all until the seqclean steps (2nd on trimming)
    """
    input:
        expand("data/fastq_trimmed/{sample}_{end}.trimmed2.fastq.gz",
            sample= SAMPLES,
            end= ENDS)



rule until_trimmomatic2:
    """
    Run all until the 2nd trimmomatic step (3rd on trimming)
    """
    input:
        expand("data/fastq_trimmed/{sample}_{end}.trimmed3.fastq.gz",
            sample= SAMPLES,
            end= ENDS)



rule until_pairer:
    """
    Run all until the re-pairing process step (4th on trimming and previous
    to assembly)
    """
    input:
        expand("data/fastq_trimmed/{sample}_{end}.trimmed4.fastq.gz",
            sample= SAMPLES,
            end= ENDS)



rule only_trinity_per_sample:
    """
    Only generate the per sample assemblies
    """
    input:
        expand("data/assembly/{sample}.fasta", 
            sample= SAMPLES)



rule only_trinity_together:
    """
    Generate the combined assembly
    """
    input:
        "data/assembly/all.fasta"



rule only_trinity_normalized_per_sample:
    input:
        expand("data/assembly/{sample}.fasta",
            sample= SAMPLES)



rule only_trinity_normalized_together:
    input:
       "data/assembly/all_normalized.fasta"



rule trimmomatic1:
    """
    Run trimmomatic on paired end mode to eliminate Illumina adaptors and remove
    low quality regions and reads. 
    """
    input:
        forward = "data/fastq_raw/{sample}_1.fastq.gz",
        reverse = "data/fastq_raw/{sample}_2.fastq.gz"
    output:
        forward  = "data/fastq_trimmed/{sample}_1.trimmed1.fastq.gz",
        reverse  = "data/fastq_trimmed/{sample}_2.trimmed1.fastq.gz",
        unpaired = "data/fastq_trimmed/{sample}_U.trimmed1.fastq.gz"
    params:
        illumina_clip = "{wildcards.ADAPTORS}:2:30:10",
        unpaired_1    = "data/fastq_trimmed/{sample}_3.trimmed1.fastq.gz",
        unpaired_2    = "data/fastq_trimmed/{sample}_4.trimmed1.fastq.gz"
    log:
        trimlog  = "data/fastq_trimmed/{sample}.trimmed1_PE.log.gz",
        messages = "data/fastq_trimmed/{sample}.trimmed1_PE.results" 
    threads:
        4 # Raise it to prevent excessive RAM usage
    shell:
        """
        mkdir -p data/fastq_trimmed/
        
        {trimmomatic} PE                            \
            -threads 4                              \
            -phred33                                \
            -trimlog >( gzip -9 > {log.trimlog} )   \
            {input.forward}                         \
            {input.reverse}                         \
            {output.forward}                        \
            {params.unpaired_1}                     \
            {output.reverse}                        \
            {params.unpaired_2}                     \
            HEADCROP:13                             \
            ILLUMINACLIP:{params.illumina_clip}     \
            MINLEN:31                               \
            AVGQUAL:10                              \
            MINLEN:31                               \
            TRAILING:19                             \
            MINLEN:31                               \
            TOPHRED33                               \
            2> {log.messages}
            
            gzip -dc {params.unpaired_1} {params.unpaired_2} |
            gzip -9 > {output.unpaired}
            
            rm {params.unpaired_1} {params.unpaired_2}
        """
        


rule seqclean:
    """
    Run seqclean to remove additional contaminants.
    Parity of the reads is lost
    """
    input:
        fastq  = "data/fastq_trimmed/{sample}_{end}.trimmed1.fastq.gz",
        univec = "data/univec/univec.fasta",
        nhr    = "data/univec/univec.fasta.nhr",
        nin    = "data/univec/univec.fasta.nin",
        nsq    = "data/univec/univec.fasta.nsq"
    output:
        "data/fastq_trimmed/{sample}_{end}.trimmed2.fastq.gz"
    log:
        "data/fastq_trimmed/{sample}_{end}.trimmed2.log"
    params:
        fasta_pre  = "data/fastq_trimmed/{sample}_{end}.trimmed1.fasta",
        qual_pre   = "data/fastq_trimmed/{sample}_{end}.trimmed1.qual",
        fasta_post = "data/fastq_trimmed/{sample}_{end}.trimmed2.fasta",
        qual_post  = "data/fastq_trimmed/{sample}_{end}.trimmed2.qual",
        clean      = "data/fastq_trimmed/{sample}_{end}.trimmed2.cln"
    threads:
        16 # Adjust it to threads/2 +1 to prevent two or more seqcleans working in parallel
    shell:
        """
        # Decompress fastqs and convert it to fasta and qual.
        # I could use tee instead of decompressing twice
        gzip -dc {input.fastq}                          |
        parallel --pipe --keep-order --max-lines 4000  \
        python scripts/fastq_to_fasta.py > {params.fasta_pre} &
        
        gzip -dc {input.fastq}                          |
        parallel --pipe --keep-order --max-lines 4000  \
        python scripts/fastq_to_qual.py  > {params.qual_pre}  &
        wait
        
        # Check if fasta file is empty (and therefore the )
        if [ ! -s {params.fasta_pre} ] ; then
            rm  {params.fasta_pre}   \
                {params.qual_pre}
            touch {output}
            exit 0
        fi
            
        # Trim with seqcelan
        seqclean                    \
            {params.fasta_pre}      \
            -c 16                   \
            -o {params.fasta_post}  \
            -r {params.clean}       \
            -v {input.univec}       \
            -M                      \
            -l 31                   \
        > {log} 2>&1

        # Clean folders
        rm -rf                  \
            cleaning_*          \
            err_seqcl*.log      \
            outparts_cln.sort   \
            *.cidx              \
            formatdb.log
        
        # Generate the new qual file and rename the output
        cln2qual {params.clean} {params.qual_pre}
        mv {params.qual_pre}.clean {params.qual_post}
        
        # Generate the new fastq file
        python scripts/fasta_qual_to_fastq.py   \
            {params.fasta_post}                 \
            {params.qual_post}                  |
        gzip -9 > {output}
        
        # Clean
        rm  {params.fasta_pre}  \
            {params.qual_pre}   \
            {params.fasta_post} \
            {params.qual_post}  \
            {params.clean}      \
            seqcl*.log
        """



rule trimmomatic2:
    """
    Run trimmomatic once again, this time on single end mode, to remove low
    quality bases
    """
    input:
        "data/fastq_trimmed/{sample}_{end}.trimmed2.fastq.gz"
    output:
        "data/fastq_trimmed/{sample}_{end}.trimmed3.fastq.gz"
    log:
        trimlog  = "data/fastq_trimmed/{sample}_{end}.trimmed3.log.gz",
        messages = "data/fastq_trimmed/{sample}_{end}.trimmed3.results" 

    threads:
        8 # To prevent excessive RAM usage
    shell:
        """
        {trimmomatic} SE                            \
            -threads 4                              \
            -phred33                                \
            -trimlog >( gzip -9 > {log.trimlog} )   \
            {input}                                 \
            {output}                                \
            MINLEN:31                               \
            AVGQUAL:10                              \
            MINLEN:31                               \
            TRAILING:19                             \
            MINLEN:31                               \
            TOPHRED33                               \
            2> {log.messages}
        """



rule pairer:
    """
    Obtain once again the paired end structure in FASTQ files.
    """
    input:
        forward=  "data/fastq_trimmed/{sample}_1.trimmed3.fastq.gz",
        reverse=  "data/fastq_trimmed/{sample}_2.trimmed3.fastq.gz",
        unpaired= "data/fastq_trimmed/{sample}_U.trimmed3.fastq.gz"
    output:
        forward=  "data/fastq_trimmed/{sample}_1.trimmed4.fastq.gz",
        reverse=  "data/fastq_trimmed/{sample}_2.trimmed4.fastq.gz",
        unpaired= "data/fastq_trimmed/{sample}_U.trimmed4.fastq.gz"
    params:
        aux_1= "data/fastq_trimmed/{sample}_3.trimmed4.fastq.gz",
        aux_2= "data/fastq_trimmed/{sample}_4.trimmed4.fastq.gz"
    threads:
        4 # Heavy ram usage. If not sure, use al threads for security 
        
    shell:
        """
        python scripts/pairer.py                \
            <( gzip -dc    {input.forward}  )   \
            <( gzip -dc    {input.reverse}  )   \
            <( gzip -dc    {input.forward}  )   \
            <( gzip -dc    {input.reverse}  )   \
            >( gzip -9  >  {output.forward} )   \
            >( gzip -9  >  {output.reverse} )   \
            >( gzip -1  >  {params.aux_1}   )   \
            >( gzip -1  >  {params.aux_2}   )
        
        gzip -dc {input.unpaired} {params.aux_1} {params.aux_2} |
        gzip -9 > {output.unpaired}
        
        rm {params.aux_1} {params.aux_2}
        """



rule trinity_per_sample:
    """
    Run Trinity
    """
    input:
        forward=  "data/fastq_trimmed/{sample}_1.trimmed4.fastq.gz",
        reverse=  "data/fastq_trimmed/{sample}_2.trimmed4.fastq.gz",
        unpaired= "data/fastq_trimmed/{sample}_U.trimmed4.fastq.gz"
    output:
        assembly= "data/assembly/{sample}.fasta"
    params:
        out_tmp= "data/assembly/{sample}_trinity",
        memory=  "20G",
        threads= "24"
    log:
        "data/assembly/{sample}.log"
    threads:
        24
    shell:
        """
        mkdir -p data/assembly
        
        {trinity}                                   \
            --seqType       fq                      \
            --max_memory    {params.memory}         \
            --left          {input.forward}         \
            --right         {input.reverse}         \
            --single        {input.unpaired}        \
            --CPU           {params.threads}        \
            --full_cleanup                          \
            --output        {params.out_tmp}        \
        > {log} 2>&1
        
        mv {params.out_tmp}.Trinity.fasta {output.assembly}
        """



rule merge_for_trinity_together:
    """
    Put all the samples in one file (one per type of end)
    """
    input:
        forward=  expand("data/fastq_trimmed/{sample}_1.trimmed4.fastq.gz",
            sample= SAMPLES),
        reverse=  expand("data/fastq_trimmed/{sample}_2.trimmed4.fastq.gz",
            sample= SAMPLES),
        unpaired= expand("data/fastq_trimmed/{sample}_U.trimmed4.fastq.gz",
            sample= SAMPLES)
    output:
        forward=  "data/fastq_trimmed/all_1.trimmed4.fastq.gz",
        reverse=  "data/fastq_trimmed/all_2.trimmed4.fastq.gz",
        unpaired= "data/fastq_trimmed/all_U.trimmed4.fastq.gz"
    threads:
        3
    shell:
        """
        gzip -dc {input.forward}  | gzip -9 > {output.forward}  &
        gzip -dc {input.reverse}  | gzip -9 > {output.reverse}  &
        gzip -dc {input.unpaired} | gzip -9 > {output.unpaired} &
        wait
        """
        


rule trinity_together:
    """
    Run Trinity with all samples together
    """
    input:
        forward=  "data/fastq_trimmed/all_1.trimmed4.fastq.gz",
        reverse=  "data/fastq_trimmed/all_2.trimmed4.fastq.gz",
        unpaired= "data/fastq_trimmed/all_U.trimmed4.fastq.gz"
    output:
        assembly= "data/assembly/all.fasta"
    params:
        out_tmp=  "data/assembly/all_trinity",
        memory=   "20G",
        threads=  "24"
    log:
        "data/assembly/all.log"
    threads:
        24
    shell:
        """
        mkdir -p data/assembly
        
        {trinity}                                   \
            --seqType       fq                      \
            --max_memory    {params.memory}         \
            --left          {input.forward}         \
            --right         {input.reverse}         \
            --single        {input.unpaired}        \
            --CPU           {params.threads}        \
            --full_cleanup                          \
            --output        {params.out_tmp}        \
        > {log} 2>&1
        
        mv {params.out_tmp}.Trinity.fasta {output.assembly}
        """



rule trinity_normalized_per_sample:
    """
    Run Trinity with the normalization step in each sample
    """
    input:
        forward=  "data/fastq_trimmed/{sample}_1.trimmed4.fastq.gz",
        reverse=  "data/fastq_trimmed/{sample}_2.trimmed4.fastq.gz",
        unpaired= "data/fastq_trimmed/{sample}_U.trimmed4.fastq.gz"
    output:
        assembly= "data/assembly/{sample}_normalized.fasta"
    params:
        out_tmp= "data/assembly/{sample}_normalized_trinity",
        memory=  "20G",
        threads= "24"
    log:
        "data/assembly/{sample}_normalized.log"
    threads:
        24
    shell:
        """
        mkdir -p data/assembly
        
        {trinity}                                   \
            --seqType       fq                      \
            --normalize_reads                       \
            --max_memory    {params.memory}         \
            --left          {input.forward}         \
            --right         {input.reverse}         \
            --single        {input.unpaired}        \
            --CPU           {params.threads}        \
            --full_cleanup                          \
            --output        {params.out_tmp}        \
        > {log} 2>&1
        
        mv {params.out_tmp}.Trinity.fasta {output.assembly}
        """



rule trinity_normalized_together:
    """
    Run Trinity with the normalization step in the combined assembly
    """
    input:
        forward=  "data/fastq_trimmed/all_1.trimmed4.fastq.gz",
        reverse=  "data/fastq_trimmed/all_2.trimmed4.fastq.gz",
        unpaired= "data/fastq_trimmed/all_U.trimmed4.fastq.gz"
    output:
        assembly= "data/assembly/all_normalized.fasta"
    params:
        out_tmp=  "data/assembly/all_normalized_trinity",
        memory=   "20G",
        threads=  "24"
    log:
        "data/assembly/all_normalized.log"
    threads:
        24
    shell:
        """
        mkdir -p data/assembly
        
        {trinity}                                   \
            --normalize_reads                       \
            --seqType       fq                      \
            --max_memory    {params.memory}         \
            --left          {input.forward}         \
            --right         {input.reverse}         \
            --single        {input.unpaired}        \
            --CPU           {params.threads}        \
            --full_cleanup                          \
            --output        {params.out_tmp}        \
        > {log} 2>&1
        
        mv {params.out_tmp}.Trinity.fasta {output.assembly}
        """

