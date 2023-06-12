rule all:
    input:
        "/mnt/c/BGAproject/data/normcounts_log2cpm_NP.csv",
        "/mnt/c/BGAproject/data/normcounts_log2cpm_PE.csv",
        "/mnt/c/BGAproject/data/GENELIST.CSV",
        "/mnt/c/BGAproject/data/NP_Pathtable.csv",
        "/mnt/c/BGAproject/data/GENELIST_PE.CSV",
        "/mnt/c/BGAproject/data/PE_Pathtable.csv",
        "/mnt/c/BGAproject/data/BP_NP.png",
        "/mnt/c/BGAproject/data/BP_PE.png"

rule run_rscript:
    input:
        rscript = "/mnt/c/BGAproject/script/rscript.R",
        data_path = "/mnt/c/BGAproject/data/"
    output:
        "/mnt/c/BGAproject/data/normcounts_log2cpm_NP.csv",
        "/mnt/c/BGAproject/data/normcounts_log2cpm_PE.csv",
        "/mnt/c/BGAproject/data/GENELIST.CSV",
        "/mnt/c/BGAproject/data/NP_Pathtable.csv",
        "/mnt/c/BGAproject/data/GENELIST_PE.CSV",
        "/mnt/c/BGAproject/data/PE_Pathtable.csv",
        "/mnt/c/BGAproject/data/BP_NP.png",
        "/mnt/c/BGAproject/data/BP_PE.png"
    shell:
        "Rscript {input.rscript}"
