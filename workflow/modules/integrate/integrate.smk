# ==============================================================================
# MODULE: integrate
# Merge every sample's Seurat object, then Harmony on the RNA and ADT PCAs,
# weighted nearest neighbours (WNN) and a WNN UMAP.
# ==============================================================================

INT_SCRIPTS = os.path.join(MODULES_DIR, "integrate", "scripts")
INT_ENV = r_env("integrate")


rule merge:
    """Merge every usable sample's Seurat object; cell names become {sample}_{barcode}."""
    input:
        seurat_objects=[P.result("preprocess", "construct", d, "objects", f"{s}_Seurat.qs2") for d, s in USABLE],
    output:
        merged_seurat=P.result("integrate", "merge", "merged_seurat.qs2"),
    resources:
        merge_jobs=1,
    log:
        P.log("integrate", "merge"),
    benchmark:
        P.benchmark("integrate", "merge")
    conda:
        INT_ENV.conda
    container:
        INT_ENV.container
    script:
        os.path.join(INT_SCRIPTS, "merge.R")


rule harmony:
    """Harmony on the RNA and ADT PCAs, then WNN graphs (wsnn_harmony) and the wnn.umap_harmony UMAP."""
    input:
        merged_seurat=rules.merge.output.merged_seurat,
    output:
        harmonized_seurat=CTX.integrated_object,
    params:
        harmony=config["integrate"]["harmony"],
        wnn=config["integrate"]["wnn"],
        isotype_pattern=config["preprocess"]["dsb"]["isotype_pattern"],
    log:
        P.log("integrate", "harmony"),
    benchmark:
        P.benchmark("integrate", "harmony")
    conda:
        INT_ENV.conda
    container:
        INT_ENV.container
    script:
        os.path.join(INT_SCRIPTS, "harmony.R")
