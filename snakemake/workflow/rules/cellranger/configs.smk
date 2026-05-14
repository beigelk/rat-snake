checkpoint make_configs:
    input:
        gex=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_gex.tsv",
        csp=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_csp.tsv",
        vdj=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_vdj.tsv",
    output:
        config=RES_IDX['cellranger/configs']['outdir'] + "{dataset}/configfiles/{sample}_config.csv",
        skipped=RES_IDX['cellranger/configs']['outdir'] + "{dataset}/tracking/{sample}_record.tsv",
    run:
        import csv
        import re

        def load_tsv(path):
            index = {}
            with open(path) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    index[row["SampleName"]] = {
                        # "fastq_id": row["FullSampleName"],
                        "fastq_id": re.sub(r'_S\d+$', '', row["FullSampleName"]),
                        "fastqs":   row["InputDir"]
                    }
            return index

        gex_tsv = load_tsv(input.gex)
        csp_tsv = load_tsv(input.csp)
        vdj_tsv = load_tsv(input.vdj)

        tsv_map     = {"gex": gex_tsv, "csp": csp_tsv, "vdj": vdj_tsv}
        feature_map = {"gex": "Gene Expression", "csp": "Antibody Capture", "vdj": "VDJ-T"}

        sample_key = wildcards.sample
        sample     = SAMPLES[wildcards.dataset][sample_key]

        libraries = []
        skipped   = []

        for modality, feature in feature_map.items():

            if modality not in sample:
                skipped.append({"sample": sample_key, "modality": modality, "result": "SKIPPED", "reason": "not in SAMPLES"})
                continue

            mod_id = sample[modality]
            tsv    = tsv_map[modality]

            if mod_id not in tsv:
                skipped.append({"sample": sample_key, "modality": modality, "result": "SKIPPED", "reason": f"ID '{mod_id}' not in samplesheet"})
                continue

            entry = tsv[mod_id]
            libraries.append({
                "fastq_id":      entry["fastq_id"],
                "fastqs":        entry["fastqs"],
                "lanes":         "any",
                "feature_types": feature
            })

        with open(output.skipped, "w") as f:
            f.write("sample\tmodality\tresult\treason\n")
            for row in skipped:
                f.write(f"{row['sample']}\t{row['modality']}\t{row['result']}\t{row['reason']}\n")

        with open(output.config, "w") as out:
            out.write("[gene-expression]\n")
            out.write(f"reference,{config['cellranger']['gex_reference']}\n")
            out.write("create-bam,false\n")
            out.write("check-library-compatibility,false\n\n")

            out.write("[feature]\n")
            out.write(f"reference,{config['cellranger']['csp_reference']}\n\n")

            out.write("[vdj]\n")
            out.write(f"reference,{config['cellranger']['vdj_reference']}\n\n")

            out.write("[libraries]\n")
            out.write("fastq_id,fastqs,lanes,feature_types\n")
            for lib in libraries:
                out.write(
                    f"{lib['fastq_id']},"
                    f"{lib['fastqs']},"
                    f"{lib['lanes']},"
                    f"{lib['feature_types']}\n"
                )

# rule make_configs:
#     input:
#         gex=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_gex.tsv",
#         csp=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_csp.tsv",
#         vdj=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_vdj.tsv",
#     output:
#         RES_IDX['cellranger/configs']['outdir'] + "cellranger_configs/{dataset}/{sample}_config.csv",
#     run:
#         import csv

#         def load_tsv(path):
#             """Returns dict keyed by FullSampleName -> {fastq_id, fastqs}"""
#             index = {}
#             with open(path) as f:
#                 reader = csv.DictReader(f, delimiter="\t")
#                 for row in reader:
#                     index[row["SampleName"]] = {
#                         "fastq_id": row["FullSampleName"],
#                         "fastqs":   row["InputDir"]
#                     }
#             return index

#         gex_tsv = load_tsv(input.gex)
#         csp_tsv = load_tsv(input.csp)
#         vdj_tsv = load_tsv(input.vdj)

#         tsv_map = {"gex": gex_tsv, "csp": csp_tsv, "vdj": vdj_tsv}

#         feature_map = {
#             "gex": "Gene Expression",
#             "csp": "Antibody Capture",
#             "vdj": "VDJ"
#         }

#         sample_key = wildcards.sample
#         sample = SAMPLES[wildcards.dataset][sample_key]

#         libraries = []

#         for modality, feature in feature_map.items():
#             if modality not in sample:
#                 continue

#             mod_id = sample[modality]  # e.g. '25_13_1205_S25'
#             tsv = tsv_map[modality]

#             print("MODALITY: " + modality)
#             print(tsv)

#             if mod_id not in tsv:
#                 raise ValueError(f"{modality} ID '{mod_id}' not found in {modality} samplesheet")

#             entry = tsv[mod_id]
#             libraries.append({
#                 "fastq_id":     entry["fastq_id"],
#                 "fastqs":       entry["fastqs"],
#                 "lanes":        "any",
#                 "feature_types": feature
#             })

#         with open(output[0], "w") as out:
#             out.write("[gene-expression]\n")
#             out.write(f"reference,{config['cellranger']['gex_reference']}\n")
#             out.write("create-bam,false\n\n")

#             out.write("[feature]\n")
#             out.write(f"reference,{config['cellranger']['csp_reference']}\n\n")

#             out.write("[vdj]\n")
#             out.write(f"reference,{config['cellranger']['vdj_reference']}\n\n")

#             out.write("[libraries]\n")
#             for lib in libraries:
#                 out.write(
#                     f"{lib['fastq_id']},"
#                     f"{lib['fastqs']},"
#                     f"{lib['lanes']},"
#                     f"{lib['feature_types']}\n"
#                 )

# # rule make_configs:
# #     input:
# #         gex=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_gex.tsv",
# #         csp=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_csp.tsv",
# #         vdj=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_vdj.tsv",
# #     output:
# #         RES_IDX['cellranger/configs']['outdir'] + "cellranger_configs/{dataset}/{sample}_config.csv",
# #     run:
# #         import csv
# #         from collections import defaultdict

# #         meta = defaultdict(dict)

# #         def load_tsv(path, modality):
# #             with open(path) as f:
# #                 reader = csv.DictReader(f, delimiter="\t")
# #                 for row in reader:
# #                     meta[row["SampleName"]][modality] = {
# #                         "fastq_id": row["FullSampleName"],
# #                         "fastqs": row["InputDir"]
# #                     }

# #         load_tsv(input.gex, "GEX")
# #         load_tsv(input.csp, "CSP")
# #         load_tsv(input.vdj, "VDJ")

# #         sample = wildcards.sample

# #         if sample not in meta:
# #             raise ValueError(f"Sample {sample} not found")

# #         libraries = []
# #         feature_map = {
# #             "GEX": "Gene Expression",
# #             "CSP": "Antibody Capture",
# #             "VDJ": "VDJ"
# #         }

# #         for modality, feature in feature_map.items():
# #             if modality not in meta[sample]:
# #                 continue
# #             entry = meta[sample][modality]
# #             libraries.append({
# #                 "fastq_id": entry["fastq_id"],
# #                 "fastqs": entry["fastqs"],
# #                 "lanes": "any",
# #                 "feature_types": feature
# #             })

# #         with open(output[0], "w") as out:
# #             out.write("[gene-expression]\n")
# #             out.write(f"reference,{config['cellranger']['gex_reference']}\n")
# #             out.write("create-bam,false\n\n")

# #             out.write("[feature]\n")
# #             out.write(f"reference,{config['cellranger']['csp_reference']}\n\n")

# #             out.write("[vdj]\n")
# #             out.write(f"reference,{config['cellranger']['vdj_reference']}\n\n")

# #             out.write("[libraries]\n")
# #             for lib in libraries:
# #                 out.write(
# #                     f"{lib['fastq_id']},"
# #                     f"{lib['fastqs']},"
# #                     f"{lib['lanes']},"
# #                     f"{lib['feature_types']}\n"
# #                 )

# # # rule make_configs:
# # #     input:
# # #         # no longer needed (SAMPLES replaces TSV parsing)
# # #         gex=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_gex.tsv",
# # #         csp=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_csp.tsv",
# # #         vdj=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_vdj.tsv",

# # #     output:
# # #         RES_IDX['cellranger/configs']['outdir'] + "cellranger_configs/{dataset}/{sample_key}_config.csv",

# # #     run:
# # #         dataset = wildcards.dataset
# # #         sample_key = wildcards.sample_key

# # #         # NEW: direct access
# # #         try:
# # #             sample = SAMPLES[dataset][sample_key]
# # #         except KeyError:
# # #             raise ValueError(f"Sample {sample_key} not found in dataset {dataset}")

# # #         libraries = []

# # #         MOD_MAP = {
# # #             "gex": "Gene Expression",
# # #             "csp": "Antibody Capture",
# # #             "vdj": "VDJ"
# # #         }

# # #         FILE_MAP = {
# # #             "gex": config['cellranger']['gex_reference'],
# # #             "csp": config['cellranger']['csp_reference'],
# # #             "vdj": config['cellranger']['vdj_reference']
# # #         }

# # #         for modality, feature in MOD_MAP.items():

# # #             if modality not in sample:
# # #                 continue

# # #             entry = sample[modality]

# # #             libraries.append({
# # #                 "fastq_id": entry["fastq_id"],
# # #                 "fastqs": entry["fastqs"],
# # #                 "lanes": "any",
# # #                 "feature_types": feature
# # #             })

# # #         with open(output[0], "w") as out:

# # #             out.write("[gene-expression]\n")
# # #             out.write(f"reference,{config['cellranger']['gex_reference']}\n")
# # #             out.write("create-bam,false\n\n")

# # #             out.write("[feature]\n")
# # #             out.write(f"reference,{config['cellranger']['csp_reference']}\n\n")

# # #             out.write("[vdj]\n")
# # #             out.write(f"reference,{config['cellranger']['vdj_reference']}\n\n")

# # #             out.write("[libraries]\n")

# # #             for lib in libraries:
# # #                 out.write(
# # #                     f"{lib['fastq_id']},"
# # #                     f"{lib['fastqs']},"
# # #                     f"{lib['lanes']},"
# # #                     f"{lib['feature_types']}\n"
# # #                 )

# # # # rule make_configs:
# # # #     input:
# # # #         gex=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_gex.tsv",
# # # #         csp=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_csp.tsv",
# # # #         vdj=RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_vdj.tsv",

# # # #     output:
# # # #         RES_IDX['cellranger/configs']['outdir'] + "cellranger_configs/{dataset}/{sample}_config.csv",

# # # #     run:
# # # #         import csv
# # # #         from collections import defaultdict

# # # #         gex_path = input.gex
# # # #         csp_path = input.csp
# # # #         vdj_path = input.vdj

# # # #         sample = wildcards.sample

# # # #         meta = defaultdict(dict)

# # # #         def load_tsv(path, modality):
# # # #             with open(path) as f:
# # # #                 reader = csv.DictReader(f, delimiter="\t")
# # # #                 for row in reader:
# # # #                     meta[row["SampleName"]][modality] = {
# # # #                         "fastq_id": row["FullSampleName"],
# # # #                         "fastqs": row["InputDir"]
# # # #                     }

# # # #         load_tsv(gex_path, "GEX")
# # # #         load_tsv(csp_path, "CSP")
# # # #         load_tsv(vdj_path, "VDJ")

# # # #         if sample not in meta:
# # # #             raise ValueError(f"Sample {sample} not found")

# # # #         libraries = []

# # # #         for modality in ["GEX", "CSP", "VDJ"]:
# # # #             if modality not in meta[sample]:
# # # #                 continue

# # # #             entry = meta[sample][modality]

# # # #             feature = {
# # # #                 "GEX": "Gene Expression",
# # # #                 "CSP": "Antibody Capture",
# # # #                 "VDJ": "VDJ"
# # # #             }[modality]

# # # #             libraries.append({
# # # #                 "fastq_id": entry["fastq_id"],
# # # #                 "fastqs": entry["fastqs"],
# # # #                 "lanes": "any",
# # # #                 "feature_types": feature
# # # #             })
# # # #         print(output)
# # # #         # with open(output[0], "w") as out:
# # # #         #     out.write("[libraries]\n")
# # # #         #     for lib in libraries:
# # # #         #         out.write(
# # # #         #             f"{lib['fastq_id']},"
# # # #         #             f"{lib['fastqs']},"
# # # #         #             f"{lib['lanes']},"
# # # #         #             f"{lib['feature_types']}\n"
# # # #         #         )
# # # #         with open(output[0], "w") as out:
# # # #             # gene-expression section
# # # #             out.write("[gene-expression]\n")
# # # #             out.write(f"reference,{config['cellranger']['gex_reference']}\n")
# # # #             out.write("create-bam,false\n")
# # # #             out.write("\n")

# # # #             # feature section
# # # #             out.write("[feature]\n")
# # # #             out.write(f"reference,{config['cellranger']['vdj_reference']}\n")
# # # #             out.write("\n")

# # # #             # vdj section
# # # #             out.write("[vdj]\n")
# # # #             out.write(f"reference,{config['cellranger']['vdj_reference']}\n")
# # # #             out.write("\n")

# # # #             # libraries section
# # # #             out.write("[libraries]\n")
# # # #             for lib in libraries:
# # # #                 out.write(
# # # #                     f"{lib['fastq_id']},"
# # # #                     f"{lib['fastqs']},"
# # # #                     f"{lib['lanes']},"
# # # #                     f"{lib['feature_types']}\n"
# # # #                 )