rule make_samplesheet:
    output:
        tsv = RES_IDX['setup/samplesheet']['outdir'] + "samplesheet_{dataset}_{modality}.tsv"
    run:
        import os
        import glob
        import csv
        import re
        from collections import defaultdict

        ds = config["datasets"][wildcards.dataset]

        filetype = ds["filetype"]
        modality = wildcards.modality.upper()

        data = defaultdict(lambda: defaultdict(dict))
        full_sample_names = {}
        sample_dirs = {}
        all_lanes = set()

        files = []

        # -------------------------------------------------
        # FIX: store actual directory path (value of input_dirs)
        # -------------------------------------------------
        for dir_name, d in ds["input_dirs"].items():
            files.extend([(f, d) for f in glob.glob(os.path.join(d, f"*.{filetype}"))])

        if not files:
            raise ValueError("No FASTQs found")

        mod_tag = modality

        for f, dir_path in files:
            name = os.path.basename(f)

            if mod_tag not in name:
                continue

            # -------------------------------------------------
            # isolate AFTER modality
            # -------------------------------------------------
            try:
                after_mod = name.split(mod_tag, 1)[1]
            except:
                continue

            after_mod = after_mod.lstrip("-_")

            # -------------------------------------------------
            # sample key
            # -------------------------------------------------
            lane_match = re.search(r"_L\d{3}_", after_mod)
            if not lane_match:
                continue

            sample = after_mod[:lane_match.start()]
            sample = sample.strip("-_")

            # -------------------------------------------------
            # FULL SAMPLE NAME
            # -------------------------------------------------
            lane_boundary = re.search(r"_L\d{3}_R[12]", name)
            if not lane_boundary:
                continue

            full_sample = name[:lane_boundary.start()].rstrip("-_")

            if sample not in full_sample_names:
                full_sample_names[sample] = full_sample

            # store actual directory path
            if sample not in sample_dirs:
                sample_dirs[sample] = dir_path

            # -------------------------------------------------
            # lane + read
            # -------------------------------------------------
            lane = re.search(r"(L\d{3})", after_mod)
            read = re.search(r"(R[12])", after_mod)

            if not lane or not read:
                continue

            lane = lane.group(1)
            read = read.group(1)

            data[sample][lane][read] = f
            all_lanes.add(lane)

        # -------------------------------------------------
        # WRITE TSV
        # -------------------------------------------------
        sorted_lanes = sorted(all_lanes)

        header = ["SampleName", "FullSampleName", "InputDir"]
        for lane in sorted_lanes:
            header += [f"{lane}_R1", f"{lane}_R2"]

        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)

        with open(output.tsv, "w", newline="") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(header)

            for sample in sorted(data.keys()):
                row = [
                    sample,
                    full_sample_names.get(sample, ""),
                    sample_dirs.get(sample, "")
                ]

                lanes_dict = data[sample]

                for lane in sorted_lanes:
                    reads = lanes_dict.get(lane, {})
                    row.append(reads.get("R1", ""))
                    row.append(reads.get("R2", ""))

                writer.writerow(row)