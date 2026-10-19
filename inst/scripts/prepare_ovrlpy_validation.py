# Assisted-by: OpenAI Codex
"""Prepare the archived Xenium/vendor and Baysor validation inputs.

Development-time only. This script never runs during package build, check,
installation, examples, tests, vignette rendering, or package loading.

Usage: python3 prepare_ovrlpy_validation.py INPUT_DIR WORK_DIR
Requires Python 3, numpy, and pyarrow. No network access is performed.

Place these exact public archive members in INPUT_DIR, using the basenames below:
  Official 10x Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip:
    experiment.xenium, transcripts.parquet, cells.parquet,
    cell_feature_matrix/features.tsv.gz
  Zenodo record 16736589 v2, baysor_xenium-brain.zip:
    the ZIP itself, plus baysor_xenium-brain/segmentation.csv,
    segmentation_cell_stats.csv, segmentation_log.log,
    segmentation_params.dump.toml
The ZIP members may be extracted to INPUT_DIR with their basenames. Source:
https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_MultiSection_1/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs.zip
https://zenodo.org/records/16736589 (DOI 10.5281/zenodo.16736589)
The source archive S3 version ID and each required member hash are in the
bundled validation record. No file name, version, or preprocessing is inferred.
This script reads the saved segmentations; it does not rerun Baysor.
"""
import csv
import gzip
import hashlib
import json
import pathlib
import sys

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv
import pyarrow.parquet as pq


def fail(message):
    raise SystemExit(message)


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def record_values():
    path = pathlib.Path(__file__).resolve().parents[1] / "extdata" / "validation" / "ovrlpy_validation_summary.csv"
    if not path.is_file():
        fail(f"Missing bundled validation record: {path}")
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    result = {row["key"]: row["value"] for row in rows}
    if len(result) != len(rows):
        fail("Duplicate validation-record keys")
    return result


def verify_inputs(input_dir, expected):
    names = ("experiment.xenium", "transcripts.parquet", "cells.parquet",
             "features.tsv.gz", "segmentation.csv", "segmentation_cell_stats.csv",
             "segmentation_log.log", "segmentation_params.dump.toml")
    observed = {}
    for name in names:
        path = input_dir / name
        if not path.is_file():
            fail(f"Missing public archive member: {path}")
        actual = digest(path)
        if actual != expected[f"sha256_{name}"]:
            fail(f"SHA-256 mismatch for {name}: {actual}")
        observed[f"sha256_{name}"] = actual
        print(f"Verified SHA-256: {name}", flush=True)
    archive = input_dir / "baysor_xenium-brain.zip"
    if not archive.is_file():
        fail(f"Missing deposited Baysor archive: {archive}")
    archive_md5 = digest(archive, "md5")
    if archive_md5 != expected["zenodo_baysor_archive_md5"]:
        fail("Deposited Baysor ZIP MD5 mismatch")
    observed["zenodo_baysor_archive_md5"] = archive_md5
    info = json.loads((input_dir / "experiment.xenium").read_text())
    if (info.get("region_name") != expected["experiment"] or
            info.get("analysis_sw_version") != expected["vendor_analysis_software"] or
            info.get("panel_name") != expected["source_panel"] or
            info.get("panel_num_targets_predesigned") != 248):
        fail("Official experiment metadata disagrees with validation record")
    log = (input_dir / "segmentation_log.log").read_text()
    if (f"Run {expected['baysor_run_id']}" not in log or
            expected["baysor_software"] not in log):
        fail("Saved Baysor log disagrees with run identity/version")
    return observed


def source_inputs(input_dir, work_dir, expected):
    with gzip.open(input_dir / "features.tsv.gz", "rt") as handle:
        genes = {row[1].encode() for row in csv.reader(handle, delimiter="\t")
                 if row[2] == "Gene Expression"}
    if len(genes) != 248:
        fail("Source feature dictionary does not contain 248 Gene Expression targets")
    feature_set = pa.array(sorted(genes), type=pa.binary())
    source = pq.ParquetFile(input_dir / "transcripts.parquet")
    if source.num_row_groups != 63:
        fail("Unexpected source Parquet row-group count")
    # Read row groups explicitly: this source file's aggregate row-count footer
    # reports one million although its 63 groups contain over 62 million rows.
    counts = {}
    total = 0
    minx = miny = float("inf")
    maxx = maxy = float("-inf")
    universe_path = work_dir / "universe.ids.u64le"
    with universe_path.open("wb") as out:
        for group in range(source.num_row_groups):
            tab = source.read_row_group(group, columns=["transcript_id", "feature_name", "qv", "x_location", "y_location"])
            keep = pc.and_(pc.is_in(tab["feature_name"], value_set=feature_set),
                           pc.greater_equal(tab["qv"], 20)).to_numpy(zero_copy_only=False)
            ids = tab["transcript_id"].to_numpy()[keep].astype("<u8", copy=False)
            ids.tofile(out)
            total += len(ids)
            x = tab["x_location"].to_numpy()[keep]
            y = tab["y_location"].to_numpy()[keep]
            if len(x):
                minx, maxx = min(minx, float(x.min())), max(maxx, float(x.max()))
                miny, maxy = min(miny, float(y.min())), max(maxy, float(y.max()))
                tiles = np.stack((np.floor(x / 500).astype(np.int16),
                                  np.floor(y / 500).astype(np.int16)), axis=1)
                unique, sizes = np.unique(tiles, axis=0, return_counts=True)
                for (tx, ty), size in zip(unique, sizes):
                    key = (int(tx), int(ty))
                    counts[key] = counts.get(key, 0) + int(size)
    if total != int(expected["source_universe_size"]):
        fail(f"Source universe size differs: {total}")
    eligible = sorted((abs(n - 200000), ty, tx, n)
                      for (tx, ty), n in counts.items()
                      if 100000 <= n <= 300000 and tx * 500 >= minx
                      and (tx + 1) * 500 <= maxx and ty * 500 >= miny
                      and (ty + 1) * 500 <= maxy)
    if len(eligible) != int(expected["roi_eligible_tiles"]) or not eligible:
        fail("Prospective source-only tile eligibility differs")
    _, ty, tx, selected_count = eligible[0]
    if (tx, ty, selected_count) != (int(expected["roi_tile_x"]),
                                    int(expected["roi_tile_y"]),
                                    int(expected["roi_universe_size"])):
        fail("Prospective tile selection differs")
    universe = np.memmap(universe_path, dtype="<u8", mode="r+")
    universe.sort()
    if np.any(universe[1:] <= universe[:-1]):
        fail("Source universe has duplicate original IDs")
    whole_digest = hashlib.sha256(memoryview(universe).cast("B")).hexdigest()
    if whole_digest != expected["source_universe_sorted_u64le_sha256"]:
        fail("Sorted source-universe exact-ID digest differs")
    del universe
    xmin, xmax = tx * 500, (tx + 1) * 500
    ymin, ymax = ty * 500, (ty + 1) * 500
    if (xmin, xmax, ymin, ymax) != tuple(int(expected[key]) for key in
            ("roi_xmin_um", "roi_xmax_um_exclusive", "roi_ymin_um", "roi_ymax_um_exclusive")):
        fail("ROI bounds disagree with validation record")
    roi_ids = []
    assigned_cells = set()
    with (work_dir / "roi_vendor.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("molecule_id", "feature_id", "x_location", "y_location", "z_location", "qv", "vendor_cell_id"))
        for group in range(source.num_row_groups):
            tab = source.read_row_group(group, columns=["transcript_id", "feature_name", "qv", "x_location", "y_location", "z_location", "cell_id"])
            x, y = tab["x_location"].to_numpy(), tab["y_location"].to_numpy()
            eligible_rows = pc.and_(pc.is_in(tab["feature_name"], value_set=feature_set),
                                    pc.greater_equal(tab["qv"], 20)).to_numpy(zero_copy_only=False)
            selected = tab.filter(pa.array(eligible_rows & (x >= xmin) & (x < xmax)
                                           & (y >= ymin) & (y < ymax))).to_pydict()
            for i, molecule_id in enumerate(selected["transcript_id"]):
                cell = selected["cell_id"][i]
                if cell != -1:
                    assigned_cells.add(cell)
                roi_ids.append(molecule_id)
                writer.writerow((molecule_id, selected["feature_name"][i].decode(),
                                 selected["x_location"][i], selected["y_location"][i],
                                 selected["z_location"][i], selected["qv"][i], cell))
    if len(roi_ids) != selected_count or len(set(roi_ids)) != selected_count:
        fail("ROI source IDs/count differ")
    sorted_roi = np.array(sorted(roi_ids), dtype="<u8")
    roi_digest = hashlib.sha256(memoryview(sorted_roi).cast("B")).hexdigest()
    if roi_digest != expected["roi_sorted_u64le_sha256"]:
        fail("ROI exact-ID digest differs")
    cells = pq.read_table(input_dir / "cells.parquet").to_pydict()
    if (len(cells["cell_id"]) != int(expected["vendor_full_census_cells"]) or
            len(set(cells["cell_id"])) != len(cells["cell_id"])):
        fail("Original vendor cell census differs")
    selected_cells = []
    for i, cell in enumerate(cells["cell_id"]):
        if (xmin <= cells["x_centroid"][i] < xmax and
                ymin <= cells["y_centroid"][i] < ymax) or cell in assigned_cells:
            selected_cells.append(i)
    with (work_dir / "roi_vendor_cells.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("cell_id", "x_centroid", "y_centroid"))
        for i in selected_cells:
            writer.writerow((cells["cell_id"][i], cells["x_centroid"][i],
                             cells["y_centroid"][i]))
    if len(selected_cells) != int(expected["vendor_census_cells"]):
        fail("Vendor ROI census differs")
    print(f"Source replay: U={total}, ROI={selected_count}, vendor census={len(selected_cells)}", flush=True)
    return sorted_roi, whole_digest, roi_digest


def baysor_inputs(input_dir, work_dir, expected, sorted_roi):
    universe = np.memmap(work_dir / "universe.ids.u64le", dtype="<u8", mode="r")
    seen = np.zeros(len(universe), dtype=np.bool_)
    with (work_dir / "roi_vendor.tsv").open() as handle:
        vendor = list(csv.DictReader(handle, delimiter="\t"))
    by_id = {int(row["molecule_id"]): row for row in vendor}
    roi_assignments = {}
    assigned = unassigned = native_rows = 0
    columns = pacsv.ConvertOptions(
        include_columns=["transcript_id", "cell_id", "gene", "qv", "cell", "is_noise"],
        column_types={"transcript_id": pa.uint64(), "cell_id": pa.int32(),
                      "gene": pa.string(), "qv": pa.float32(), "cell": pa.string(),
                      "is_noise": pa.bool_()})
    reader = pacsv.open_csv(input_dir / "segmentation.csv",
                            read_options=pacsv.ReadOptions(block_size=8 << 20),
                            convert_options=columns)
    for batch in reader:
        ids = batch.column("transcript_id").to_numpy()
        native_rows += len(ids)
        index = np.searchsorted(universe, ids)
        within = index < len(universe)
        member = np.zeros(len(ids), dtype=np.bool_)
        member[within] = universe[index[within]] == ids[within]
        positions = index[member]
        if np.any(seen[positions]):
            fail("Duplicate saved Baysor source-universe molecule ID")
        seen[positions] = True
        noise = batch.column("is_noise").to_numpy(zero_copy_only=False)[member]
        unassigned += int(noise.sum())
        assigned += len(noise) - int(noise.sum())
        roi_index = np.searchsorted(sorted_roi, ids)
        in_range = roi_index < len(sorted_roi)
        roi_member = np.zeros(len(ids), dtype=np.bool_)
        roi_member[in_range] = sorted_roi[roi_index[in_range]] == ids[in_range]
        if np.any(roi_member):
            selected = batch.filter(pa.array(roi_member)).to_pydict()
            for i, molecule_id in enumerate(selected["transcript_id"]):
                source = by_id[molecule_id]
                if (selected["gene"][i] != source["feature_id"] or
                        selected["cell_id"][i] != int(source["vendor_cell_id"]) or
                        np.float32(selected["qv"][i]) != np.float32(source["qv"])):
                    fail(f"Saved Baysor source fields differ for original ID {molecule_id}")
                if molecule_id in roi_assignments:
                    fail(f"Duplicate saved Baysor ROI molecule ID {molecule_id}")
                cell = selected["cell"][i]
                if selected["is_noise"][i]:
                    if cell not in (None, ""):
                        fail(f"Noise molecule has a forced Baysor cell: {molecule_id}")
                    roi_assignments[molecule_id] = None
                else:
                    if not cell:
                        fail(f"Assigned molecule lacks a Baysor cell: {molecule_id}")
                    roi_assignments[molecule_id] = cell
    if (native_rows != int(expected["baysor_native_rows"]) or
            assigned != int(expected["baysor_full_assigned"]) or
            unassigned != int(expected["baysor_full_unassigned"]) or
            not seen.all() or len(roi_assignments) != len(sorted_roi)):
        fail("Saved Baysor full-universe or ROI accounting differs")
    with (work_dir / "roi_baysor.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("molecule_id", "status", "cell_id"))
        for row in vendor:
            molecule_id = int(row["molecule_id"])
            cell = roi_assignments[molecule_id]
            writer.writerow((molecule_id, "unassigned" if cell is None else "assigned",
                             "" if cell is None else cell))
    assigned_cells = {cell for cell in roi_assignments.values() if cell is not None}
    with (input_dir / "segmentation_cell_stats.csv").open() as handle:
        stats = list(csv.DictReader(handle))
    if (len(stats) != int(expected["baysor_full_census_cells"]) or
            len({row["cell"] for row in stats}) != len(stats)):
        fail("Baysor native census count differs")
    xmin, xmax = int(expected["roi_xmin_um"]), int(expected["roi_xmax_um_exclusive"])
    ymin, ymax = int(expected["roi_ymin_um"]), int(expected["roi_ymax_um_exclusive"])
    selected_cells = [row for row in stats
                      if (xmin <= float(row["x"]) < xmax and
                          ymin <= float(row["y"]) < ymax) or row["cell"] in assigned_cells]
    with (work_dir / "roi_baysor_cells.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("cell_id", "x", "y"))
        for row in selected_cells:
            writer.writerow((row["cell"], row["x"], row["y"]))
    if len(selected_cells) != int(expected["baysor_census_cells"]):
        fail("Baysor ROI census differs")
    print(f"Baysor replay: native rows={native_rows}, full U={assigned + unassigned}, ROI census={len(selected_cells)}", flush=True)
    return native_rows, assigned, unassigned


def main():
    if len(sys.argv) != 3:
        fail("Usage: python3 prepare_ovrlpy_validation.py INPUT_DIR WORK_DIR")
    input_dir, work_dir = (pathlib.Path(arg).expanduser().resolve() for arg in sys.argv[1:])
    if not input_dir.is_dir():
        fail(f"Missing input directory: {input_dir}")
    package_root = pathlib.Path(__file__).resolve().parents[2]
    if work_dir == package_root or package_root in work_dir.parents:
        fail("WORK_DIR must be outside the package source/installation tree")
    work_dir.mkdir(parents=True, exist_ok=True)
    expected = record_values()
    observed = verify_inputs(input_dir, expected)
    sorted_roi, whole_digest, roi_digest = source_inputs(input_dir, work_dir, expected)
    native_rows, assigned, unassigned = baysor_inputs(input_dir, work_dir, expected, sorted_roi)
    manifest = list(observed.items()) + [("source_universe_size", expected["source_universe_size"]),
                ("source_universe_sorted_u64le_sha256", whole_digest),
                ("roi_universe_size", str(len(sorted_roi))),
                ("roi_sorted_u64le_sha256", roi_digest),
                ("baysor_native_rows", str(native_rows)),
                ("baysor_full_assigned", str(assigned)),
                ("baysor_full_unassigned", str(unassigned))]
    with (work_dir / "input_verification.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("key", "value"))
        writer.writerows(manifest)
    print("Preparation and source checks passed; run reproduce_ovrlpy_validation.R", flush=True)


if __name__ == "__main__":
    main()
