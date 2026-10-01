"""Quantify the 45 specified ROIs in the 15 final mouse image panels.

Only images named in data/manifest.json are read; folders are never scanned for
additional experiments. The output is recalculated from these final images. It
does not compare against old plotted values, infer animal IDs, or run hypothesis
tests. The display-color algorithm retains the 2026-09-14 numerical operations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np
from PIL import Image, ImageDraw

BASE = Path(__file__).resolve().parent


class AnalysisError(ValueError):
    """Invalid or inconsistent input to the final-15-image analysis."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _finite_number(value, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise AnalysisError(f"{label} must be a finite number; got {value!r}.")
    return value


def validate_measurement(rgb, box, lo, hi) -> None:
    if not isinstance(rgb, np.ndarray) or rgb.ndim != 3 or rgb.shape[2] != 3:
        raise AnalysisError("Source image must be an H x W x 3 RGB array.")
    if rgb.dtype != np.uint8:
        raise AnalysisError("Source image must use 8-bit RGB values (uint8).")
    if not isinstance(box, (list, tuple)) or len(box) != 4:
        raise AnalysisError("ROI box must be [x0, y0, x1, y1].")
    if any(isinstance(v, bool) or not isinstance(v, int) for v in box):
        raise AnalysisError("ROI coordinates must be integers.")
    x0, y0, x1, y1 = box
    height, width = rgb.shape[:2]
    if not (0 <= x0 < x1 <= width and 0 <= y0 < y1 <= height):
        raise AnalysisError(f"ROI {box!r} is empty or outside the {width} x {height} image.")
    _finite_number(lo, "scale_min")
    _finite_number(hi, "scale_max")
    if not 0 <= lo < hi:
        raise AnalysisError("Color scale must satisfy 0 <= scale_min < scale_max.")


def extract_palette(rgb):
    """Locate the source image's color bar and keep every second median row."""
    if rgb.shape[0] < 1160 or rgb.shape[1] < 1700:
        raise AnalysisError("The fixed-layout algorithm requires at least 1700 x 1160 pixels.")
    search = rgb[60:1160, 1450:1700].astype(np.int16)
    saturated = (np.ptp(search, axis=2) > 60) & (search.max(axis=2) > 100)
    counts = saturated.sum(axis=0)
    good = np.flatnonzero(counts > 600)
    if not len(good):
        raise AnalysisError("Could not detect a color bar in the specified search window.")
    runs = np.split(good, np.flatnonzero(np.diff(good) > 1) + 1)
    run = max(runs, key=lambda r: (len(r), counts[r].sum()))
    x0, x1 = int(run[0]) + 1450, int(run[-1]) + 1451
    bar = rgb[60:1160, x0:x1].astype(np.int16)
    rows = np.flatnonzero((np.ptp(bar, axis=2) > 60).sum(axis=1) >= max(3, (x1-x0)//2))
    if not len(rows):
        raise AnalysisError("Detected no usable color-bar rows.")
    y0, y1 = int(rows[0]) + 60, int(rows[-1]) + 61
    palette = np.median(rgb[y0:y1, x0:x1], axis=1)[::2].astype(np.float32)
    if len(palette) < 2:
        raise AnalysisError("At least two color-bar rows are required.")
    return palette, [x0, y0, x1, y1]


def measure(rgb, box, lo, hi):
    """Return signal (a.u.), accepted pixel count/mask, and detected color bar.

    The ROI box [x0, y0, x1, y1] uses exclusive upper bounds. Chromatic pixels
    have max(RGB)-min(RGB)>25 and max(RGB)>55. Each retained color is mapped to
    the closest color-bar entry in Euclidean RGB space, accepting distances <=55.
    A linear high-to-low scale is used and accepted values are summed. The float
    dtypes, block size, first-index tie rule, and sum order are deliberately kept
    as in the original 2026-09-14 algorithm. No background subtraction, area
    normalization, or between-image normalization is performed.
    """
    validate_measurement(rgb, box, lo, hi)
    x0, y0, x1, y1 = box
    region = rgb[y0:y1, x0:x1].astype(np.int16)
    palette, bar = extract_palette(rgb)
    mask = (np.ptp(region, axis=2) > 25) & (region.max(axis=2) > 55)
    pixels = region[mask].astype(np.float32)
    lut = np.linspace(hi, lo, len(palette), dtype=np.float64)
    vals = []
    acc = []
    for start in range(0, len(pixels), 4096):
        block = pixels[start:start+4096]
        d2 = np.sum((block[:, None, :] - palette[None, :, :])**2, axis=2)
        idx = np.argmin(d2, axis=1)
        accepted = np.sqrt(d2[np.arange(len(block)), idx]) <= 55
        vals.extend(lut[idx][accepted])
        acc.extend(accepted)
    accepted_mask = np.zeros(mask.shape, bool)
    accepted_mask[mask] = acc
    return float(np.sum(vals)), int(accepted_mask.sum()), accepted_mask, bar


def validate_manifest(manifest) -> None:
    """Require exactly five groups at 2, 4, and 8 h, three ROIs per panel."""
    if not isinstance(manifest, list) or len(manifest) != 15:
        raise AnalysisError("This analysis requires exactly 15 manifest image panels.")
    required = {"id", "time_h", "group", "label", "original_filename", "image", "sha256", "scale_min", "scale_max", "rois"}
    ids, images = set(), set()
    groups_by_time = {2: set(), 4: set(), 8: set()}
    for index, spec in enumerate(manifest):
        if not isinstance(spec, dict) or not required.issubset(spec):
            raise AnalysisError(f"Manifest panel {index+1} is missing required fields.")
        for field in ("id", "group", "label", "original_filename", "image"):
            if not isinstance(spec[field], str) or not spec[field]:
                raise AnalysisError(f"Panel {index+1}: {field} must be a nonempty string.")
        if spec["id"] in ids or spec["image"] in images:
            raise AnalysisError("Manifest panel IDs and image paths must each be unique.")
        ids.add(spec["id"])
        images.add(spec["image"])
        if isinstance(spec["time_h"], bool) or spec["time_h"] not in (2, 4, 8):
            raise AnalysisError(f"{spec['id']}: time_h must be 2, 4, or 8.")
        groups = groups_by_time[spec["time_h"]]
        if spec["group"] in groups:
            raise AnalysisError(f"{spec['id']}: duplicate group at this time point.")
        groups.add(spec["group"])
        if not isinstance(spec["sha256"], str) or not re.fullmatch(r"[0-9a-f]{64}", spec["sha256"]):
            raise AnalysisError(f"{spec['id']}: sha256 must be 64 lowercase hexadecimal characters.")
        _finite_number(spec["scale_min"], "scale_min")
        _finite_number(spec["scale_max"], "scale_max")
        if not 0 <= spec["scale_min"] < spec["scale_max"]:
            raise AnalysisError(f"{spec['id']}: require 0 <= scale_min < scale_max.")
        rois = spec["rois"]
        if not isinstance(rois, list) or len(rois) != 3:
            raise AnalysisError(f"{spec['id']}: exactly three ROI records are required.")
        slots, positions = set(), set()
        for roi in rois:
            if not isinstance(roi, dict) or not {"display_slot", "source_mouse_position", "box"}.issubset(roi):
                raise AnalysisError(f"{spec['id']}: malformed ROI record.")
            slot, position = roi["display_slot"], roi["source_mouse_position"]
            if isinstance(slot, bool) or not isinstance(slot, int) or slot not in (1, 2, 3):
                raise AnalysisError(f"{spec['id']}: display_slot must be 1, 2, or 3.")
            if isinstance(position, bool) or not isinstance(position, int) or position < 1:
                raise AnalysisError(f"{spec['id']}: source_mouse_position must be a positive integer.")
            if slot in slots or position in positions:
                raise AnalysisError(f"{spec['id']}: ROI slots and source positions must be unique.")
            slots.add(slot)
            positions.add(position)
    if any(len(groups) != 5 for groups in groups_by_time.values()) or not (groups_by_time[2] == groups_by_time[4] == groups_by_time[8]):
        raise AnalysisError("Manifest must contain the same five groups once at each of 2, 4, and 8 h.")


def source_path(data_root: Path, relative: str) -> Path:
    """Resolve a manifest image without permitting an escape from data_root."""
    if Path(relative).is_absolute() or re.match(r"^[A-Za-z]:", relative) or relative.startswith("\\"):
        raise AnalysisError(f"Image path must be relative to the data directory: {relative!r}.")
    if "\\" in relative:
        raise AnalysisError("Manifest image paths must use portable forward slashes.")
    resolved = (data_root / relative).resolve()
    if not resolved.is_relative_to(data_root.resolve()):
        raise AnalysisError(f"Image path escapes the data directory: {relative!r}.")
    return resolved


def write_csv(path: Path, rows) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze(manifest_path: Path, data_root: Path, output_dir: Path, *, write_qc=False) -> dict:
    """Recalculate all 45 ROIs from only the 15 explicitly listed images."""
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise AnalysisError(f"Cannot read manifest {manifest_path}: {exc}") from exc
    validate_manifest(manifest)
    rows, summaries, checksums = [], [], []
    qc_images = []
    for spec in manifest:
        source = source_path(data_root, spec["image"])
        try:
            digest = sha256(source)
            if digest != spec["sha256"]:
                raise AnalysisError(f"Image checksum changed: {source}.")
            with Image.open(source) as image:
                rgb = np.asarray(image.convert("RGB"))
            checksums.append({"image": spec["image"], "sha256": digest, "verified": True})
            panel_values = []
            for roi in sorted(spec["rois"], key=lambda item: item["display_slot"]):
                value, pixels, _, bar = measure(rgb, roi["box"], spec["scale_min"], spec["scale_max"])
                panel_values.append(value)
                x0, y0, x1, y1 = roi["box"]
                rows.append({
                    "panel_id": spec["id"], "time_h": spec["time_h"], "group": spec["group"], "label": spec["label"],
                    "display_slot": roi["display_slot"], "source_mouse_position": roi["source_mouse_position"],
                    "original_filename": spec["original_filename"], "image": spec["image"],
                    "roi_x0": x0, "roi_y0": y0, "roi_x1": x1, "roi_y1": y1,
                    "scale_min": spec["scale_min"], "scale_max": spec["scale_max"],
                    "image_derived_integrated_signal_au": value, "accepted_pixels": pixels,
                })
        except (OSError, ValueError) as exc:
            raise AnalysisError(f"Panel {spec['id']}: {exc}") from exc
        summaries.append({
            "panel_id": spec["id"], "time_h": spec["time_h"], "group": spec["group"], "label": spec["label"],
            "n_rois": len(panel_values),
            "mean_image_derived_integrated_signal_au": float(np.mean(panel_values)),
            "sd_image_derived_integrated_signal_au": float(np.std(panel_values, ddof=1)),
        })
        if write_qc:
            overlay = Image.fromarray(rgb.copy())
            draw = ImageDraw.Draw(overlay)
            for roi in spec["rois"]:
                x0, y0, x1, y1 = roi["box"]
                draw.rectangle((x0, y0, x1-1, y1-1), outline="white", width=3)
                label = f"slot {roi['display_slot']} / source {roi['source_mouse_position']}"
                draw.text((x0+5, y0+5), label, fill="white", stroke_width=1, stroke_fill="black")
            draw.rectangle(tuple(bar), outline="white", width=2)
            name = re.sub(r"[^A-Za-z0-9_.-]", "_", spec["id"]) + ".png"
            qc_images.append((name, overlay))

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "quantification.csv", rows)
    write_csv(output_dir / "summary.csv", summaries)
    if write_qc:
        qc_dir = output_dir / "roi_qc"
        qc_dir.mkdir(exist_ok=True)
        for name, overlay in qc_images:
            overlay.save(qc_dir / name)
    report = {
        "panels_analyzed": len(manifest), "rois_quantified": len(rows),
        "images_checksum_verified": len(checksums),
        "time_points_h": [2, 4, 8], "groups_per_time_point": 5, "rois_per_panel": 3,
        "algorithm": "20260914 nearest-RGB color-bar mapping",
        "signal_unit": "image-derived integrated signal (a.u.)",
        "recalculated_from_final15_images": True,
        "prior_figure_numbers_compared": False,
        "animal_identity_validated": False,
        "hypothesis_tests_performed": False,
        "summary_definition": "Arithmetic mean and sample SD (ddof=1) of three ROIs per panel.",
        "slot_definition": "Display slot and source mouse position are panel positions, not animal IDs.",
        "zero_definition": "No displayed pixels passed the fixed chromatic/color-distance criteria.",
        "qc_images_written": len(qc_images), "source_image_checksums": checksums,
    }
    (output_dir / "verification.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=BASE / "data/manifest.json", help="Manifest containing exactly 15 specified image panels.")
    parser.add_argument("--data-root", type=Path, default=BASE / "data", help="Base directory for manifest image paths; defaults to repository data/.")
    parser.add_argument("--output-dir", type=Path, default=BASE / "outputs", help="CSV/JSON output directory; defaults to repository outputs/.")
    qc = parser.add_mutually_exclusive_group()
    qc.add_argument("--write-qc", action="store_true", help="Explicitly write one ROI overlay per panel (15 additional PNGs in outputs/roi_qc/).")
    qc.add_argument("--skip-qc", action="store_true", help="Write no ROI overlay images (the default).")
    args = parser.parse_args(argv)
    try:
        report = analyze(args.manifest, args.data_root, args.output_dir, write_qc=args.write_qc)
    except (AnalysisError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({key: value for key, value in report.items() if key != "source_image_checksums"}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
