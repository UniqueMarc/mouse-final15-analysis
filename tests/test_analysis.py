"""Scientific invariants and strict final-15 input-scope checks."""

from copy import deepcopy
from pathlib import Path
import unittest

import numpy as np

from analyze_mice import AnalysisError, measure, source_path, validate_manifest


def fixture_manifest():
    return [{
        "id": f"{time}_{group}", "time_h": time, "group": group, "label": group,
        "original_filename": f"{time}_{group}.tif", "image": f"images/{time}/{group}.tif",
        "sha256": "0" * 64, "scale_min": 10, "scale_max": 20,
        "rois": [{"display_slot": slot, "source_mouse_position": slot,
                  "box": [(slot-1)*20, 0, slot*20, 20]} for slot in (1, 2, 3)],
    } for time in (2, 4, 8) for group in ("A", "B", "C", "D", "E")]


class PixelAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.image = np.full((1200, 1720, 3), 40, dtype=np.uint8)
        rows = np.arange(1001)
        cls.colors = np.column_stack((np.full(1001, 255), rows % 256, rows // 256)).astype(np.uint8)
        cls.image[100:1101, 1500:1520] = cls.colors[:, None, :]

    def test_grayscale_roi_has_zero_accepted_pixels_and_zero_signal(self):
        value, count, mask, _ = measure(self.image, [0, 0, 20, 20], 10, 20)
        self.assertEqual((value, count), (0, 0))
        self.assertFalse(mask.any())

    def test_scale_endpoints_sum_and_signal_outside_roi_is_excluded(self):
        rgb = self.image.copy()
        rgb[3, 2] = self.colors[0]       # Scale maximum: 20.
        rgb[4, 3] = self.colors[-1]      # Scale minimum: 10.
        rgb[5, 25] = self.colors[0]      # Outside the requested ROI.
        value, count, mask, _ = measure(rgb, [0, 0, 20, 20], 10, 20)
        self.assertEqual((value, count, int(mask.sum())), (30, 2, 2))

    def test_color_outside_the_rgb_distance_limit_is_excluded(self):
        rgb = self.image.copy()
        rgb[3, 2] = [0, 255, 0]
        value, count, _, _ = measure(rgb, [0, 0, 20, 20], 10, 20)
        self.assertEqual((value, count), (0, 0))

    def test_invalid_roi_is_rejected(self):
        for box in ([0, 0, 0, 20], [-1, 0, 20, 20], [20, 0, 1, 20], [0, 0, 1721, 20]):
            with self.subTest(box=box), self.assertRaisesRegex(AnalysisError, "ROI"):
                measure(self.image, box, 10, 20)

    def test_invalid_scales_are_rejected(self):
        for lo, hi in ((-1, 20), (20, 10), (10, 10), (float("nan"), 20), (10, float("inf"))):
            with self.subTest(lo=lo, hi=hi), self.assertRaises(AnalysisError):
                measure(self.image, [0, 0, 20, 20], lo, hi)

    def test_missing_color_bar_is_rejected(self):
        with self.assertRaisesRegex(AnalysisError, "color bar"):
            measure(np.full_like(self.image, 40), [0, 0, 20, 20], 10, 20)


class ScopeTests(unittest.TestCase):
    def test_fifteen_panels_with_three_rois_and_five_groups_are_accepted(self):
        validate_manifest(fixture_manifest())

    def test_more_or_fewer_than_fifteen_panels_are_rejected(self):
        records = fixture_manifest()
        for invalid in (records[:-1], records + [deepcopy(records[0])]):
            with self.assertRaisesRegex(AnalysisError, "exactly 15"):
                validate_manifest(invalid)

    def test_extra_roi_is_rejected(self):
        records = fixture_manifest()
        records[0]["rois"].append({"display_slot": 4, "source_mouse_position": 4, "box": [60, 0, 80, 20]})
        with self.assertRaisesRegex(AnalysisError, "exactly three"):
            validate_manifest(records)

    def test_duplicate_group_and_changed_group_set_are_rejected(self):
        for new_group in ("A", "X"):
            records = fixture_manifest()
            records[1]["group"] = new_group
            with self.subTest(group=new_group), self.assertRaises(AnalysisError):
                validate_manifest(records)

    def test_absolute_paths_and_directory_escape_are_rejected(self):
        root = Path("data").resolve()
        for invalid in ("../other/file.tif", "/outside/file.tif", "C:/private/file.tif", "images\\file.tif"):
            with self.subTest(path=invalid), self.assertRaises(AnalysisError):
                source_path(root, invalid)


if __name__ == "__main__":
    unittest.main()
