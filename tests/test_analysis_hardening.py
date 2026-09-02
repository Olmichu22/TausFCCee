import logging
import unittest

from modules.analysis_hardening import (
    make_event_id,
    remap_prediction_keys,
    resolve_worker_result,
    has_any_detector_signal,
    split_prediction_keys,
)


class _FailedFuture:
    def result(self):
        raise ValueError("synthetic worker failure")


class AnalysisHardeningTest(unittest.TestCase):
    def test_boundary_event_ids_do_not_collide(self):
        identities = {
            make_event_id(file_id, event_number)
            for file_id in (0, 1)
            for event_number in (0, 1999)
        }
        self.assertEqual(len(identities), 4)

    def test_997_files_times_2000_events_are_unique_logically(self):
        # Cantor pairing is injective. Check every requested campaign pair by
        # asserting that decoding-independent diagonal intervals never repeat.
        seen = set()
        for file_id in range(997):
            for event_number in range(2000):
                identity = make_event_id(file_id, event_number)
                self.assertNotIn(identity, seen)
                seen.add(identity)
        self.assertEqual(len(seen), 997 * 2000)

    def test_more_than_2000_events_per_file(self):
        identities = {
            make_event_id(file_id, event_number)
            for file_id in range(3)
            for event_number in (0, 1999, 2000, 10_000)
        }
        self.assertEqual(len(identities), 12)

    def test_prediction_keys_are_split_without_event_limit(self):
        predictions = {(0, 5000): "a", (1, 0): "b", (2, 9000): "c"}
        chunks = split_prediction_keys(predictions, [["f0", "f1"], ["f2"]])
        self.assertEqual(chunks[0], {make_event_id(0, 5000): "a", make_event_id(1, 0): "b"})
        self.assertEqual(chunks[1], {make_event_id(0, 9000): "c"})

    def test_prediction_keys_remap_after_sharding(self):
        predictions = {(0, 5000): "a", (1, 0): "b", (2, 9000): "c"}
        self.assertEqual(
            remap_prediction_keys(predictions, {0: 0, 2: 1}),
            {(0, 5000): "a", (1, 9000): "c"},
        )

    def test_detector_signal_checks_every_mc_particle(self):
        stats = {1: {"n_track": 1, "n_ecal": 0}, 2: {"n_track": 0, "n_ecal": 0}}
        self.assertTrue(has_any_detector_signal(stats))
        self.assertFalse(has_any_detector_signal({2: {"n_track": 0, "n_ecal": 0}}))

    def test_worker_exception_is_not_hidden(self):
        with self.assertRaisesRegex(RuntimeError, "worker 7 failed"):
            resolve_worker_result(_FailedFuture(), 7, logging.getLogger("test"))


if __name__ == "__main__":
    unittest.main()
