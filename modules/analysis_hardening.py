"""Physics-neutral helpers for robust event identity and worker failures."""


def event_key(source_file_id: int, event_in_file: int) -> tuple[int, int]:
    """Return the authoritative composite event identity."""
    if source_file_id < 0 or event_in_file < 0:
        raise ValueError("event identity components must be non-negative")
    return int(source_file_id), int(event_in_file)


def make_event_id(source_file_id: int, event_in_file: int) -> int:
    """Return a collision-free integer identity for a pair of non-negative IDs.

    This is the Cantor pairing function. It has no events-per-file assumption.
    """
    if source_file_id < 0 or event_in_file < 0:
        raise ValueError("event identity components must be non-negative")
    diagonal = source_file_id + event_in_file
    return diagonal * (diagonal + 1) // 2 + event_in_file


def split_prediction_keys(predictions, file_chunks):
    """Split ``(global_file_id, event_in_file)`` keys into worker-local keys."""
    chunks = []
    file_offset = 0
    for chunk in file_chunks:
        file_limit = file_offset + len(chunk)
        subset = {
            make_event_id(file_id - file_offset, event_in_file): value
            for (file_id, event_in_file), value in predictions.items()
            if file_offset <= file_id < file_limit
        }
        chunks.append(subset)
        file_offset = file_limit
    return chunks


def remap_prediction_keys(predictions, new_position_by_old):
    """Remap tuple prediction keys after file sharding."""
    remapped = {}
    for (old_file_id, event_in_file), value in predictions.items():
        new_file_id = new_position_by_old.get(old_file_id)
        if new_file_id is not None:
            remapped[(new_file_id, event_in_file)] = value
    return remapped


def resolve_worker_result(future, worker_id, logger):
    """Return a worker result or log and fail the complete analysis loudly."""
    try:
        return future.result()
    except Exception as exc:
        logger.exception("Worker %d lanzó excepción; el análisis es inválido", worker_id)
        raise RuntimeError(f"worker {worker_id} failed; analysis aborted") from exc


def has_any_detector_signal(mc_stats) -> bool:
    """Return whether any MC particle has any nonzero detector-signal count."""
    return any(any(count != 0 for count in stats.values()) for stats in mc_stats.values())
