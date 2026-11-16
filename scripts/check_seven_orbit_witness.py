#!/usr/bin/env python3
"""Replay a seven-orbit mixture using only supplied exact rational scalars."""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import re
import sys


SUPPORT = (92, 3614, 4909, 5596, 6711, 7952, 8163)
CELLS = ("++", "+-", "-+", "--")
ZERO = Fraction(0)
RATIONAL = re.compile(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?\Z")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rational(value, label):
    require(isinstance(value, str) and RATIONAL.fullmatch(value) is not None,
            label + ": expected an integer or numerator/positive-denominator string")
    return Fraction(value)


def sequence(value, length, label):
    require(isinstance(value, list) and len(value) == length,
            label + ": incorrect list length")
    return value


def interval(value, label):
    value = sequence(value, 2, label)
    lower, upper = [rational(v, label) for v in value]
    require(lower <= upper, label + ": reversed interval")
    return [lower, upper]


def intervals(value, label):
    return [interval(v, label + "[" + str(k) + "]")
            for k, v in enumerate(sequence(value, 3, label))]


def exact_tree(value, label):
    """Normalize a nested rational-only saved result, never success flags."""
    if isinstance(value, dict):
        return {k: exact_tree(v, label + "." + k) for k, v in value.items()}
    if isinstance(value, list):
        return [exact_tree(v, label + "[" + str(k) + "]")
                for k, v in enumerate(value)]
    return rational(value, label)


def strings(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {k: strings(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [strings(v) for v in value]
    return value


def containment(predicted, target, label):
    slacks = [predicted[0] - target[0], target[1] - predicted[1]]
    require(all(s >= ZERO for s in slacks),
            label + ": predicted interval is not wholly contained in target")
    return slacks


def check(document):
    require(isinstance(document, dict), "expected a JSON object")
    require(type(document.get("schema_version")) is int
            and document["schema_version"] == 1, "unsupported schema version")
    require(rational(document["unlisted_weight"], "unlisted_weight") == ZERO,
            "all omitted orbit weights must be zero")
    raw_rows = sequence(document["orbits"], 7, "orbits")
    rows = {}
    for item in raw_rows:
        require(isinstance(item, dict), "orbit row must be an object")
        idx = item["id"]
        require(type(idx) is int and idx in SUPPORT and idx not in rows,
                "unexpected or duplicate occupied orbit ID")
        name = "orbit " + str(idx)
        weight = rational(item["weight"], name + " weight")
        require(weight > ZERO, name + ": occupied weight must be positive")
        signs = sequence(item["face_signs"], 2, name + " face_signs")
        require(all(type(v) is int and v in (-1, 1) for v in signs),
                name + ": face signs must be integer -1 or 1")
        require(item["archive_qualification"] == "qualified",
                name + ": archived scalar enclosure is not qualified")
        require(type(item["initial_ground_rank"]) is int
                and item["initial_ground_rank"] == 1,
                name + ": expected archived initial ground rank one")
        bands = intervals(item["band_intervals"], name + " bands")
        down = interval(item["downward_interval"], name + " downward")
        require(down == [ZERO, ZERO],
                name + ": downward interval must be exactly [0,0]")
        rows[idx] = (weight, signs, bands, down)
    require(set(rows) == set(SUPPORT), "incorrect occupied support")
    normalization = sum((rows[idx][0] for idx in SUPPORT), ZERO)
    require(normalization == 1, "weights do not sum to one exactly")

    targets = document["targets"]
    require(isinstance(targets, dict), "targets must be an object")
    target_bands = intervals(targets["bands"], "target bands")
    target_moments = intervals(targets["moments"], "target moments")
    target_cells = targets["joint_populations"]
    require(isinstance(target_cells, dict) and set(target_cells) == set(CELLS),
            "incorrect target joint-population labels")
    target_cells = {k: interval(target_cells[k], "target population " + k)
                    for k in CELLS}

    bands = [[ZERO, ZERO] for _ in range(3)]
    moments = [ZERO, ZERO, ZERO]
    cells = {k: ZERO for k in CELLS}
    downward = [ZERO, ZERO]
    for idx in SUPPORT:
        weight, (s, t), orbit_bands, orbit_down = rows[idx]
        label = ("+" if s == 1 else "-") + ("+" if t == 1 else "-")
        cells[label] += weight
        for k, value in enumerate((s, t, s * t)):
            moments[k] += weight * value
        for endpoint in range(2):
            downward[endpoint] += weight * orbit_down[endpoint]
            for k in range(3):
                bands[k][endpoint] += weight * orbit_bands[k][endpoint]
    require(all(v >= ZERO for v in cells.values())
            and sum(cells.values(), ZERO) == 1,
            "joint populations must be nonnegative and normalized")
    recovered_moments = [
        cells["++"] + cells["+-"] - cells["-+"] - cells["--"],
        cells["++"] - cells["+-"] + cells["-+"] - cells["--"],
        cells["++"] - cells["+-"] - cells["-+"] + cells["--"],
    ]
    require(moments == recovered_moments, "moments and joint populations disagree")
    slacks = {
        "bands": [containment(v, target_bands[k], "band " + str(k + 1))
                  for k, v in enumerate(bands)],
        "moments": [containment([v, v], target_moments[k], "moment " + str(k))
                    for k, v in enumerate(moments)],
        "joint_populations": {
            k: containment([cells[k], cells[k]], target_cells[k], "population " + k)
            for k in CELLS
        },
    }
    require(downward == [ZERO, ZERO], "mixture downward mass is not exactly zero")
    actual = {
        "normalization": normalization,
        "bands": bands,
        "moments": moments,
        "joint_populations": cells,
        "downward_interval": downward,
        "containment_slacks": slacks,
    }
    saved = document["saved_witness"]
    require(isinstance(saved, dict), "saved_witness must be an object")
    for key, value in actual.items():
        require(exact_tree(saved[key], "saved_witness." + key) == value,
                "recomputed " + key + " differs from saved witness")
    saved_weights = {}
    for item in sequence(saved["weights"], 7, "saved weights"):
        require(isinstance(item, dict), "saved weight row must be an object")
        idx = item["id"]
        require(type(idx) is int and idx in SUPPORT and idx not in saved_weights,
                "unexpected or duplicate saved weight ID")
        saved_weights[idx] = rational(item["weight"], "saved weight")
    require(saved_weights == {idx: rows[idx][0] for idx in SUPPORT},
            "exported weights differ from saved witness")
    return strings(dict(
        actual,
        status="PASS",
        occupied_orbits=list(SUPPORT),
        weights=[{"id": idx, "weight": rows[idx][0]} for idx in SUPPORT],
        unlisted_weight=ZERO,
        minimum_downward_mass=ZERO,
        checks={
            "exact_support_and_normalization": True,
            "one_weight_vector_for_all_observables": True,
            "complete_interval_containment": True,
            "every_occupied_downward_interval_is_zero": True,
            "saved_witness_exactly_reproduced": True,
        },
        minimum_reason="Nonnegative unbroadened spectral mass and an attained zero witness.",
        scope="Exact replay of supplied scalar enclosures for the specified local probe.",
    ))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=Path(__file__).resolve().parents[1] / "data" / "seven_orbit_witness.json")
    parser.add_argument("--output", type=Path,
                        help="write deterministic JSON here instead of standard output")
    args = parser.parse_args()
    try:
        document = json.loads(args.input.read_text(encoding="utf-8"))
        result = check(document)
        text = json.dumps(result, indent=2, sort_keys=True) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(text, encoding="utf-8")
        else:
            sys.stdout.write(text)
    except (OSError, ValueError, TypeError, KeyError, ZeroDivisionError) as exc:
        # Paths and local environment details are deliberately absent from check reports.
        detail = "file access failed" if isinstance(exc, OSError) else str(exc)
        print("FAIL: " + detail, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
