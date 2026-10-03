"""Retained CLASS -> separate native DESI probe; no child launcher.

The reviewed parent driver supplies exclusive output/capture, a single job,
deadline and exact source/tool pins. Calling this module is numerical adapter
work and requires that grant. No CLASS invocation or replacement ruler occurs.
"""
import json
from pathlib import Path

import bao
import run as class_run
import theory
import transport

BUILD = "6d495efb166006c6ce651359a366af5d87686ce516ecd7f0f49eed4e8d8be07e"


class PreparationRefusal(ValueError):
    def __init__(self, cause, prefix):
        super().__init__(str(cause)[:2048])
        self.primary_cause = type(cause).__name__
        self.error = class_run.failure(cause)
        self.earned_prefix = prefix


def prepare(attempt, pins_identity, deadline_check):
    """Form four exact ordered means; retain each earned case on later failure.

    The parent must check deadline_check before entry and throughout its own
    capture. No scientific array is serialized into an input-authority file.
    """
    prefix = {"consumed_sources": [], "completed_predictions": [],
              "native_input": None, "dataset": None}
    try:
        deadline_check()
        PINS = Path(pins_identity["path"])
        class_run.require(PINS.is_absolute() and ".." not in PINS.parts, "absolute source-pins path")
        class_run.require(type(pins_identity) is dict and pins_identity.get("path") == str(PINS) and
                          type(pins_identity.get("bytes")) is int and
                          0 < pins_identity["bytes"] <= 65536 and
                          type(pins_identity.get("sha256")) is str and
                          transport.SHA.fullmatch(pins_identity["sha256"]),
                          "positive exact source-pins authority required")
        raw, pin = class_run.file_bytes(PINS, 65536, pins_identity)
        prefix["consumed_sources"].append(pin)
        pins = json.loads(raw, object_pairs_hook=class_run.pairs)
        class_run.require(pins["schema_version"] == 1, "source-pins version")
        authority = {p["role"]: p for p in pins["source_pins"]}
        class_run.require(len(authority) == len(pins["source_pins"]), "duplicate source role")

        def consume(role):
            deadline_check()
            expected = authority[role]
            raw, actual = class_run.file_bytes(expected["path"], expected["bytes"], expected)
            prefix["consumed_sources"].append({"role": role, **actual})
            return raw, actual

        record_raw, record_pin = consume("class_record")
        record = json.loads(record_raw, object_pairs_hook=class_run.pairs)
        reference_raw, reference_pin = consume("class_reference")
        reference = json.loads(reference_raw, object_pairs_hook=class_run.pairs)
        class_run.require(record["reference"]["sha256"] == reference_pin["sha256"] and
                          len(record["cases"]) == 4 and len(pins["cases"]) == 4,
                          "retained record/reference/case identity")
        _, binary = consume("class_binary")
        input_root = Path(authority["desi_mean"]["path"]).parents[2]
        class_run.require(all(authority[role]["path"] == str(input_root / fixed["relative_path"])
                              for role, fixed in zip(("desi_mean", "desi_covariance"), bao.INPUTS)),
                          "released input root/path relationship")
        prefix["dataset"] = bao.read_data(input_root)
        for index, declared in enumerate(pins["cases"]):
            case = record["cases"][index]
            case_id = declared["id"]
            class_run.require(case_id == bao.CASE_IDS[index] == case["id"] and
                              case["status"] == "completed", "retained case order/status")
            parameters = theory.class_parameters(case["ns"])
            class_run.require(parameters == {**reference["parameters"], "n_s": case["ns"]},
                              "same declared CLASS parameter point")
            ini_raw, ini = consume(case_id + "/input")
            case_dir = Path(ini["path"]).parent
            class_run.require(ini_raw == theory.render_ini(parameters, case_dir / "class"),
                              "exact retained INI bytes differ from physical parameter state")
            expected_argv = [binary["path"], ini["path"]]
            sources = {"input_ini": ini, "binary": binary,
                       "retained_record": record_pin, "scientific_reference": reference_pin}
            if case["precision"]:
                _, precision = consume("class_precision")
                sources["supplied_precision"] = precision
                expected_argv.append(precision["path"])
            class_run.require(case["process"]["argv"] == expected_argv and
                              case["ns"] == reference["cases"][index]["ns"] and
                              type(case["precision"]) is bool and
                              case["precision"] == reference["cases"][index]["precision"],
                              "actual CLASS invocation/precision ancestry")
            background_raw, sources["class_background.dat"] = consume(case_id + "/background")
            _, sources["class_thermodynamics.dat"] = consume(case_id + "/thermodynamics")
            stdout_raw, sources["stdout"] = consume(case_id + "/stdout")
            # Only the full retained emitted table is passed to the adapter.
            table = theory.parse_table(background_raw, max_rows=100000)
            state = {"background": table, "drag": theory.predicted_drag(stdout_raw),
                     "parameters": parameters, "source_identities": sources}
            prediction = bao.predict(state, prefix["dataset"], case_id)
            prefix["completed_predictions"].append(prediction)
            class_run.write_new(Path(attempt) / (case_id + ".prediction.json"),
                                class_run.json_bytes(prediction))
            del table, state, background_raw
        deadline_check()
        raw = transport.make_input(prefix["dataset"], prefix["completed_predictions"])
        prefix["native_input"] = class_run.write_new(Path(attempt) / "native.input", raw)
        return prefix
    except BaseException as exc:
        raise PreparationRefusal(exc, prefix) from exc


def terminal_admission(raw, prepared):
    """Parent retains raw stdout/stderr and revokes acceptance on identity drift."""
    value, errors = None, []
    try:
        value = transport.admit_output(raw, prepared["dataset"],
                                       prepared["completed_predictions"], BUILD)
    except BaseException as exc:
        errors.append({"stage": "native-output-admission", **class_run.failure(exc)})
    try:
        transport.verify_data(prepared["dataset"])
    except BaseException as exc:
        errors.append({"stage": "terminal-DESI-source", **class_run.failure(exc)})
    for pin in prepared["consumed_sources"]:
        try:
            class_run.file_bytes(pin["path"], pin["bytes"], pin)
        except BaseException as exc:
            errors.append({"stage": "terminal-CLASS-source", "path": pin["path"],
                           **class_run.failure(exc)})
    if errors:
        exc = PreparationRefusal(ValueError("terminal source identity failure"),
                                 {"native_output": value, "prepared": prepared})
        exc.error["independent_source_failures"] = errors
        raise exc
    return value
