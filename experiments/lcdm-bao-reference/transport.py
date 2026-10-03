"""Experiment-specific native BAO transport/admission draft; no process launch.

Intended for the reviewed bao adapter after SDK/source/artifact admission. The
parent owns compilation, bounded fresh child capture, exact consumed raw input
and output receipts, and source/runtime checks. No command comes from this input.
"""
import json
import math
import re

import bao

INTERFACE = "external-class-bao-gaussian-native/v1"
MAX_INPUT = 16384
MAX_OUTPUT = 32768
POLICY = {"maximum_forward_sensitivity": 1e-8, "maximum_matrix_elements": 169,
          "maximum_models": 4, "maximum_input_bytes": MAX_INPUT,
          "maximum_library_payload_bytes": 1048576}
SHA = re.compile(r"[0-9a-f]{64}\Z")
DENSITY_STATUS = {"finite", "outside_support", "invalid_input", "unsupported_domain",
                  "numerical_failure", "incompatible_metadata"}
NUMERICAL_STATUS = {"ok", "invalid_input", "nonfinite_input", "overflow", "work_limit",
                    "outside_domain", "singular", "not_positive_definite",
                    "conditioning_budget_exceeded"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_data(data):
    """Re-read the original two buffers and compare the actual parsed objects.

    This is called before input formation and again by the outer caller after
    the numerical child. A supplied identity alone cannot bless changed values.
    """
    sources = data.get("source_identities")
    require(type(sources) is list and len(sources) == 2, "missing exact BAO input authority")
    raw = []
    for source, pin in zip(sources, bao.INPUTS):
        require(source.get("role") == pin["role"] and source.get("bytes") == pin["bytes"] and
                source.get("sha256") == pin["sha256"], "BAO source identity differs")
        value, current = bao.read_pinned(source["path"], pin)
        require(all(current[k] == source[k] for k in current), "BAO consumed source identity drift")
        raw.append(value)
    actual = bao.parse_data(*raw)
    require(actual["schema"] == data["schema"] and actual["rows"] == data["rows"] and
            actual["covariance"] == data["covariance"], "BAO parsed object differs from consumed source")


def make_input(data, predictions):
    """Whole 13-vector batch; covariance is consumed once by the native owner."""
    verify_data(data)
    require(type(predictions) is list and 1 <= len(predictions) <= 4, "BAO case batch size")
    lines = ["CLASS_BAO_FIXED_MEAN_V1", f"13 {len(predictions)}",
             bao.INPUTS[0]["sha256"] + " " + bao.INPUTS[1]["sha256"], " ".join(bao.IDS),
             " ".join(repr(bao.finite(x["observed"], "BAO observation")) for x in data["rows"])]
    lines.append(" ".join(repr(bao.finite(x, "BAO covariance")) for row in data["covariance"] for x in row))
    for index, prediction in enumerate(predictions):
        require(prediction["schema"] == "desi-dr2-class-predictions/v1" and
                prediction["case_id"] == bao.CASE_IDS[index], "CLASS case order differs")
        require(type(prediction["common_state_sha256"]) is str and
                SHA.fullmatch(prediction["common_state_sha256"]), "missing CLASS state identity")
        require(prediction["class_source_identities"] is not None, "missing CLASS consumed output identities")
        rows = prediction["rows"]
        require(len(rows) == 13, "CLASS BAO mean dimension")
        values = []
        for actual, row in zip(data["rows"], rows):
            require(all(row[k] == actual[k] for k in actual), "CLASS BAO observation/order identity")
            value = bao.finite(row["predicted"], "CLASS BAO prediction")
            require(value > 0 and row["residual"] == actual["observed"] - value,
                    "CLASS BAO prediction/residual domain")
            values.append(value)
        lines.extend([prediction["case_id"] + " " + prediction["common_state_sha256"],
                      " ".join(repr(x) for x in values)])
    raw = ("\n".join(lines) + "\n").encode("ascii")
    require(len(raw) <= MAX_INPUT, "native BAO input byte bound")
    return raw


def pairs(items):
    result = {}
    for key, value in items:
        require(key not in result, "duplicate native BAO JSON key")
        result[key] = value
    return result


def closed(value, keys, label):
    require(type(value) is dict and set(value) == set(keys), "closed " + label)


def integer(value, lower, upper, label):
    require(type(value) is int and lower <= value <= upper, "integer domain " + label)
    return value


def nullable_status(value, choices, label):
    require(value is None or type(value) is str and value in choices, "status domain " + label)


def close_components(actual, expected, components):
    # Internal identity check only: allow rounding of exposed binary64 terms.
    # This is not a likelihood error allocation, refinement or physics gate.
    scale = max(1.0, *(abs(x) for x in components))
    return abs(actual - expected) <= 64 * math.ulp(scale)


def admit_output(raw, data, predictions, sdk_build_id):
    """Validate earned groups, retaining complete mixed/failed parsed objects.

    The returned object can be retained even if status=='refused'; it must not
    qualify a required batch. Raw capture is an independent caller obligation.
    """
    require(type(sdk_build_id) is str and SHA.fullmatch(sdk_build_id), "positive actual SDK build identity")
    require(type(raw) is bytes and len(raw) <= MAX_OUTPUT, "native BAO output bytes")
    value = json.loads(raw.decode("ascii"), object_pairs_hook=pairs,
                       parse_constant=lambda _: (_ for _ in ()).throw(ValueError("nonfinite native JSON")))
    closed(value, {"schema_version", "interface_id", "status", "stage", "sdk_build_id",
                   "preparation_status", "preparation_numerical_status", "arithmetic_id",
                   "producer_policy", "ordered_ids", "observed", "work", "rows", "error"}, "native output")
    require(type(value["schema_version"]) is int and value["schema_version"] == 1 and
            value["interface_id"] == INTERFACE and value["sdk_build_id"] == sdk_build_id,
            "native BAO interface/SDK identity")
    require(value["status"] in ("accepted", "refused") and
            value["stage"] in ("transport", "preparation", "evaluation", "complete"), "native disposition")
    require(value["ordered_ids"] == list(bao.IDS), "native BAO exact coordinate order")
    closed(value["producer_policy"], POLICY, "native producer policy")
    for key, expected in POLICY.items():
        require(type(value["producer_policy"][key]) is type(expected) and
                value["producer_policy"][key] == expected, "native policy value/type")
    nullable_status(value["preparation_status"], DENSITY_STATUS, "preparation")
    nullable_status(value["preparation_numerical_status"], NUMERICAL_STATUS, "preparation numerical")
    require((value["preparation_status"] is None) == (value["preparation_numerical_status"] is None),
            "paired preparation statuses")
    require(value["arithmetic_id"] is None or value["arithmetic_id"] == "F02/binary64-legacy/v1",
            "native arithmetic identity")
    require((value["arithmetic_id"] is None) == (value["preparation_status"] is None),
            "arithmetic/preparation availability identity")
    closed(value["work"], {"preparation_attempts", "evaluation_attempts"}, "native work")
    preparations = integer(value["work"]["preparation_attempts"], 0, 1, "preparation attempts")
    evaluations = integer(value["work"]["evaluation_attempts"], 0, len(predictions), "evaluation attempts")
    require(value["preparation_status"] is None or preparations == 1, "preparation causal work")
    if value["observed"] is not None:
        require(type(value["observed"]) is list and len(value["observed"]) == 13, "observed echo dimension")
        require([bao.finite(x, "native observation") for x in value["observed"]] ==
                [x["observed"] for x in data["rows"]], "native observed echo")
    require(type(value["rows"]) is list and len(value["rows"]) <= len(predictions), "native row prefix bound")
    attempted, stopped = 0, False
    for index, row in enumerate(value["rows"]):
        closed(row, {"index", "id", "common_state_sha256", "input_complete", "attempted", "accepted",
                     "density_status", "numerical_status", "payload"}, "native row")
        require(type(row["index"]) is int and row["index"] == index and row["id"] == predictions[index]["case_id"]
                and row["common_state_sha256"] == predictions[index]["common_state_sha256"], "native model order/identity")
        require(all(type(row[k]) is bool for k in ("input_complete", "attempted", "accepted")), "native row Boolean types")
        nullable_status(row["density_status"], DENSITY_STATUS, "density")
        nullable_status(row["numerical_status"], NUMERICAL_STATUS, "density numerical")
        require((row["density_status"] is None) == (row["numerical_status"] is None), "paired row statuses")
        if row["attempted"]:
            require(not stopped and row["input_complete"] and preparations == 1 and value["observed"] is not None
                    and value["preparation_status"] == "finite" and value["preparation_numerical_status"] == "ok",
                    "native evaluation causal prefix")
            attempted += 1
        else:
            stopped = True
            require(not row["accepted"] and row["density_status"] is None and row["payload"] is None,
                    "unattempted native group unavailable")
        if row["density_status"] is not None:
            require(row["attempted"], "native status before actual evaluation")
        elif row["attempted"]:
            require(value["error"] is not None and value["stage"] == "evaluation", "failed call has explicit cause")
        payload = row["payload"]
        if row["accepted"]:
            require(row["attempted"] and row["input_complete"] and row["density_status"] == "finite" and
                    row["numerical_status"] == "ok", "accepted native numerical status")
            closed(payload, {"predictions", "residuals", "quadratic", "log_determinant", "normalization",
                             "log_density", "backward_residual", "estimated_forward_sensitivity"}, "native earned density")
            for field in ("predictions", "residuals"):
                require(type(payload[field]) is list and len(payload[field]) == 13, "native vector dimension")
                require([bao.finite(x, "native " + field) for x in payload[field]] ==
                        [x["predicted" if field == "predictions" else "residual"] for x in predictions[index]["rows"]],
                        "native vector echo differs")
            q, ld, norm, logp, backward, forward = [bao.finite(payload[k], "native " + k) for k in
                ("quadratic", "log_determinant", "normalization", "log_density", "backward_residual", "estimated_forward_sensitivity")]
            require(q >= 0 and backward >= 0 and 0 <= forward <= POLICY["maximum_forward_sensitivity"],
                    "native density/diagnostic domain")
            require(close_components(norm, 13 * math.log(2 * math.pi), [norm]) and
                    close_components(logp, -0.5 * (q + ld + norm), [q, ld, norm, logp]),
                    "native density component identity")
        else:
            require(payload is None, "failed native group must not expose default scalar values")
    require(attempted == evaluations, "native evaluation work identity")
    if value["error"] is not None:
        closed(value["error"], {"kind", "stage"}, "native error")
        require(value["error"]["kind"] == "native_refusal" and value["error"]["stage"] == value["stage"], "native error disposition")
    if value["status"] == "accepted":
        require(value["stage"] == "complete" and value["error"] is None and
                value["preparation_status"] == "finite" and value["preparation_numerical_status"] == "ok" and
                value["arithmetic_id"] == "F02/binary64-legacy/v1" and value["observed"] is not None and
                len(value["rows"]) == len(predictions) and evaluations == len(predictions) and
                all(x["accepted"] for x in value["rows"]), "complete native batch acceptance")
    else:
        require(value["stage"] != "complete", "refused native batch cannot claim complete stage")
    return value
