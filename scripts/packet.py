"""Small packet contract shared by validation and execution."""
import hashlib
import json
import math
from pathlib import Path

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def load(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    def number(value):
        result = float(value)
        if not math.isfinite(result):
            raise ValueError(f"Nonfinite JSON number: {value}")
        return result

    return json.loads(Path(path).read_text(), object_pairs_hook=pairs,
                      parse_float=number, parse_constant=number)


def within(root, relative):
    root = Path(root).resolve()
    part = Path(relative)
    if part.is_absolute() or ".." in part.parts or not part.parts:
        raise ValueError(f"Expected a relative path inside {root}: {relative}")
    target = (root / part).resolve()
    if not target.is_relative_to(root):
        raise ValueError(f"Path leaves {root}: {relative}")
    return target


def read_packet(folder):
    folder = Path(folder).resolve()
    schema = load(ROOT / "schemas/experiment.schema.json")
    Draft202012Validator.check_schema(schema)
    packet = load(folder / "experiment.json")
    Draft202012Validator(schema).validate(packet)
    if folder.name != packet["id"] or not (folder / "README.md").is_file():
        raise ValueError("Packet needs a matching directory ID and README.md")
    execution = packet["execution"]
    for item in packet["inputs"]:
        if not item["path"].startswith("data/"):
            raise ValueError("Input paths must live under ignored data/")
        within(ROOT, item["path"])
    origin = packet["origin"]
    design = None
    if origin["kind"] == "prospector_candidate":
        candidate = within(folder, origin["snapshot"])
        if sha256(candidate) != origin["sha256"]:
            raise ValueError("Prospector snapshot hash differs")
        design = load(candidate)
        if design.get("schema_version") != 1 or design.get("kind") != "candidate_design":
            raise ValueError("Not a Prospector candidate-design snapshot")
    if packet["status"] == "blocked":
        return packet, None
    request = within(folder, execution["request"])
    content = load(request)
    if content.get("schema_version") != 2 or content.get("operation") != execution["operation"]:
        raise ValueError("Request schema/operation does not match packet")
    # Irreducible owns the scientific request schema and its runtime validation.
    if design is not None:
        if (design.get("kind") != "candidate_design"
                or design.get("readiness") != "ready_for_consumer_review"
                or design.get("blockers") != [] or design.get("unknowns") != []
                or not isinstance(design.get("review_id"), str) or not design["review_id"]):
            raise ValueError("Candidate is not ready for consumer review")
        consumer = design.get("consumer", {})
        proposed = design.get("executable_request") or {}
        if (consumer.get("repository") != "sprajs/irreducible"
                or consumer.get("revision") != execution["engine_revision"]
                or consumer.get("dirty") is not False
                or proposed.get("sha256") != sha256(request)
                or proposed.get("consumer_validation") != "passed_at_inspected_revision"):
            raise ValueError("Candidate engine/request identity differs; review a new design first")
    return packet, request


def verify_inputs(packet):
    for item in packet["inputs"]:
        path = within(ROOT, item["path"])
        if (not path.is_file() or path.stat().st_size != item["bytes"]
                or sha256(path) != item["sha256"]):
            raise ValueError(f"Input missing or changed: {item['path']}")


def validate_all():
    paths = sorted((ROOT / "experiments").glob("*/experiment.json"))
    if not paths:
        raise ValueError("No experiment packets")
    for path in paths:
        read_packet(path.parent)
    return len(paths)


if __name__ == "__main__":
    print(f"Validated {validate_all()} packets (structure and identity, not science).")
