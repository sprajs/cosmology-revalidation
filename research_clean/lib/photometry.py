"""Strict reader for the released SNANA text photometry tables."""


def read_raisin(path):
    header = {}
    columns = None
    rows = []
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        if not line or line.startswith("#") or ":" not in line:
            continue
        key, value = line.split(":", 1)
        value = value.strip()
        if key == "VARLIST":
            columns = value.split()
        elif key == "OBS":
            assert columns is not None
            values = value.split()
            assert len(values) == len(columns)
            rows.append(
                dict(
                    zip(columns, values),
                    line_number=lineno,
                    observation_index=len(rows),
                )
            )
        elif key not in header:
            header[key] = value
    assert len(rows) == int(header["NOBS"])
    return header, rows
