"""A stand-in for Ensembl's REST service, on localhost.

The fetch tests elsewhere inject a fake transport, which proves the parsing and
the URL construction. This goes one step further and serves the responses over
HTTP from a real socket, so the whole chain -- argparse, urllib, the FASTA
writer, the design pipeline -- runs without anything mocked out. What remains
untested is only whether the real Ensembl answers in the shape recorded here.

The sequence served is synthetic. There is no real zebrafish sequence anywhere
in this repository, including in the tests.
"""

from __future__ import annotations

import http.server
import json
import random
import threading
import urllib.parse

#: Shaped like a real ``/lookup/symbol`` response for myo7aa, with the fields
#: this tool reads. Coordinates are reduced so a test locus stays small.
LOOKUP = {
    "myo7aa": {
        "id": "ENSDARG00000099871",
        "display_name": "myo7aa",
        "description": "myosin VIIAa [Source:ZFIN;Acc:ZDB-GENE-020709-1]",
        "assembly_name": "GRCz11",
        "seq_region_name": "18",
        "start": 3001,
        "end": 5000,
        "strand": -1,
        "biotype": "protein_coding",
        "object_type": "Gene",
        "Transcript": [
            {"id": "ENSDART00000172076", "is_canonical": 0},
            {"id": "ENSDART00000148103", "is_canonical": 1},
        ],
    },
    "myo7ab": {
        "id": "ENSDARG00000000002",
        "display_name": "myo7ab",
        "description": "myosin VIIAb [Source:ZFIN]",
        "assembly_name": "GRCz11",
        "seq_region_name": "7",
        "start": 3001,
        "end": 5000,
        "strand": 1,
        "biotype": "protein_coding",
        "object_type": "Gene",
        "Transcript": [{"id": "ENSDART00000000002", "is_canonical": 1}],
    },
}


def _sequence(seed: int, length: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _fasta(name: str, seq: str) -> str:
    body = "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60))
    return f">{name}\n{body}\n"


class _Handler(http.server.BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    def log_message(self, *args) -> None:  # noqa: ANN002 - quiet test output
        pass

    def _send(self, body: str, content_type: str) -> None:
        payload = body.encode()
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def _not_found(self, message: str) -> None:
        payload = json.dumps({"error": message}).encode()
        self.send_response(404)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def do_GET(self) -> None:  # noqa: N802 - http.server's spelling
        path = urllib.parse.urlsplit(self.path).path
        if path.startswith("/lookup/symbol/"):
            symbol = path.rsplit("/", 1)[-1]
            record = LOOKUP.get(symbol)
            if record is None:
                self._not_found(f"No valid lookup found for symbol {symbol}")
                return
            self._send(json.dumps(record), "application/json")
            return
        if path.startswith("/sequence/region/"):
            region = urllib.parse.unquote(path.rsplit("/", 1)[-1])
            locus, span = region.split(":")[0], region.split(":")[1]
            start, end = (int(x) for x in span.split(".."))
            seq = _sequence(hash(locus) & 0xFFFF, end - start + 1)
            self._send(_fasta(f"{locus}:{start}:{end}", seq), "text/x-fasta")
            return
        if path.startswith("/sequence/id/"):
            transcript = path.rsplit("/", 1)[-1]
            seq = _sequence(sum(transcript.encode()) % 9999, 1200)
            self._send(_fasta(f"{transcript} cdna", seq), "text/x-fasta")
            return
        self._not_found(f"no handler for {path}")


class StubEnsembl:
    """Serve the stub on an ephemeral port; yields its base URL."""

    def __init__(self) -> None:
        self._server: http.server.ThreadingHTTPServer | None = None
        self._thread: threading.Thread | None = None

    def __enter__(self) -> str:
        self._server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), _Handler)
        self._thread = threading.Thread(target=self._server.serve_forever, daemon=True)
        self._thread.start()
        host, port = self._server.server_address[:2]
        return f"http://{host}:{port}"

    def __exit__(self, *exc) -> None:  # noqa: ANN002
        if self._server is not None:
            self._server.shutdown()
            self._server.server_close()
        if self._thread is not None:
            self._thread.join(timeout=5)
