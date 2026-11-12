"""A minimal HTTP server that honours Range requests.

Python's ``http.server`` ignores ``Range`` and returns the whole file with a
200, which danioprime treats as an error -- correctly, since the point of the
remote-variant feature is to read a few kilobytes of a 214 MB file. Testing
that path therefore needs a server that behaves like the real ones.

Used as a context manager::

    with RangeServer(directory) as base_url:
        ...
"""

from __future__ import annotations

import http.server
import threading
from pathlib import Path


class _Handler(http.server.SimpleHTTPRequestHandler):
    """Serves files, implementing single-range byte requests and HEAD."""

    def log_message(self, *args) -> None:  # noqa: ANN002 - silence the test output
        pass

    def _resolve(self) -> Path | None:
        path = Path(self.translate_path(self.path))
        return path if path.is_file() else None

    def do_HEAD(self) -> None:  # noqa: N802 - http.server naming
        target = self._resolve()
        if target is None:
            self.send_error(404)
            return
        self.send_response(200)
        self.send_header("Content-Length", str(target.stat().st_size))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()

    def do_GET(self) -> None:  # noqa: N802
        target = self._resolve()
        if target is None:
            self.send_error(404)
            return
        size = target.stat().st_size
        header = self.headers.get("Range")
        if not header or not header.startswith("bytes="):
            self.send_response(200)
            self.send_header("Content-Length", str(size))
            self.send_header("Accept-Ranges", "bytes")
            self.end_headers()
            self.wfile.write(target.read_bytes())
            return

        spec = header.split("=", 1)[1].split(",")[0].strip()
        first, _, last = spec.partition("-")
        try:
            start = int(first) if first else 0
            stop = int(last) if last else size - 1
        except ValueError:
            self.send_error(400, "malformed Range")
            return
        if start >= size:
            self.send_response(416)
            self.send_header("Content-Range", f"bytes */{size}")
            self.end_headers()
            return
        stop = min(stop, size - 1)
        with target.open("rb") as fh:
            fh.seek(start)
            payload = fh.read(stop - start + 1)
        self.send_response(206)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Range", f"bytes {start}-{stop}/{size}")
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()
        self.wfile.write(payload)


class RangeServer:
    """Serve *directory* over HTTP with range support, on an ephemeral port."""

    def __init__(self, directory: str | Path):
        self.directory = str(directory)
        self._server: http.server.ThreadingHTTPServer | None = None
        self._thread: threading.Thread | None = None

    def __enter__(self) -> str:
        directory = self.directory

        class Bound(_Handler):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, directory=directory, **kwargs)

        self._server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Bound)
        self._thread = threading.Thread(target=self._server.serve_forever, daemon=True)
        self._thread.start()
        host, port = self._server.server_address[:2]
        return f"http://{host}:{port}"

    def __exit__(self, *exc: object) -> None:
        if self._server is not None:
            self._server.shutdown()
            self._server.server_close()
        if self._thread is not None:
            self._thread.join(timeout=5)
