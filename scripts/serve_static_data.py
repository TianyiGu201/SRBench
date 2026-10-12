"""Local, read-only object storage simulator with a narrow CORS policy."""

import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path


class DataHandler(SimpleHTTPRequestHandler):
    allowed_origin = ""

    def end_headers(self):
        if self.headers.get("Origin") == self.allowed_origin:
            self.send_header("Access-Control-Allow-Origin", self.allowed_origin)
            self.send_header("Vary", "Origin")
            self.send_header("Access-Control-Allow-Methods", "GET, HEAD, OPTIONS")
            self.send_header("Access-Control-Allow-Headers", "Content-Type")
        super().end_headers()

    def do_OPTIONS(self):
        self.send_response(204)
        self.end_headers()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--port", type=int, default=7801)
    parser.add_argument("--origin", default="http://127.0.0.1:7800")
    args = parser.parse_args()
    DataHandler.allowed_origin = args.origin
    server = ThreadingHTTPServer(("127.0.0.1", args.port), partial(DataHandler, directory=str(args.root.resolve())))
    print(f"Serving {args.root} on http://127.0.0.1:{args.port}; allowed origin {args.origin}", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
