"""Allow ``python -m danioprime``."""

from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
