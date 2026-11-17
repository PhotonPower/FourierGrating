"""``python -m fouriergrating``: prints the compact reference for language-model agents (LLM_REFERENCE.md)."""

import sys


def main() -> int:
    from . import llm_reference

    sys.stdout.write(llm_reference())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
