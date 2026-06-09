#!/usr/bin/env python3
"""Driver: ejecuta todos los escenarios y muestra PASS/FAIL.

Uso:
  python3 run.py              # corre todos los escenarios
  python3 run.py 01 03        # solo escenarios cuyo name empieza por 01 o 03
  python3 run.py --keep       # no borra los archivos intermedios despues
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from framework import run_scenario, parse_tsv
from scenarios import ALL_SCENARIOS


def main(argv: list[str]) -> int:
    filters = [a for a in argv if not a.startswith("--")]
    work = HERE / "work"
    work.mkdir(exist_ok=True)

    scenarios = ALL_SCENARIOS
    if filters:
        scenarios = [s for s in scenarios if any(s.name.startswith(f) for f in filters)]

    width = max(len(s.name) for s in scenarios) + 2
    passed = 0
    failed = 0
    for s in scenarios:
        try:
            ok, errors, out = run_scenario(s, work)
        except Exception as e:  # noqa: BLE001
            print(f"{s.name:<{width}} ERROR  {e}")
            failed += 1
            continue
        if ok:
            print(f"{s.name:<{width}} PASS   {s.description}")
            passed += 1
        else:
            print(f"{s.name:<{width}} FAIL   {s.description}")
            for line in errors:
                print(line)
            # Mostrar todas las filas producidas para ayudar a diagnosticar
            _, rows = parse_tsv(out)
            print(f"  -> {len(rows)} fila(s) producida(s) en {out}")
            for row in rows:
                short = {
                    k: row.get(k, "")
                    for k in (
                        "Positions",
                        "Gene",
                        "Base Changes",
                        "AA Changes",
                        "Variant Type",
                        "Change Type",
                        "Reference Codon",
                        "MNV Codon",
                        "Event Class",
                        "Event Components",
                    )
                }
                print(f"     {short}")
            failed += 1

    print()
    print(f"Total: {passed} PASS, {failed} FAIL")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
