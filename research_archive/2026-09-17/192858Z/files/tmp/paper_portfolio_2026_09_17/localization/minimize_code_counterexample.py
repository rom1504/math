"""Bounded search for a smaller exact general-code obstruction."""
import importlib.util
import json
from pathlib import Path
import numpy as np

source = Path("computations/paper_localization_2026_09_17_code_noise.py")
spec = importlib.util.spec_from_file_location("code_noise", source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
rng = np.random.default_rng(17092026)
trials = 0
found = None
for dimension in range(6, 16):
    for rank in range(2, min(dimension - 2, 8) + 1):
        for _ in range(2000):
            trials += 1
            columns = rng.integers(1, 1 << rank, dimension - 1, dtype=np.int64)
            columns = np.concatenate((columns, [np.bitwise_xor.reduce(columns)]))
            found = module.examine(columns, rank, 4)
            if found:
                break
        if found:
            break
    if found:
        break
    print("completed dimension", dimension, "trials", trials, flush=True)
result = {"seed": 17092026, "trials": trials, "counterexample": found}
Path("tmp/paper_portfolio_2026_09_17/localization/code_noise_small.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2), flush=True)
