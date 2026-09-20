# External GLMY executable

The local copy of `GLMY.exe` comes from the supplied analysis code. It is excluded
from Git by default. The executable's source code, version, build instructions,
coefficient field, and redistribution license were not supplied. Record these
details before claiming full source-level reproducibility or redistributing it.

`python/homology.py` preserves the original input dialogue (vertices, weighted
directed edges, `#`, dimension parameter `4`, `y`) and expects `homology.json`.
Every invocation records the executable SHA-256 hash. Do not replace it with an
arbitrary executable expecting a different format.
