# Research workflow

This repository is for a long-running mathematical research project. The
current problem, notation, proved results, failed approaches, and open targets
are maintained in [`ledger.md`](ledger.md). The purpose of this README is only
to preserve the research method so that it is not forgotten during a long
sequence of agent runs.

## Research loop

1. Read this README and the latest entries in `ledger.md` before starting a
   research wave.
2. Based on the ledger, think of ten distinct ideas that could solve the
   problem or materially advance it. Check that they do not merely repeat an
   approach the ledger has already falsified.
3. Rank the ideas by promise, testability, and independence from one another.
4. Select three ideas and spawn three subagents in parallel, giving each one a
   concrete mathematical target and the relevant parts of the ledger or
   supporting artifacts.
5. Check the agents approximately every five minutes. Look for concrete
   progress such as a new identity, a reduced lemma, a proof, a counterexample,
   or a useful computation.
6. Stop an agent when it is clearly stuck: it repeatedly retraces the same
   failed argument, ignores a known obstruction, cannot state a testable next
   step, or makes no concrete progress. Difficulty alone is not evidence that
   an agent is stuck.
7. Give an agent at most twenty minutes by default. Let it continue longer only
   when it is clearly promising—for example, when it has produced a new
   checkable result and is close to finishing a well-defined derivation or
   verification. Reassess any extension at the next five-minute check.
8. Independently inspect every result. Try to falsify it, verify all
   assumptions and constants, and distinguish a proof from numerical evidence
   or an unproved proposal.
9. Update `ledger.md` after every wave and after every material proof,
   falsification, correction, or new open target. Preserve failed approaches
   and their precise obstructions so later agents do not unknowingly repeat
   them.
10. Rerank the best surviving ideas, replace failed ideas with informed
    successors, and start another three-agent wave. Continue this cycle until
    the problem is solved or a genuine external blocker is identified.

## Ledger discipline

Keep the ledger as the durable source of truth. Clearly label statements as:

- **Verified** — the proof has been reconstructed and checked.
- **Pending audit** — plausible, but not yet fully checked.
- **Numerical** — supported by computation only.
- **Falsified** — defeated by a precise gap or counterexample.
- **Open target** — a concrete statement whose proof would advance the
  problem.

Move a claim to **Verified** only after checking definitions, assumptions,
normalizations, constants, signs, factors of two, asymptotic uniformity, and
dependency claims. Whenever practical, perform an independent derivation or a
falsifying computation.

## Markdown conventions

Use GitHub-compatible Markdown throughout the active files:

- use GitHub's backtick-protected inline-math form when an expression contains
  Markdown-sensitive TeX characters such as underscores, asterisks, escaped
  braces, spacing commands, or norm bars;
- use fenced `math` blocks for display mathematics;
- inside a list item, keep short formulas inline and use a protected
  `\displaystyle` expression when a separate indented math fence would be
  interpreted as a code block;
- keep blank lines around display equations, headings, tables, and lists;
- verify that equations, links, tables, and lists render correctly before a
  checkpoint is committed.

## Git checkpoints

Commit and push regularly so the remote repository preserves the research
history. A good default is one coherent commit after each completed agent wave,
with immediate checkpoints for decisive proofs, counterexamples, or
corrections. Review the diff, update the ledger, and reread this README before
committing. Commit messages must accurately distinguish verified mathematics
from pending or numerical work.

Do not begin a new wave from memory alone: return to this README and the ledger
first.

## Python environment

Use the repository-local virtual environment for every Python computation:

```bash
source .venv/bin/activate
```

The direct dependencies and tested versions are recorded in
[`requirements.txt`](requirements.txt). To recreate the environment:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip setuptools wheel
.venv/bin/python -m pip install -r requirements.txt
```

Run artifact scripts with the activated environment or explicitly with
`.venv/bin/python`. Do not install project dependencies into the system Python.
The `.venv/` directory and generated Python caches are intentionally excluded
from Git.
