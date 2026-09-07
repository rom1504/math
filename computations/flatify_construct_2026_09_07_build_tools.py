"""Build the finite phase experiment's tools from tracked C++ sources.

Requires g++ with C++17 support. Products are regenerable, not research
inputs. A fresh project-local directory avoids relying on an old binary.
"""
from pathlib import Path
import subprocess
import tempfile


def build_tool(shared=False):
    root = Path(__file__).resolve().parents[1]
    (root / 'tmp').mkdir(exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix='flatify_replay_build_', dir=str(root / 'tmp')))
    source = root / 'computations' / (
        'flatify_construct_2026_09_07_bridge12_eval.cpp' if shared
        else 'exact_fixed_signing_gray.cpp')
    output = directory / ('bridge12_eval.so' if shared else 'exact_gray')
    command = ['g++', '-O3', '-std=c++17']
    if shared:
        command += ['-shared', '-fPIC']
    subprocess.run(command + [str(source), '-o', str(output)], check=True)
    return output
