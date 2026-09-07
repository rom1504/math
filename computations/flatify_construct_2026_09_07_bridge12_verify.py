"""Independent full-parent Gray-code replay of the bounded bridge probe."""
import json,subprocess
from pathlib import Path
from flatify_construct_2026_09_07_build_tools import build_tool
root=Path(__file__).resolve().parents[1]
source=root/'computations/results/flatify_construct_2026_09_07_bridge12_probe.json'
record=json.loads(source.read_text()); A=record['best']['parent_matrix']
data=str(len(A))+'\n'+'\n'.join(' '.join(map(str,row)) for row in A)+'\n'
completed=subprocess.run([str(build_tool())],input=data,text=True,capture_output=True,check=True)
answer=json.loads(completed.stdout);assert answer['cap']==record['best']['cap']
answer['source']=str(source.relative_to(root));answer['status']='independent exhaustive full-parent Gray-code replay'
(root/'computations/results/flatify_construct_2026_09_07_bridge12_verify.json').write_text(json.dumps(answer,indent=2)+'\n')
print(json.dumps(answer))
