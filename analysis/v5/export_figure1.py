"""Export the existing disclosure-controlled Brief Report Figure 1 at 170 mm.

This credential-free route uses the repository's rounded display summaries.
It preserves their numeric precision; it does not refit clinical analyses.
The submission export was replayed from the original full-precision producer.
"""
from pathlib import Path
import argparse,importlib.util,sys,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
ROOT=Path(__file__).resolve().parents[2]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 sys.path.insert(0,str(ROOT/'analysis'));import reproduce_aggregate as base
 spec=importlib.util.spec_from_file_location('v5_figure_gate',ROOT/'tools/small_cell_gate.py');g=importlib.util.module_from_spec(spec);sys.modules[spec.name]=g;spec.loader.exec_module(g)
 hits,_,_=g.scan_paths([ROOT/'outputs/aggregate'],include_restricted=True)
 if hits:raise RuntimeError('Source disclosure gate failed')
 original=plt.Figure.savefig
 def save(fig,path,**kwargs):
  if not getattr(fig,'_v5_size_applied',False):
   scale=(170/25.4)/fig.get_size_inches()[0]
   for t in fig.findobj(Text):t.set_fontsize(t.get_fontsize()*scale)
   fig.set_size_inches(fig.get_size_inches()*scale);fig._v5_size_applied=True
  return original(fig,path,**kwargs)
 plt.Figure.savefig=save
 try:base.plot(a.out)
 finally:plt.Figure.savefig=original
 # The manuscript owns its caption. Exclude the historical producer's extra
 # prose companion from this figure-only export rather than changing formats
 # to evade the data gate.
 caption=a.out/'Figure1_caption.txt'
 if caption.exists():
  caption.unlink()
 hits,_,_=g.scan_paths([a.out],include_restricted=True)
 if hits:raise RuntimeError('Export disclosure gate failed')
 print(json.dumps({'status':'PASS','width_mm':170,'dpi':600,'input_precision':'repository display summaries'}))
if __name__=='__main__':main()
