#!/usr/bin/env python3
from pathlib import Path
p=Path('source/eom_ev_mod.f90')
text=p.read_text(encoding='utf-8')
old='subroutine eom3_KV_st_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf, &\n       yax,yay,yaz,fst,timesub,k)'
new='subroutine eom3_KV_st_ev(ipoint,yxx,yyy,yzz,yst,yvx,yvy,yvz,yvl,yve,ycf, &\n       yax,yay,yaz,fevlocal,fst,timesub,k)'
if text.count(old)!=1:
    raise SystemExit(f'expected one eom3 signature, got {text.count(old)}')
p.write_text(text.replace(old,new,1),encoding='utf-8')
print('fixed eom3_KV_st_ev signature')
