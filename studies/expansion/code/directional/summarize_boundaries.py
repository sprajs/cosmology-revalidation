"""Summarize realized optimizer-boundary occupancy from saved seeded fits."""
import pathlib
import pandas as pd

p=pathlib.Path(__file__).resolve().parents[2]/'runs/directional'
x=pd.read_csv(p/'synthetic-realizations.csv')
out=x.assign(qbound=abs(x.q)>2.999,jbound=abs(x.j)>19.999).groupby(['model','zmax']).agg(n=('q','size'),qbound=('qbound','sum'),jbound=('jbound','sum'),success=('success','sum')).reset_index()
out.to_csv(p/'synthetic-boundaries.csv',index=False)
