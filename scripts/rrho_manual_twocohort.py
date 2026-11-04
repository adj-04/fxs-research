import pandas as pd, numpy as np
from scipy.stats import hypergeom
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R="results/"
def rk(f):
    d=pd.read_csv(f).dropna(subset=["Symbol","pvalue","log2FoldChange"]); d=d[d.Symbol!=""]
    d["m"]=-np.log10(d.pvalue)*np.sign(d.log2FoldChange); d["a"]=d.m.abs()
    return d.sort_values("a",ascending=False).drop_duplicates("Symbol").set_index("Symbol").m
p=rk(R+"DE/DESeq2_FXS_vs_Control_all.csv"); v=rk(R+"validation/GSE218676_validation_DEA_all.csv")
g=p.index.intersection(v.index); n=len(g); print("shared",n)
o1=p[g].sort_values(ascending=False,kind="stable").index; o2=v[g].sort_values(ascending=False,kind="stable").index
r2=pd.Series(np.arange(n),index=o2)[o1].values  # position in list2 of each gene in list1 order
step=n//100; cuts=np.arange(step,n-step+1,step)
M=np.zeros((len(cuts),len(cuts)))
for i,a in enumerate(cuts):
    pos=np.sort(r2[:a])
    for j,b in enumerate(cuts):
        ov=np.searchsorted(pos,b)
        pv=hypergeom.sf(ov-1,n,b,a)
        M[i,j]=-np.log10(max(pv,np.finfo(float).tiny))
print("peak",M.max().round(2))
pd.DataFrame(M,index=cuts,columns=cuts).to_csv("results/validation/RRHO_manual_matrix_twocohort.csv")
fig,ax=plt.subplots(figsize=(7,6))
im=ax.imshow(M.T,origin="lower",cmap="inferno",extent=[cuts[0],cuts[-1],cuts[0],cuts[-1]],aspect="auto")
fig.colorbar(im,ax=ax,label="-log10(p)")
ax.set_xlabel("Organoid rank cutoff (GSE146339)"); ax.set_ylabel("NPC/B-cell/iPSC rank cutoff (GSE218676)")
ax.set_title(f"Manual RRHO cross-check\nPeak -log10(p) = {M.max():.2f} | {n:,} jointly-measured genes",fontsize=11)
fig.tight_layout(); fig.savefig("figures/validation/RRHO_manual_twocohort.png",dpi=300)
