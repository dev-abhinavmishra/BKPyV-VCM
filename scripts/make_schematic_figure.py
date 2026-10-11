"""Model schematic figure (Figure 1 in CLINICAL_PAPER.md).

Five coupled layers of the 23-dim model — graft, urinary, intracellular,
immune+drugs, NCCR quasi-species — rendered in the paper's ink/paper style.
Output: outputs/figures/model_schematic.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(10.5, 5.6), dpi=170)
ax.set_xlim(0, 10.5); ax.set_ylim(0, 5.6); ax.axis("off")
fig.patch.set_facecolor("#faf8f2")
INK="#1c1a14"; TEAL="#0e5d56"; RUST="#a4442a"; FAINT="#8f897b"

def box(x,y,w,h,fc,ec,label,sub,ls="-"):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.06",
                 fc=fc,ec=ec,lw=1.3,ls=ls))
    ax.text(x+0.14,y+h-0.26,label,fontsize=8.6,family="monospace",
            color=ec,fontweight="bold")
    ax.text(x+0.14,y+h-0.62,sub,fontsize=7.6,color=INK,va="top",linespacing=1.55)

box(0.3,3.4,3.1,1.9,"#0e5d560d",TEAL,"GRAFT (kidney)",
    "C healthy -> I infected -> D dead\nV free virions - L latent reservoir\ncell-to-cell spread (V-independent)")
box(6.9,3.6,3.3,1.5,"#8f897b10",FAINT,"URINARY (Funk topology)",
    "Cu -> Iu -> Vu\ndrained + bidirectional cross-feed",ls=(0,(5,3)))
box(0.3,0.7,3.1,1.9,"#0e5d560d",TEAL,"INSIDE THE INFECTED CELL",
    "S-phase gate -> T-antigen\nDNA synthesis -> capsid\nCC relaxes - TAg never drives cycle")
box(6.9,1.6,3.3,1.7,"#a4442a10",RUST,"IMMUNE + DRUGS",
    "T_naive -> T_eff (primed, capped)\ntac ng/mL -> blunts T-cell control\nsir ng/mL -> removes S-phase perm.")
box(3.75,0.35,3.0,1.15,"#8f897b10",FAINT,"NCCR QUASI-SPECIES",
    "archetype -> rearranged\nkidney pool F_rr - urinary pool F_rr_u",ls=(0,(5,3)))

def arr(x1,y1,x2,y2,c=INK):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",
                 mutation_scale=11,color=c,lw=1.2,shrinkA=2,shrinkB=2))
arr(3.45,4.3,6.85,4.3); ax.text(4.9,4.45,"kidney->urine drainage",fontsize=7,color=FAINT,family="monospace")
arr(1.85,3.35,1.85,2.65); ax.text(2.0,3.0,"infects",fontsize=7,color=FAINT,family="monospace")
arr(6.9,2.5,3.5,4.0); ax.text(5.35,2.95,"immune control + PK effects",fontsize=7,color=FAINT,family="monospace")
arr(5.25,1.55,5.25,2.9); ax.text(5.4,2.15,"rr fraction of both pools",fontsize=7,color=FAINT,family="monospace")
ax.text(0.3,0.12,"23 state variables - ~80 parameters - rates per day - append-only index order",
        fontsize=7.5,color=FAINT,family="monospace")
ax.set_title("BKPyV-VCM - five coupled layers", fontsize=11, color=INK, loc="left", pad=8)
plt.tight_layout()
plt.savefig("outputs/figures/model_schematic.png",facecolor="#faf8f2",bbox_inches="tight")
print("wrote outputs/figures/model_schematic.png")
