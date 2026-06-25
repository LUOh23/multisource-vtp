import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

fig, ax = plt.subplots(1, 1, figsize=(12, 5.5))
ax.set_xlim(0, 12)
ax.set_ylim(0, 6)
ax.axis('off')
ax.set_facecolor('#FFFFFF')
fig.patch.set_facecolor('#FFFFFF')

# ═══════════ COLORS ═══════════
NAVY      = '#003366'
BLUE      = '#007BFF'
PALE_GRAY = '#F5F7FA'
WHITE     = '#FFFFFF'
TEXT_DARK = '#333333'
TEXT_MID  = '#666666'

# ═══════════ LAYOUT ═══════════
left_x   = 1.9
center_x = 6.0
right_x  = 10.1

mw = 3.0   # module width
mh = 0.65  # module height
gap = 0.12  # gap between modules

def module(cx, y, w, h, title, subtitle=None, bg=PALE_GRAY, tc=TEXT_DARK, border=NAVY, lw=1.5):
    r = FancyBboxPatch((cx-w/2, y-h/2), w, h, boxstyle='round,pad=0.12',
                       facecolor=bg, edgecolor=border, linewidth=lw, zorder=5)
    ax.add_patch(r)
    if subtitle:
        ax.text(cx, y+0.10, title, ha='center', va='center', fontsize=10, fontweight='bold', color=tc, zorder=7)
        ax.text(cx, y-0.14, subtitle, ha='center', va='center', fontsize=8.5, color=TEXT_MID, zorder=7)
    else:
        ax.text(cx, y, title, ha='center', va='center', fontsize=10, fontweight='bold', color=tc, zorder=7)

def dark_module(cx, y, w, h, title, subtitle=None):
    """Navy background, white text."""
    module(cx, y, w, h, title, subtitle, bg=NAVY, tc=WHITE, border=NAVY, lw=1.5)

def arrow(x, y1, y2, color=NAVY, lw=1.5):
    ax.annotate('', xy=(x, y2+mh/2+0.02), xytext=(x, y1-mh/2-0.02),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw), zorder=3)

# ═══════════ LEFT: Traditional ═══════════
ax.text(left_x, 5.35, 'Traditional Single-Source VTP', ha='center', va='center',
        fontsize=12, fontweight='bold', color=NAVY)
ax.text(left_x, 5.05, 'Repeated k times independently', ha='center', va='center',
        fontsize=9, fontstyle='italic', color=TEXT_MID)

y = 4.55
module(left_x, y, mw, mh, 'Input', 'Mesh M, source vertex s')
arrow(left_x, y, y-gap)

y -= gap + mh
dark_module(left_x, y, mw, mh, 'Reset  clear()', 'Reset all vertices & edges')
arrow(left_x, y, y-gap)

y -= gap + mh
module(left_x, y, mw, mh, 'Init  init_one_source(s)', 'distance(s)=0, create windows')
arrow(left_x, y, y-gap)

y -= gap + mh
module(left_x, y, mw, mh, 'Propagate', 'Wavefront: Rule 1 + Rule 2')

y -= gap + mh
module(left_x, y, mw, mh, 'Output', 'd(v) = d_g(v_s, v)')

# k-times annotation
ax.annotate('', xy=(left_x+mw/2+0.45, 4.25), xytext=(left_x+mw/2+0.45, 2.0),
            arrowprops=dict(arrowstyle='<->', color='#CC3333', lw=1.8, ls='--'), zorder=2)
ax.text(left_x+mw/2+0.65, 3.12, 'k', fontsize=11, fontweight='bold', color='#CC3333',
        va='center', rotation=90)
ax.text(left_x+mw/2+0.85, 3.12, 'times', fontsize=9, color='#CC3333', va='center', rotation=90)

# ═══════════ CENTER: Innovation ═══════════
iw, ih = 2.4, 2.8
cy = 3.65
r = FancyBboxPatch((center_x-iw/2, cy-ih/2), iw, ih, boxstyle='round,pad=0.25',
                   facecolor=BLUE, edgecolor=BLUE, linewidth=2.0, zorder=10)
ax.add_patch(r)
ax.text(center_x, cy+1.05, 'KEY INNOVATION', ha='center', va='center',
        fontsize=10, fontweight='bold', color=WHITE, zorder=11)
ax.text(center_x, cy+0.35, 'Single clear() call', ha='center', va='center',
        fontsize=11, fontweight='bold', color=WHITE, zorder=11)
ax.text(center_x, cy-0.15, 'for ALL k sources', ha='center', va='center',
        fontsize=11, fontweight='bold', color=WHITE, zorder=11)
ax.text(center_x, cy-0.85, 'Eliminates k-1 redundant', ha='center', va='center',
        fontsize=9, color='#CCE0FF', zorder=11)
ax.text(center_x, cy-1.15, 'mesh topology rebuilds', ha='center', va='center',
        fontsize=9, color='#CCE0FF', zorder=11)

# ═══════════ RIGHT: Multi-Source (Ours) ═══════════
ax.text(right_x, 5.35, 'Multi-Source VTP  (This Paper)', ha='center', va='center',
        fontsize=12, fontweight='bold', color=NAVY)
ax.text(right_x, 5.05, 'Single pass, k sources merged', ha='center', va='center',
        fontsize=9, fontstyle='italic', color=TEXT_MID)

y = 4.55
module(right_x, y, mw, mh, 'Input', 'Mesh M, source set S = {s1, ..., sk}')
arrow(right_x, y, y-gap)

y -= gap + mh
dark_module(right_x, y, mw, mh, 'Reset  clear()', 'Reset once for all k sources')
# ONLY ONCE annotation
ax.text(right_x+mw/2+0.15, y+0.18, 'ONLY', fontsize=7, fontweight='bold', color='#CC3333', ha='left')
ax.text(right_x+mw/2+0.15, y-0.02, 'ONCE', fontsize=7, fontweight='bold', color='#CC3333', ha='left')
arrow(right_x, y, y-gap)

y -= gap + mh
# Init with blue border accent
r = FancyBboxPatch((right_x-mw/2, y-mh/2), mw, mh, boxstyle='round,pad=0.12',
                   facecolor=WHITE, edgecolor=BLUE, linewidth=2.0, zorder=5)
ax.add_patch(r)
ax.text(right_x, y+0.10, 'Init  init_all_sources(S)', ha='center', va='center',
        fontsize=10, fontweight='bold', color=NAVY, zorder=7)
ax.text(right_x, y-0.14, 'Batch init + Deduplication', ha='center', va='center',
        fontsize=8.5, color=TEXT_MID, zorder=7)
arrow(right_x, y, y-gap)

y -= gap + mh
dark_module(right_x, y, mw, mh, 'Propagate  (unchanged)', 'Wavefront: Rule 1 + Rule 2')
ax.text(right_x, y-mh/2-0.16, 'Mutual suppression occurs here', ha='center', va='top', fontsize=7.5, fontstyle='italic', color='#CC3333', fontfamily='Times New Roman')
arrow(right_x, y, y-gap)

y -= gap + mh
module(right_x, y, mw, mh, 'Output', 'd(v) = min_j d_g(v_sj, v)')

# ═══════════ BOTTOM: Speedup ═══════════
ax.text(6.0, 0.55, r'Acceleration:  $\eta(k) \to k$  (Property 2)',
        ha='center', va='center', fontsize=11, fontweight='bold', color=NAVY,
        fontfamily='Times New Roman')
ax.text(6.0, 0.22, r'$T_{\rm multi} \approx$ const  (independent of $k$)',
        ha='center', va='center', fontsize=9, fontstyle='italic', color=TEXT_MID,
        fontfamily='Times New Roman')

plt.tight_layout(pad=0.3)
plt.savefig('e:/vtp_geodesic/论文截图/fig1_algorithm_flow.png', dpi=300,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('Fig.2 saved')
