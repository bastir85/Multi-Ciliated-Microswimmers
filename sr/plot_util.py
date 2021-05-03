import numpy as np

degree_sign= u'\N{DEGREE SIGN}'

def set_angles_x(ax, angles, ax_in_rad=False):
    if ax_in_rad:
        ax.set_xticks([a/180.*np.pi for a in angles])
    else:
        ax.set_xticks([a for a in angles])
    ax.set_xticklabels([u"{0:.0f}{1}".format(a, degree_sign) for a in angles])

def set_angles_y(ax, angles, ax_in_rad=False):
    if ax_in_rad:
        ax.set_yticks([a/180.*np.pi for a in angles])
    else:
        ax.set_yticks([a for a in angles])
    ax.set_yticklabels([u"{0:.0f}{1}".format(a, degree_sign) for a in angles])

REVT_FULL=7.05826
REVT_HALF=3.40457

NJP_FULL=6.0 #tw is 6.17804in (a bit less)
NJP_SMALL=NJP_FULL*0.75

TEX_TEXTWIDTH = 5.72##6.224 #inch 
ratio = 4/3.
FIG_FULL = (TEX_TEXTWIDTH,TEX_TEXTWIDTH*ratio)
FIG_HALF = (TEX_TEXTWIDTH/2.0,TEX_TEXTWIDTH/2.0)
