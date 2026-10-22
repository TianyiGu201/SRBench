"""Plot representative predicted and reference spatial expression maps."""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch
from matplotlib.cm import ScalarMappable
import pickle
import gzip
from matplotlib.colors import Normalize, LinearSegmentedColormap
from PIL import Image

pkl_rngs_XFuse = {'Xe_hsa_breast_c': {'pd_row_rng': (110, 803), 'pd_col_rng': (110, 1060)}, 'Xe_hsa_colon_h': {'pd_row_rng': (137, 958), 'pd_col_rng': (138, 1281)}, 'Xe_hsa_heart_h': {'pd_row_rng': (118, 961), 'pd_col_rng': (110, 455)}, 'Xe_hsa_liver_h': {'pd_row_rng': (216, 1171), 'pd_col_rng': (181, 1614)}, 'Xe_hsa_skin_h': {'pd_row_rng': (181, 1297), 'pd_col_rng': (173, 1614)}, 'Xe_hsa_kidney_h': {'pd_row_rng': (135, 497), 'pd_col_rng': (134, 1266)}, 'ViHD_mmu_kidney_h': {'pd_row_rng': (109, 939), 'pd_col_rng': (119, 857), 'gt_row_rng': (4, 834), 'gt_col_rng': (9, 747)}, 'ViHD_hsa_colon_c': {'pd_row_rng': (105, 873), 'pd_col_rng': (126, 958), 'gt_row_rng': (4, 772), 'gt_col_rng': (3, 835)}, 'ViHD_mmu_embryo_h': {'pd_row_rng': (109, 924), 'pd_col_rng': (114, 952), 'gt_row_rng': (0, 815), 'gt_col_rng': (0, 838)}, 'ViHD_mmu_brain_h': {'pd_row_rng': (117, 907), 'pd_col_rng': (113, 740), 'gt_row_rng': (7, 797), 'gt_col_rng': (89, 716)}, 'ViHD_hsa_pancreas_h': {'pd_row_rng': (118, 955), 'pd_col_rng': (119, 955), 'gt_row_rng': (0, 837), 'gt_col_rng': (1, 837)}, 'cos_hsa_prostate_c': {'pd_row_rng': (14, 146), 'pd_col_rng': (14, 146)}, 'cos_hsa_sarcoma_c': {'pd_row_rng': (46, 441), 'pd_col_rng': (46, 441)}}
need_pad = ['Xe_hsa_heart_h', 'Xe_hsa_liver_h', 'Xe_hsa_skin_h', 'Xe_hsa_kidney_h']
pkl_rngs_iSTAR = {'ViHD_mmu_brain_h': {'pd_row_rng': (19, 808), 'pd_col_rng': (17, 643), 'gt_row_rng': (7, 796), 'gt_col_rng': (89, 715)}, 'ViHD_mmu_embryo_h': {'pd_row_rng': (10, 825), 'pd_col_rng': (13, 851), 'gt_row_rng': (0, 815), 'gt_col_rng': (0, 838)}, 'ViHD_mmu_kidney_h': {'pd_row_rng': (9, 847), 'pd_col_rng': (9, 764), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 755)}, 'ViHD_hsa_colon_c': {'pd_row_rng': (0, 776), 'pd_col_rng': (17, 854), 'gt_row_rng': (0, 776), 'gt_col_rng': (0, 837)}, 'ViHD_hsa_pancreas_h': {'pd_row_rng': (15, 853), 'pd_col_rng': (15, 853), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 838)}}
pkl_rngs_scstGCN = {'ViHD_mmu_brain_h': {'pd_row_rng': (19, 806), 'pd_col_rng': (17, 643), 'gt_row_rng': (8, 795), 'gt_col_rng': (89, 715)}, 'ViHD_mmu_embryo_h': {'pd_row_rng': (10, 813), 'pd_col_rng': (13, 841), 'gt_row_rng': (6, 809), 'gt_col_rng': (5, 833)}, 'ViHD_mmu_kidney_h': {'pd_row_rng': (9, 846), 'pd_col_rng': (9, 755), 'gt_row_rng': (0, 837), 'gt_col_rng': (5, 751)}, 'ViHD_hsa_colon_c': {'pd_row_rng': (0, 767), 'pd_col_rng': (17, 851), 'gt_row_rng': (5, 772), 'gt_col_rng': (2, 836)}, 'ViHD_hsa_pancreas_h': {'pd_row_rng': (15, 853), 'pd_col_rng': (15, 853), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 838)}}
pkl_rngs_TESLA = {'ViHD_mmu_brain_h': {'pd_row_rng': (22, 812), 'pd_col_rng': (19, 646), 'gt_row_rng': (7, 797), 'gt_col_rng': (89, 716)}, 'ViHD_mmu_kidney_h': {'pd_row_rng': (11, 849), 'pd_col_rng': (11, 767), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 756)}, 'ViHD_mmu_embryo_h': {'pd_row_rng': (12, 827), 'pd_col_rng': (17, 855), 'gt_row_rng': (0, 815), 'gt_col_rng': (0, 838)}, 'ViHD_hsa_colon_c': {'pd_row_rng': (2, 779), 'pd_col_rng': (19, 857), 'gt_row_rng': (0, 777), 'gt_col_rng': (0, 838)}, 'ViHD_hsa_pancreas_h': {'pd_row_rng': (18, 856), 'pd_col_rng': (18, 856), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 838)}}
pkl_rngs_SuperST_t = {'ViHD_mmu_kidney_h': {'pd_row_rng': (9, 847), 'pd_col_rng': (9, 765), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 756)}, 'ViHD_mmu_embryo_h': {'pd_row_rng': (11, 826), 'pd_col_rng': (14, 852), 'gt_row_rng': (0, 815), 'gt_col_rng': (0, 838)}, 'ViHD_mmu_brain_h': {'pd_row_rng': (8, 798), 'pd_col_rng': (13, 640), 'gt_row_rng': (7, 797), 'gt_col_rng': (89, 716)}, 'ViHD_hsa_colon_c': {'pd_row_rng': (0, 777), 'pd_col_rng': (14, 852), 'gt_row_rng': (0, 777), 'gt_col_rng': (0, 838)}, 'ViHD_hsa_pancreas_h': {'pd_row_rng': (14, 852), 'pd_col_rng': (15, 853), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 838)}}
pkl_rngs_SuperST_t_g = {'ViHD_mmu_kidney_h': {'pd_row_rng': (10, 848), 'pd_col_rng': (11, 767), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 756)}, 'ViHD_mmu_embryo_h': {'pd_row_rng': (12, 827), 'pd_col_rng': (16, 854), 'gt_row_rng': (0, 815), 'gt_col_rng': (0, 838)}, 'ViHD_mmu_brain_h': {'pd_row_rng': (8, 798), 'pd_col_rng': (15, 642), 'gt_row_rng': (7, 797), 'gt_col_rng': (89, 716)}, 'ViHD_hsa_colon_c': {'pd_row_rng': (0, 777), 'pd_col_rng': (13, 851), 'gt_row_rng': (0, 777), 'gt_col_rng': (0, 838)}, 'ViHD_hsa_pancreas_h': {'pd_row_rng': (14, 852), 'pd_col_rng': (15, 853), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 838)}}
pkl_rngs = {'iSTAR': {'ViHD_mmu_brain_h': {'pd_row_rng': (19, 808), 'pd_col_rng': (17, 643), 'gt_row_rng': (7, 796), 'gt_col_rng': (89, 715)}, 'ViHD_mmu_embryo_h': {'pd_row_rng': (10, 825), 'pd_col_rng': (13, 851), 'gt_row_rng': (0, 815), 'gt_col_rng': (0, 838)}, 'ViHD_mmu_kidney_h': {'pd_row_rng': (9, 847), 'pd_col_rng': (9, 764), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 755)}, 'ViHD_hsa_colon_c': {'pd_row_rng': (0, 776), 'pd_col_rng': (17, 854), 'gt_row_rng': (0, 776), 'gt_col_rng': (0, 837)}, 'ViHD_hsa_pancreas_h': {'pd_row_rng': (15, 853), 'pd_col_rng': (15, 853), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 838)}}, 'scstGCN': {'ViHD_mmu_brain_h': {'pd_row_rng': (19, 806), 'pd_col_rng': (17, 643), 'gt_row_rng': (8, 795), 'gt_col_rng': (89, 715)}, 'ViHD_mmu_embryo_h': {'pd_row_rng': (10, 813), 'pd_col_rng': (13, 841), 'gt_row_rng': (6, 809), 'gt_col_rng': (5, 833)}, 'ViHD_mmu_kidney_h': {'pd_row_rng': (9, 846), 'pd_col_rng': (9, 755), 'gt_row_rng': (0, 837), 'gt_col_rng': (5, 751)}, 'ViHD_hsa_colon_c': {'pd_row_rng': (0, 767), 'pd_col_rng': (17, 851), 'gt_row_rng': (5, 772), 'gt_col_rng': (2, 836)}, 'ViHD_hsa_pancreas_h': {'pd_row_rng': (15, 853), 'pd_col_rng': (15, 853), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 838)}}, 'TESLA': {'ViHD_mmu_brain_h': {'pd_row_rng': (22, 812), 'pd_col_rng': (19, 646), 'gt_row_rng': (7, 797), 'gt_col_rng': (89, 716)}, 'ViHD_mmu_kidney_h': {'pd_row_rng': (11, 849), 'pd_col_rng': (11, 767), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 756)}, 'ViHD_mmu_embryo_h': {'pd_row_rng': (12, 827), 'pd_col_rng': (17, 855), 'gt_row_rng': (0, 815), 'gt_col_rng': (0, 838)}, 'ViHD_hsa_colon_c': {'pd_row_rng': (2, 779), 'pd_col_rng': (19, 857), 'gt_row_rng': (0, 777), 'gt_col_rng': (0, 838)}, 'ViHD_hsa_pancreas_h': {'pd_row_rng': (18, 856), 'pd_col_rng': (18, 856), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 838)}}, 'XFuse': {'Xe_hsa_breast_c': {'pd_row_rng': (110, 803), 'pd_col_rng': (110, 1060)}, 'Xe_hsa_colon_h': {'pd_row_rng': (137, 958), 'pd_col_rng': (138, 1281)}, 'Xe_hsa_heart_h': {'pd_row_rng': (118, 961), 'pd_col_rng': (110, 455)}, 'Xe_hsa_liver_h': {'pd_row_rng': (216, 1171), 'pd_col_rng': (181, 1614)}, 'Xe_hsa_skin_h': {'pd_row_rng': (181, 1297), 'pd_col_rng': (173, 1614)}, 'Xe_hsa_kidney_h': {'pd_row_rng': (135, 497), 'pd_col_rng': (134, 1266)}, 'ViHD_mmu_kidney_h': {'pd_row_rng': (109, 939), 'pd_col_rng': (119, 857), 'gt_row_rng': (4, 834), 'gt_col_rng': (9, 747)}, 'ViHD_hsa_colon_c': {'pd_row_rng': (105, 873), 'pd_col_rng': (126, 958), 'gt_row_rng': (4, 772), 'gt_col_rng': (3, 835)}, 'ViHD_mmu_embryo_h': {'pd_row_rng': (109, 924), 'pd_col_rng': (114, 952), 'gt_row_rng': (0, 815), 'gt_col_rng': (0, 838)}, 'ViHD_mmu_brain_h': {'pd_row_rng': (117, 907), 'pd_col_rng': (113, 740), 'gt_row_rng': (7, 797), 'gt_col_rng': (89, 716)}, 'ViHD_hsa_pancreas_h': {'pd_row_rng': (118, 955), 'pd_col_rng': (119, 955), 'gt_row_rng': (0, 837), 'gt_col_rng': (1, 837)}}, 'SuperST_t': {'ViHD_mmu_kidney_h': {'pd_row_rng': (9, 847), 'pd_col_rng': (9, 765), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 756)}, 'ViHD_mmu_embryo_h': {'pd_row_rng': (11, 826), 'pd_col_rng': (14, 852), 'gt_row_rng': (0, 815), 'gt_col_rng': (0, 838)}, 'ViHD_mmu_brain_h': {'pd_row_rng': (8, 798), 'pd_col_rng': (13, 640), 'gt_row_rng': (7, 797), 'gt_col_rng': (89, 716)}, 'ViHD_hsa_colon_c': {'pd_row_rng': (0, 777), 'pd_col_rng': (14, 852), 'gt_row_rng': (0, 777), 'gt_col_rng': (0, 838)}, 'ViHD_hsa_pancreas_h': {'pd_row_rng': (14, 852), 'pd_col_rng': (15, 853), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 838)}}, 'SuperST_t_g': {'ViHD_mmu_kidney_h': {'pd_row_rng': (10, 848), 'pd_col_rng': (11, 767), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 756)}, 'ViHD_mmu_embryo_h': {'pd_row_rng': (12, 827), 'pd_col_rng': (16, 854), 'gt_row_rng': (0, 815), 'gt_col_rng': (0, 838)}, 'ViHD_mmu_brain_h': {'pd_row_rng': (8, 798), 'pd_col_rng': (15, 642), 'gt_row_rng': (7, 797), 'gt_col_rng': (89, 716)}, 'ViHD_hsa_colon_c': {'pd_row_rng': (0, 777), 'pd_col_rng': (13, 851), 'gt_row_rng': (0, 777), 'gt_col_rng': (0, 838)}, 'ViHD_hsa_pancreas_h': {'pd_row_rng': (14, 852), 'pd_col_rng': (15, 853), 'gt_row_rng': (0, 838), 'gt_col_rng': (0, 838)}}}

class NumpyCompatUnpickler(pickle.Unpickler):

    def find_class(self, module, name):
        if module.startswith('numpy._core'):
            module = module.replace('numpy._core', 'numpy.core')
        return super().find_class(module, name)

def load_exp_matrix(file_path: str):
    if file_path.endswith('.pkl.gz'):
        with gzip.open(file_path, 'rb') as f:
            return NumpyCompatUnpickler(f).load()
    else:
        with open(file_path, 'rb') as f:
            return NumpyCompatUnpickler(f).load()

def load_he_hw(he_path: str):
    with Image.open(he_path) as im:
        w, h = im.size
    return (h, w)

def build_gt_hist_by_he_grid(coords_yx: np.ndarray, pd_h: int, pd_w: int, he_h: int, he_w: int):
    """
    coords_yx: (N,2) float, columns [y, x] in HE pixel coordinate system
    returns hist: (pd_h, pd_w) counts (zeros where uncovered)
    """
    hist = np.zeros((pd_h, pd_w), dtype=np.float32)
    if coords_yx is None or coords_yx.size == 0:
        return hist
    patch_h = he_h / float(pd_h)
    patch_w = he_w / float(pd_w)
    y = coords_yx[:, 0]
    x = coords_yx[:, 1]
    yi = np.floor(y / patch_h).astype(np.int64)
    xi = np.floor(x / patch_w).astype(np.int64)
    yi = np.clip(yi, 0, pd_h - 1)
    xi = np.clip(xi, 0, pd_w - 1)
    np.add.at(hist, (yi, xi), 1.0)
    return hist

def plt_qq(pred, truth, method: str, *, sat: float=0.8, inner: float=0.6, show_legend: bool=True, legend_loc: str='upper left', file_name: str | None=None, s: float=1, alpha: float=1.0, dpi: int=300, min_ticks: int=5):
    """
    QQ scatter: x=Prediction, y=Ground Truth.
    
    Guarantee (AFTER Ctrl+A and checking width in Illustrator):
    - Total object width: 28 mm
    - Axis labels: 7 pt
    - Tick labels: 5 pt
    
    Key principle: Create at larger size, then scale to 28mm in AI.
    Font sizes are pre-compensated so they become correct after scaling.
    """
    MM_PER_INCH = 25.4
    FINAL_WIDTH_MM = 28.0
    ORIGINAL_WIDTH_MM = 140.0
    SCALE = FINAL_WIDTH_MM / ORIGINAL_WIDTH_MM
    LEFT_MARGIN_MM = 31
    RIGHT_MARGIN_MM = 7.5
    TOP_MARGIN_MM = 5
    BOTTOM_MARGIN_MM = 25
    PLOT_WIDTH_MM = ORIGINAL_WIDTH_MM - LEFT_MARGIN_MM - RIGHT_MARGIN_MM
    PLOT_HEIGHT_MM = PLOT_WIDTH_MM
    TOTAL_HEIGHT_MM = PLOT_HEIGHT_MM + TOP_MARGIN_MM + BOTTOM_MARGIN_MM
    TARGET_AXIS_LABEL_PT = 7.0
    TARGET_TICK_PT = 5.0
    AXIS_LABEL_FONTSIZE = TARGET_AXIS_LABEL_PT / SCALE
    TICK_FONTSIZE = TARGET_TICK_PT / SCALE
    LEGEND_FONTSIZE = TARGET_TICK_PT / SCALE
    FIG_WIDTH_IN = ORIGINAL_WIDTH_MM / MM_PER_INCH
    FIG_HEIGHT_IN = TOTAL_HEIGHT_MM / MM_PER_INCH
    fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN), dpi=dpi)
    left_frac = LEFT_MARGIN_MM / ORIGINAL_WIDTH_MM
    bottom_frac = BOTTOM_MARGIN_MM / TOTAL_HEIGHT_MM
    width_frac = PLOT_WIDTH_MM / ORIGINAL_WIDTH_MM
    height_frac = PLOT_HEIGHT_MM / TOTAL_HEIGHT_MM
    ax = fig.add_axes([left_frac, bottom_frac, width_frac, height_frac])
    pred = np.asarray(pred, dtype=float).ravel()
    truth = np.asarray(truth, dtype=float).ravel()
    if pred.shape != truth.shape:
        raise ValueError(f'pred and truth must have the same shape. Got {pred.shape} vs {truth.shape}')
    if sat <= 0:
        raise ValueError('sat must be > 0')
    if inner < 0:
        raise ValueError('inner must be >= 0')
    diff = pred - truth
    color_values = np.zeros_like(diff, dtype=float)
    mask_nonzero = truth > 0
    color_values[mask_nonzero] = diff[mask_nonzero] / truth[mask_nonzero]
    mask_zero = truth == 0
    if np.any(mask_zero):
        color_values[mask_zero & (pred > 0)] = sat
        color_values[mask_zero & (pred < 0)] = -sat
        color_values[mask_zero & (pred == 0)] = 0
    color_values = np.clip(color_values, -sat, sat)
    mid_pos = (inner + sat) / 2.0 / sat
    inner_norm = inner / sat if sat != 0 else 0.0
    colors = ['#0d47a1', '#42a5f5', '#d7191c', '#d7191c', '#d7191c', '#ffeb3b', '#ffa000']
    positions = [0.0, 0.5 - mid_pos / 2, 0.5 - inner_norm / 2, 0.5, 0.5 + inner_norm / 2, 0.5 + mid_pos / 2, 1.0]
    cmap = LinearSegmentedColormap.from_list('custom_diverging', list(zip(positions, colors)), N=256)
    norm = Normalize(vmin=-sat, vmax=sat, clip=True)
    sc = ax.scatter(pred, truth, c=color_values, cmap=cmap, norm=norm, s=s, alpha=alpha, linewidths=0, rasterized=True)
    max_value = float(max(np.max(pred, initial=0.0), np.max(truth, initial=0.0)))
    ax.set_aspect('equal', adjustable='box')
    if max_value > 0:
        num_intervals = max(min_ticks - 1, 4)
        rough_interval = max_value / num_intervals
        nice_interval = max(1, int(np.ceil(rough_interval)))
        tick_positions = np.arange(0, max_value + nice_interval, nice_interval, dtype=int)
        while len(tick_positions) < min_ticks and nice_interval > 1:
            nice_interval -= 1
            tick_positions = np.arange(0, max_value + nice_interval, nice_interval, dtype=int)
        display_max = float(tick_positions[-1])
    else:
        tick_positions = np.array([0], dtype=int)
        display_max = 0.0
    ax.set_xlim(0, display_max)
    ax.set_ylim(0, display_max)
    ax.set_xticks(tick_positions)
    ax.set_yticks(tick_positions)
    ax.plot([0, display_max], [0, display_max], ls='--', c='k', alpha=0.7, linewidth=1)
    ax.set_xlabel('Prediction', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('Ground Truth', fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(axis='both', which='both', labelsize=TICK_FONTSIZE)
    if show_legend:
        handles = [Patch(facecolor='#0d47a1', edgecolor='none', label='Underestimation (>60% below GT)'), Patch(facecolor='#d7191c', edgecolor='none', label='Within 60% of GT'), Patch(facecolor='#ffa000', edgecolor='none', label='Overestimation (>60% above GT)')]
        ax.legend(handles=handles, ncol=3, loc=legend_loc, frameon=False, fontsize=LEGEND_FONTSIZE, handlelength=1.2, handletextpad=0.5, columnspacing=1.2, borderaxespad=0.3)
    if file_name:
        plt.savefig(file_name, format='pdf', dpi=dpi, bbox_inches=None, pad_inches=0)
    plt.show()
    return (fig, ax, sc)

def plot_spatial_gene_prediction(data_path_root: str, method: str, gene_name: str, vmin=0, vmax=20, use_max=False, cmap='RdBu_r', cbar_position='right', fig_name=None, rotate_cw90: bool=False):
    """
    Plot spatial gene expression prediction with fixed physical size.

    Illustrator target (after scaling to 28 mm width):
    - Title font ≈ 6 pt
    - Colorbar tick font ≈ 4 pt
    """
    MM_PER_INCH = 25.4
    FIG_WIDTH_MM = 140
    FINAL_WIDTH_MM = 28
    SCALE = FINAL_WIDTH_MM / FIG_WIDTH_MM
    LEFT_MARGIN_MM = 6
    RIGHT_MARGIN_MM = 6
    TOP_MARGIN_FRAC = 0.1
    CBAR_WIDTH_MM_VERTICAL = 6
    CBAR_HEIGHT_MM_HORIZONTAL = 6
    CBAR_GAP_MM = 2
    TARGET_TITLE_PT = 7
    TARGET_CBAR_TICK_PT = 5
    TITLE_FONTSIZE = TARGET_TITLE_PT / SCALE
    CBAR_TICK_FONTSIZE = TARGET_CBAR_TICK_PT / SCALE
    gene_pred = load_exp_matrix(f'{data_path_root}/{method}/{gene_name}.pkl.gz')
    if rotate_cw90:
        gene_pred = np.rot90(gene_pred, -1)
        cbar_position = 'bottom'
    if use_max:
        vmin = float(np.nanmin(gene_pred))
        vmax = float(np.nanmax(gene_pred))
    norm = Normalize(vmin=vmin, vmax=vmax)
    img_h, img_w = gene_pred.shape
    img_aspect = img_w / img_h
    if cbar_position == 'right':
        available_width_mm = FIG_WIDTH_MM - LEFT_MARGIN_MM - RIGHT_MARGIN_MM - CBAR_WIDTH_MM_VERTICAL - CBAR_GAP_MM
    else:
        available_width_mm = FIG_WIDTH_MM - LEFT_MARGIN_MM - RIGHT_MARGIN_MM
    panel_width_mm = available_width_mm
    panel_width_in = panel_width_mm / MM_PER_INCH
    panel_height_in = panel_width_in / img_aspect
    panel_height_mm = panel_height_in * MM_PER_INCH
    if cbar_position == 'bottom':
        FIG_HEIGHT_MM = panel_height_mm * (1 + TOP_MARGIN_FRAC) + CBAR_HEIGHT_MM_HORIZONTAL + CBAR_GAP_MM
    else:
        FIG_HEIGHT_MM = panel_height_mm * (1 + TOP_MARGIN_FRAC)
    FIG_WIDTH_IN = FIG_WIDTH_MM / MM_PER_INCH
    FIG_HEIGHT_IN = FIG_HEIGHT_MM / MM_PER_INCH
    fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))
    left_frac = LEFT_MARGIN_MM / FIG_WIDTH_MM
    width_frac = panel_width_mm / FIG_WIDTH_MM
    height_frac = panel_height_mm / FIG_HEIGHT_MM
    bottom_frac = 1 - height_frac - TOP_MARGIN_FRAC
    ax = fig.add_axes([left_frac, bottom_frac, width_frac, height_frac])
    ax.imshow(gene_pred, cmap=cmap, norm=norm, origin='upper')
    ax.set_axis_off()
    title_text = f'{gene_name} Ground Truth' if method == 'ground_truth' else f'{method} prediction'
    ax.set_title(title_text, fontsize=TITLE_FONTSIZE, pad=10)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    if cbar_position == 'right':
        cbar_left = (LEFT_MARGIN_MM + panel_width_mm + CBAR_GAP_MM) / FIG_WIDTH_MM
        cbar_width = CBAR_WIDTH_MM_VERTICAL / FIG_WIDTH_MM
        cax = fig.add_axes([cbar_left, bottom_frac, cbar_width, height_frac])
        cbar = fig.colorbar(sm, cax=cax, orientation='vertical')
        ticks = np.linspace(vmin, vmax, 6)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f'{t:.2f}' for t in ticks])
        cbar.ax.tick_params(labelsize=CBAR_TICK_FONTSIZE, rotation=0)
        for label in cbar.ax.get_yticklabels():
            label.set_ha('left')
            label.set_va('center')
    else:
        cbar_height = CBAR_HEIGHT_MM_HORIZONTAL / FIG_HEIGHT_MM
        cbar_bottom = bottom_frac - CBAR_GAP_MM / FIG_HEIGHT_MM - cbar_height
        cax = fig.add_axes([left_frac, cbar_bottom, width_frac, cbar_height])
        cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
        ticks = np.linspace(vmin, vmax, 6)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f'{t:.2f}' for t in ticks])
        cbar.ax.tick_params(labelsize=CBAR_TICK_FONTSIZE, rotation=90)
        for label in cbar.ax.get_xticklabels():
            label.set_ha('center')
            label.set_va('top')
    if fig_name:
        plt.savefig(fig_name, format='pdf', dpi=300)
    plt.show()
    return gene_pred.shape

def dump_pkl_gz(obj, path, compresslevel=6):
    with gzip.open(path, 'wb', compresslevel=compresslevel) as f:
        pickle.dump(obj, f)

def pad_pd_by_ratio(pd, ratio, pad_value=0):
    a, b = pd.shape
    top = int(round(ratio['top'] * a))
    bottom = int(round(ratio['bottom'] * a))
    left = int(round(ratio['left'] * b))
    right = int(round(ratio['right'] * b))
    return np.pad(pd, pad_width=((top, bottom), (left, right)), mode='constant', constant_values=pad_value)

def detect_black_border_ratio(he, eps=0):
    h, w, _ = he.shape
    is_black = np.all(he <= eps, axis=2)
    top = 0
    while top < h and np.all(is_black[top, :]):
        top += 1
    bottom = 0
    while bottom < h and np.all(is_black[h - 1 - bottom, :]):
        bottom += 1
    left = 0
    while left < w and np.all(is_black[:, left]):
        left += 1
    right = 0
    while right < w and np.all(is_black[:, w - 1 - right]):
        right += 1
    return {'top': top / h, 'bottom': bottom / h, 'left': left / w, 'right': right / w}

def gene_grid_from_anndata(gene: str, gt_data, tissue_positions: pd.DataFrame, gt_row_range: tuple, gt_col_range: tuple) -> np.ndarray:
    need_cols = {'barcode', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres'}
    if not need_cols.issubset(tissue_positions.columns):
        raise ValueError(f'tissue_positions must contain columns: {need_cols}')
    if gene not in gt_data.var_names:
        lower_map = {g.lower(): g for g in gt_data.var_names}
        if gene.lower() in lower_map:
            gene = lower_map[gene.lower()]
        else:
            raise KeyError(f"Gene '{gene}' is not in gt_data.var_names.")
    df = tissue_positions.copy()
    df['barcode'] = df['barcode'].astype(str)
    gt_barcodes = pd.Index(gt_data.obs_names.astype(str))
    df = df[df['barcode'].isin(gt_barcodes)].copy()
    if df.empty:
        raise ValueError('No matching barcodes between tissue_positions and AnnData.obs_names.')
    df.set_index('barcode', inplace=True)
    df = df.reindex(gt_barcodes.intersection(df.index))
    X = gt_data[df.index, gene].X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    else:
        X = np.asarray(X)
    if X.ndim == 1:
        expr = X
    elif X.shape[1] > 1:
        expr = X.mean(axis=1)
    else:
        expr = X[:, 0]
    expr = np.asarray(expr).ravel()
    if len(expr) != len(df.index):
        raise ValueError(f'Expression length {len(expr)} does not match df.index length {len(df.index)} for gene {gene}. Check AnnData slicing / duplicated names.')
    expr_s = pd.Series(expr, index=df.index, dtype=float)
    r0, r1 = gt_row_range
    c0, c1 = gt_col_range
    H, W = (r1 - r0, c1 - c0)
    if H <= 0 or W <= 0:
        raise ValueError('gt_row_range / gt_col_range invalid (end must be greater than start).')
    grid = np.full((H, W), 1e-12, dtype=float)
    df_sub = df.loc[(df['array_row'] >= r0) & (df['array_row'] < r1) & (df['array_col'] >= c0) & (df['array_col'] < c1)]
    ar = df_sub['array_row'].astype(float).to_numpy()
    ac = df_sub['array_col'].astype(float).to_numpy()
    yv = df_sub['pxl_row_in_fullres'].astype(float).to_numpy()
    xv = df_sub['pxl_col_in_fullres'].astype(float).to_numpy()

    def _neg_corr(a, b):
        if len(a) < 2 or np.all(a == a[0]) or np.all(b == b[0]):
            return False
        c = np.corrcoef(a, b)[0, 1]
        return np.isfinite(c) and c < 0
    flip_rows = _neg_corr(ar, yv)
    flip_cols = _neg_corr(ac, xv)
    for bc, row in df_sub.iterrows():
        rr = int(row['array_row'])
        cc = int(row['array_col'])
        val = expr_s.get(bc, 1e-12)
        ri = rr - r0
        ci = cc - c0
        if flip_rows:
            ri = H - 1 - ri
        if flip_cols:
            ci = W - 1 - ci
        if 0 <= ri < H and 0 <= ci < W:
            grid[ri, ci] = float(val) if pd.notna(val) else 1e-12
    return grid

