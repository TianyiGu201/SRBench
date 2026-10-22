"""Rank FOV predictions and plot spatial overlays."""
from __future__ import annotations
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from .fov_status import mask_non_evaluable_fovs


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'ARIAL'
HIGHER_IS_BETTER = {'PCC', 'SCC', 'SSIM', 'PSNR'}
_RANGE_RE = re.compile('minx(?P<minx>-?\\d+(?:\\.\\d+)?)_maxx(?P<maxx>-?\\d+(?:\\.\\d+)?)_miny(?P<miny>-?\\d+(?:\\.\\d+)?)_maxy(?P<maxy>-?\\d+(?:\\.\\d+)?)')

def parse_range_from_name(path_or_name: str):
    """Extract minx/maxx/miny/maxy bounds from a filename or full path."""
    base = os.path.splitext(os.path.basename(path_or_name))[0]
    m = _RANGE_RE.search(base)
    if m is None:
        return None
    d = {k: float(v) for k, v in m.groupdict().items()}
    return (d['minx'], d['maxx'], d['miny'], d['maxy'])

def read_table_auto(path):
    return pd.read_csv(path, sep='\t' if path.endswith('.tsv') else ',')

def _rank_worst_invalid(series, ascending: bool):
    """Rank numeric values, assigning NaN and infinite values the worst ranks."""
    s = pd.to_numeric(series, errors='coerce').replace([np.inf, -np.inf], np.nan)
    r = s.rank(ascending=ascending, method='average', na_option='bottom')
    if r.notna().any():
        worst = float(r.max())
        r = r.fillna(worst)
    else:
        r = pd.Series(np.ones(len(s), dtype=float), index=s.index)
    return r

def _minmax_0_1(arr: np.ndarray):
    """Scale finite values to [0, 1]; map equal finite values to 0 and invalid values to 1."""
    x = np.asarray(arr, dtype=float)
    finite = np.isfinite(x)
    out = np.empty_like(x, dtype=float)
    out[:] = 1.0
    if not finite.any():
        return out
    xmin = x[finite].min()
    xmax = x[finite].max()
    if xmax == xmin:
        out[finite] = 0.0
        return out
    out[finite] = (x[finite] - xmin) / (xmax - xmin)
    return out

def prepare_rank_par_norm_fovsize_level(fov_dirs: dict, out_csv_path: str, metrics=('SCC', 'PCC', 'KLD', 'JSD'), out_metric_name='PARnorm'):
    """Calculate gene-level PAR and PARnorm within each FOV size.

    Pool gene rows over regions, rank each metric, average metric ranks to PAR,
    then min-max normalize PAR. Retain region coordinates in the combined CSV.
    Whole-FOV failures are masked after ranking.
    """
    higher_better_for_rank = {'SCC', 'PCC'}
    all_sizes_rows = []
    for fov_size, in_dir in fov_dirs.items():
        rows_size = []
        for fn in os.listdir(in_dir):
            if not fn.endswith(('.csv', '.tsv')):
                continue
            in_path = os.path.join(in_dir, fn)
            r = parse_range_from_name(in_path)
            if r is None:
                continue
            minx, maxx, miny, maxy = r
            df = read_table_auto(in_path).copy()
            base = os.path.splitext(os.path.basename(fn))[0]
            df['FOV_size'] = fov_size
            df['FOV_region'] = base
            df['minx'], df['maxx'], df['miny'], df['maxy'] = (minx, maxx, miny, maxy)
            for m in metrics:
                if m not in df.columns:
                    df[m] = np.nan
            rows_size.append(df)
        if not rows_size:
            continue
        df_size = pd.concat(rows_size, ignore_index=True)
        rank_cols = []
        for m in metrics:
            asc = m not in higher_better_for_rank
            col_rank = f'{m}_rank'
            df_size[col_rank] = _rank_worst_invalid(df_size[m], ascending=asc)
            rank_cols.append(col_rank)
        df_size['PAR'] = df_size[rank_cols].mean(axis=1)
        df_size[out_metric_name] = _minmax_0_1(df_size['PAR'].to_numpy())
        all_sizes_rows.append(df_size)
    if not all_sizes_rows:
        raise ValueError('No valid FOV data found. Check fov_dirs paths and filename range pattern.')
    df_all = mask_non_evaluable_fovs(pd.concat(all_sizes_rows, ignore_index=True))
    os.makedirs(os.path.dirname(out_csv_path) or '.', exist_ok=True)
    df_all.to_csv(out_csv_path, index=False)
    return df_all

def build_grid(regions, tol=1e-06):
    ys = sorted({round(r['cy'] / tol) * tol for r in regions})
    grid = []
    for y in ys:
        row = [r for r in regions if abs(r['cy'] - y) < tol]
        row = sorted(row, key=lambda d: d['cx'])
        grid.append(row)
    return grid

def pick_by_row_col(regions, row, col):
    grid = build_grid(regions)
    row = min(max(row, 0), len(grid) - 1)
    col = min(max(col, 0), len(grid[row]) - 1)
    return grid[row][col]

def auto_box_draw_size(regions):
    xs = sorted({r['cx'] for r in regions})
    ys = sorted({r['cy'] for r in regions})
    dx = min(np.diff(xs)) if len(xs) > 1 else np.inf
    dy = min(np.diff(ys)) if len(ys) > 1 else np.inf
    return min(dx, dy)

def find_text_position(hl, he_shape, text_pref=None, pad=2, fontsize=11, dpi=None, char_width_factor=0.6):
    """Position text beside a highlighted FOV, with padding measured in character widths."""
    H, W = he_shape[:2]
    if dpi is None:
        dpi = plt.rcParams.get('figure.dpi', 100)
    pad_px = float(pad) * char_width_factor * float(fontsize) * (float(dpi) / 72.0)
    if text_pref is None:
        candidates = [('right', 'top'), ('right', 'bottom'), ('left', 'top'), ('left', 'bottom'), ('top', 'left'), ('top', 'right'), ('bottom', 'left'), ('bottom', 'right')]
    else:
        candidates = [(text_pref['side'], text_pref['align'])]
    for side, align in candidates:
        if side == 'right':
            x = hl['maxx'] + pad_px
            if align == 'top':
                y = hl['miny']
                va = 'top'
            else:
                y = hl['maxy']
                va = 'bottom'
            ha = 'left'
        elif side == 'left':
            x = hl['minx'] - pad_px
            if align == 'top':
                y = hl['miny']
                va = 'top'
            else:
                y = hl['maxy']
                va = 'bottom'
            ha = 'right'
        elif side == 'top':
            y = hl['miny'] - pad_px
            if align == 'left':
                x = hl['minx']
                ha = 'left'
            else:
                x = hl['maxx']
                ha = 'right'
            va = 'bottom'
        else:
            y = hl['maxy'] + pad_px
            if align == 'left':
                x = hl['minx']
                ha = 'left'
            else:
                x = hl['maxx']
                ha = 'right'
            va = 'top'
        if 0 <= x <= W and 0 <= y <= H:
            return (x, y, ha, va)
    return (hl['minx'], hl['miny'], 'left', 'top')

def collect_regions_from_bigcsv(df_all: pd.DataFrame, fov_size: str, metric: str):
    """Average finite gene-level metric values within each FOV region."""
    df = df_all[df_all['FOV_size'] == fov_size].copy()
    if df.empty:
        return []
    if metric not in df.columns:
        return []
    vals = pd.to_numeric(df[metric], errors='coerce').replace([np.inf, -np.inf], np.nan)
    df['_metric_val_'] = vals
    group_cols = ['FOV_region', 'minx', 'maxx', 'miny', 'maxy']
    g = df.groupby(group_cols, dropna=False)['_metric_val_'].mean()
    regions = []
    for (region, minx, maxx, miny, maxy), v in g.items():
        regions.append(dict(minx=float(minx), maxx=float(maxx), miny=float(miny), maxy=float(maxy), cx=0.5 * (float(minx) + float(maxx)), cy=0.5 * (float(miny) + float(maxy)), value=float(v) if np.isfinite(v) else np.nan, region=region))
    return regions

def draw_fov_panel(ax, he_np, regions, metric, fov_name, highlight_index=(0, 0), box_draw_size=None, text_pref=None, norm=None, cmap=None, rotate_cw90=False, title_fontsize=None, text_fontsize=None):
    """
    Draw a single FOV panel with regions overlaid
    
    Parameters:
    -----------
    title_fontsize, text_fontsize: float
        Font sizes pre-calculated for final output size
    """
    if rotate_cw90:
        H, W = he_np.shape[:2]
        he_np = np.rot90(he_np, k=-1)
        regions = rotate_regions_cw90(regions, H)
    ax.imshow(he_np, origin='upper')
    ax.set_axis_off()
    fov_value = str(fov_name).lower().replace('fov', '')
    ax.set_title(f'FOV size: {fov_value} pixels', fontsize=title_fontsize, pad=4)
    if box_draw_size is None:
        box_draw_size = auto_box_draw_size(regions)
    half = box_draw_size / 2
    for r in regions:
        ax.add_patch(Rectangle((r['cx'] - half, r['cy'] - half), box_draw_size, box_draw_size, facecolor=cmap(norm(r['value'])), edgecolor='none'))
    hl = pick_by_row_col(regions, *highlight_index)
    ax.add_patch(Rectangle((hl['minx'], hl['miny']), hl['maxx'] - hl['minx'], hl['maxy'] - hl['miny'], fill=False, edgecolor='black', linewidth=2))
    tx, ty, ha, va = find_text_position(hl, he_np.shape, text_pref, fontsize=text_fontsize)

def is_region_valid(df_all: pd.DataFrame, fov_size: str, fov_region: str, metrics):
    """A whole FOV is non-evaluable only when all requested raw metrics are undefined."""
    group = df_all[(df_all['FOV_size'] == fov_size) & (df_all['FOV_region'] == fov_region)]
    if group.empty:
        return False
    values = group.reindex(columns=list(metrics)).apply(pd.to_numeric, errors='coerce')
    return bool(values.replace([np.inf, -np.inf], np.nan).notna().any().any())

def rotate_regions_cw90(regions, H):
    """Rotate region bounds clockwise by 90 degrees using the original image height H."""
    out = []
    for r in regions:
        minx, maxx = (r['minx'], r['maxx'])
        miny, maxy = (r['miny'], r['maxy'])
        pts = [(minx, miny), (minx, maxy), (maxx, miny), (maxx, maxy)]
        pts_rot = [(H - y, x) for x, y in pts]
        xs = [p[0] for p in pts_rot]
        ys = [p[1] for p in pts_rot]
        nr = r.copy()
        nr['minx'], nr['maxx'] = (min(xs), max(xs))
        nr['miny'], nr['maxy'] = (min(ys), max(ys))
        nr['cx'] = 0.5 * (nr['minx'] + nr['maxx'])
        nr['cy'] = 0.5 * (nr['miny'] + nr['maxy'])
        out.append(nr)
    return out

def plot_fovs(he_np, big_csv_path: str, metric: str, fov_sizes=None, highlight_index=(0, 0), box_draw_size=None, text_pref=None, vmin=None, vmax=None, fig_name=None, rotate_cw90=False):
    """
    Plot FOV panels with fixed output width and maximized horizontal fill
    
    Key features:
    - Output PDF width is always 140 mm
    - Images fill ~95% of horizontal space
    - Colorbar has fixed physical size (3 mm for vertical, 2.5 mm for horizontal)
    - Font sizes are pre-scaled for 28 mm final width in Illustrator
    """
    MM_PER_INCH = 25.4
    FIG_WIDTH_MM = 140
    FINAL_WIDTH_MM = 28
    SCALE = FINAL_WIDTH_MM / FIG_WIDTH_MM
    LEFT_MARGIN_MM = 5.2
    RIGHT_MARGIN_MM = 5.2
    CBAR_WIDTH_MM_VERTICAL = 6
    CBAR_HEIGHT_MM_HORIZONTAL = 6
    CBAR_GAP_MM = 2
    TARGET_TITLE_PT = 6
    TARGET_CBAR_TICK_PT = 5
    TITLE_FONTSIZE = TARGET_TITLE_PT / SCALE
    TEXT_FONTSIZE = 5 / SCALE
    CBAR_TICK_FONTSIZE = TARGET_CBAR_TICK_PT / SCALE
    FIG_WIDTH_IN = FIG_WIDTH_MM / MM_PER_INCH
    df_all = mask_non_evaluable_fovs(pd.read_csv(big_csv_path))
    if fov_sizes is None:
        fov_sizes = sorted(df_all['FOV_size'].dropna().unique().tolist())
    REQUIRED_METRICS = ('PCC', 'SCC', 'JSD', 'KLD')
    higher_is_better = metric in HIGHER_IS_BETTER
    cmap = plt.get_cmap('RdBu_r' if higher_is_better else 'RdBu')
    regions_all = {}
    values_all = []
    for fsz in fov_sizes:
        df_fsz = df_all[df_all['FOV_size'] == fsz]
        if df_fsz.empty:
            continue
        valid_regions = []
        for region_name in df_fsz['FOV_region'].unique():
            if is_region_valid(df_all, fsz, region_name, REQUIRED_METRICS):
                valid_regions.append(region_name)
            else:
                print(f'[skip] {fsz} / {region_name} skipped (invalid metrics)')
        if not valid_regions:
            continue
        regs = collect_regions_from_bigcsv(df_all, fsz, metric)
        regs = [r for r in regs if r['region'] in valid_regions]
        if not regs:
            continue
        regions_all[fsz] = regs
        values_all.extend([r['value'] for r in regs if np.isfinite(r['value'])])
    if not regions_all or not values_all:
        raise ValueError('No finite values for this metric in the evaluated FOVs')
    vmin = float(np.nanmin(values_all)) if vmin is None else vmin
    vmax = float(np.nanmax(values_all)) if vmax is None else vmax
    norm = Normalize(vmin=vmin, vmax=vmax)
    n = len(fov_sizes)
    nrows, ncols = (1, 1) if n == 1 else (1, n)
    if rotate_cw90:
        img_height, img_width = (he_np.shape[1], he_np.shape[0])
    else:
        img_height, img_width = he_np.shape[:2]
    img_aspect = img_width / img_height
    if nrows == 1 and ncols > 1:
        available_width_mm = FIG_WIDTH_MM - LEFT_MARGIN_MM - RIGHT_MARGIN_MM - CBAR_WIDTH_MM_VERTICAL - CBAR_GAP_MM
    else:
        available_width_mm = FIG_WIDTH_MM - LEFT_MARGIN_MM - RIGHT_MARGIN_MM
    panel_width_mm = available_width_mm / ncols
    panel_width_in = panel_width_mm / MM_PER_INCH
    panel_height_in = panel_width_in / img_aspect
    panel_height_mm = panel_height_in * MM_PER_INCH
    BOTTOM_SPACE_MM = 8.0
    if nrows == 1 and ncols > 1:
        FIG_HEIGHT_MM = panel_height_mm + BOTTOM_SPACE_MM + 5
    else:
        FIG_HEIGHT_MM = panel_height_mm + BOTTOM_SPACE_MM + CBAR_HEIGHT_MM_HORIZONTAL + CBAR_GAP_MM + 3
    FIG_HEIGHT_IN = FIG_HEIGHT_MM / MM_PER_INCH
    fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))
    left_margin_frac = LEFT_MARGIN_MM / FIG_WIDTH_MM
    panel_width_frac = panel_width_mm / FIG_WIDTH_MM
    panel_height_frac = panel_height_mm / FIG_HEIGHT_MM
    axes = []
    for i in range(ncols):
        ax_left = left_margin_frac + i * panel_width_frac
        ax_bottom = BOTTOM_SPACE_MM / FIG_HEIGHT_MM
        ax = fig.add_axes([ax_left, ax_bottom, panel_width_frac, panel_height_frac])
        axes.append(ax)
    if n == 1:
        axes = axes
    for ax, fsz in zip(axes, fov_sizes):
        if fsz not in regions_all:
            ax.set_axis_off()
            ax.set_title(f'FOV size: {fsz} pixels\nNot evaluable', fontsize=TITLE_FONTSIZE)
            continue
        draw_fov_panel(ax, he_np, regions_all[fsz], metric, fsz, highlight_index, box_draw_size, text_pref, norm, cmap, rotate_cw90=rotate_cw90, title_fontsize=TITLE_FONTSIZE, text_fontsize=TEXT_FONTSIZE)
    if nrows == 1 and ncols > 1:
        cbar_left = (LEFT_MARGIN_MM + available_width_mm + CBAR_GAP_MM) / FIG_WIDTH_MM
        cbar_width = CBAR_WIDTH_MM_VERTICAL / FIG_WIDTH_MM
        cbar_bottom = BOTTOM_SPACE_MM / FIG_HEIGHT_MM
        cbar_height = panel_height_frac
        cax = fig.add_axes([cbar_left, cbar_bottom, cbar_width, cbar_height])
        orientation = 'vertical'
    else:
        cbar_left = left_margin_frac
        cbar_width = panel_width_frac
        cbar_height = CBAR_HEIGHT_MM_HORIZONTAL / FIG_HEIGHT_MM
        cbar_bottom = BOTTOM_SPACE_MM / FIG_HEIGHT_MM - CBAR_GAP_MM / FIG_HEIGHT_MM - cbar_height
        cax = fig.add_axes([cbar_left, cbar_bottom, cbar_width, cbar_height])
        orientation = 'horizontal'
    extend = 'max' if higher_is_better else 'min'
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation=orientation, extend=extend)
    ticks = np.linspace(vmin, vmax, 6)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{t:.3f}' for t in ticks])
    if orientation == 'horizontal':
        cbar.ax.tick_params(labelsize=CBAR_TICK_FONTSIZE, rotation=90)
        for label in cbar.ax.get_xticklabels():
            label.set_ha('center')
            label.set_va('top')
    else:
        cbar.ax.tick_params(labelsize=CBAR_TICK_FONTSIZE, rotation=0)
        for label in cbar.ax.get_yticklabels():
            label.set_ha('left')
            label.set_va('center')
    if fig_name:
        plt.savefig(fig_name, format='pdf', dpi=300)
    plt.show()
    print(f'\nFigure saved: {fig_name}')
    print(f'   PDF width: {FIG_WIDTH_MM} mm (fixed)')
    print(f'   Image fill: {available_width_mm / FIG_WIDTH_MM * 100:.1f}% of width')
    print(f"   Colorbar size: {(CBAR_WIDTH_MM_VERTICAL if orientation == 'vertical' else CBAR_HEIGHT_MM_HORIZONTAL)} mm")
    print(f'   Image aspect ratio: {img_aspect:.3f}')
    print(f'   After scaling to {FINAL_WIDTH_MM} mm:')
    print(f'   - Title font: ~{TARGET_TITLE_PT} pt')
    print(f'   - Colorbar ticks: ~{TARGET_CBAR_TICK_PT} pt')

