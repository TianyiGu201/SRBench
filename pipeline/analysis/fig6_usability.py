"""Normalize usability observations and plot method comparisons."""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from typing import Optional


BIN_METRICS = {'Installation', 'Input', 'Purpose/Usage', 'Output', 'Peer review', 'Scale'}
TRI_METRICS = {'Analysis', 'Reliability', 'Parameters'}
CONT_METRICS = {'Activity', 'Response', 'PArobustness', 'SArobustness'}

def _read_usability_table(path: str) -> pd.DataFrame:
    """
    Robust reader for:
      - comma/tab separated CSV
      - whitespace-separated table (ragged)
    Assumes first column is method/algorithm name.
    Keeps missing as NaN.
    """
    try:
        df = pd.read_csv(path, index_col=0)
        if df.shape[1] >= 3:
            return df
    except Exception:
        pass
    df = pd.read_csv(path, sep=None, engine='python', index_col=0)
    return df

def _score_response_closeness(x: pd.Series, target: float=30.0) -> pd.Series:
    """
    Response: closer to 30 is better -> score in [0,1].
    score = 1 - |x-target| / max(|x-target|)
    If all finite values equal to target -> score=1 for finite.
    """
    x = pd.to_numeric(x, errors='coerce')
    dev = (x - target).abs()
    if not np.isfinite(dev.values).any():
        return pd.Series(np.nan, index=x.index)
    max_dev = np.nanmax(dev.values)
    if np.isclose(max_dev, 0):
        return pd.Series(np.where(np.isfinite(x.values), 1.0, np.nan), index=x.index)
    s = 1.0 - dev / max_dev
    return s.clip(0, 1)

def _minmax01(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors='coerce')
    if not np.isfinite(x.values).any():
        return pd.Series(np.nan, index=x.index)
    vmin, vmax = (np.nanmin(x.values), np.nanmax(x.values))
    if np.isclose(vmin, vmax):
        return pd.Series(np.nan, index=x.index)
    return (x - vmin) / (vmax - vmin)

def compute_normalized_scores(df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a 0..1 normalized score table (same shape), with rules:
      - Binary metrics: 0->0, 1->1 (NaN preserved)
      - Ternary metrics: 0->0, 0.5->0.5, 1->1 (NaN preserved)
      - Continuous:
          * Activity, PArobustness, SArobustness: min-max across methods
          * Response: closeness-to-30 score in [0,1] (no min-max after)
    """
    out = pd.DataFrame(index=df.index)
    for col in df.columns:
        if col in BIN_METRICS:
            x = pd.to_numeric(df[col], errors='coerce')
            out[col] = x.where(x.isna(), x.clip(0, 1))
        elif col in TRI_METRICS:
            x = pd.to_numeric(df[col], errors='coerce')
            out[col] = x.where(x.isna(), x.clip(0, 1))
        elif col == 'Response':
            out[col] = _score_response_closeness(df[col], target=30.0)
        elif col in {'Activity', 'PArobustness', 'SArobustness'}:
            out[col] = _minmax01(df[col])
        else:
            out[col] = _minmax01(df[col])
    return out

def _rgba(hex_color: str, alpha: float=1.0):
    r, g, b, _ = mcolors.to_rgba(hex_color)
    return (r, g, b, alpha)

def _blend(c1, c2, t: float):
    return tuple((1 - t) * np.array(c1) + t * np.array(c2))

def _make_category_levels(base_hex: str):
    """
    Produce levels using one hue:
      - lvl0: very light (for 0)
      - lvl05: mid (for 0.5)
      - base: base hue
      - lvl1: darker (for 1)
    """
    base = _rgba(base_hex, 1.0)
    white = (1, 1, 1, 1.0)
    dark = (base[0] * 0.75, base[1] * 0.75, base[2] * 0.75, 1.0)
    lvl0 = _blend(white, base, 0.08)
    lvl05 = _blend(white, base, 0.35)
    lvl1 = _blend(base, dark, 0.45)
    return (lvl0, lvl05, base, lvl1)

def _make_soft_cmap(base_hex: str, name: str):
    lvl0, lvl05, base, lvl1 = _make_category_levels(base_hex)
    return mcolors.LinearSegmentedColormap.from_list(name, [lvl0, lvl05, base, lvl1], N=256)

def _value_to_color(val: float, base_hex: str, mode: str):
    """
    mode:
      - 'bin' : {0,1}
      - 'tri' : {0,0.5,1}
      - 'cont': [0,1] continuous
      - 'rank5': 5-level discrete in [0,1]
    """
    lvl0, lvl05, base, lvl1 = _make_category_levels(base_hex)
    if not np.isfinite(val):
        return '#E6E6E6'
    if mode == 'bin':
        return lvl1 if val >= 0.5 else lvl0
    if mode == 'tri':
        if val < 0.25:
            return lvl0
        elif val < 0.75:
            return lvl05
        else:
            return lvl1
    if mode == 'rank5':
        idx = int(np.clip(np.floor(val * 5), 0, 4))
        cmap = _make_soft_cmap(base_hex, 'tmp_rank5')
        pts = np.linspace(0.12, 0.95, 5)
        return cmap(pts[idx])
    cmap = _make_soft_cmap(base_hex, 'tmp')
    return cmap(mcolors.Normalize(0, 1)(val))

def plot_usability_heatmap(df_raw: pd.DataFrame, gap: float=0.18, figsize=(10.5, 4.2), margins=(0.1, 0.78, 0.8, 0.22), save_pdf: Optional[str]=None, show: bool=True):
    """
    UPDATED:
      - overall calculation remains unchanged (uses normalized score table)
      - GitHub response:
          * Activity cells colored using RAW Activity values (already log10(...) results)
            with a "compressed" mapping to keep differences subtle but visible.
          * Response cells colored by distance-to-30 in RAW space (Response = 30 - median response time),
            so closer to 30 is darker (better), far is lighter; supports negative values naturally.
      - Activity and Response each get their own colorbar, using the original GitHub response cmap hue.
    """
    categories = [('Tutorial', ['Installation', 'Analysis', 'Reliability']), ('Function documentation', ['Input', 'Purpose/Usage', 'Parameters', 'Output']), ('Peer review', ['Peer review']), ('GitHub response', ['Activity', 'Response']), ('Downstream', ['Scale']), ('Stability', ['PArobustness', 'SArobustness'])]
    cat_color = {'Tutorial': '#4C78A8', 'Function documentation': '#F2A541', 'Peer review': '#7F7F7F', 'GitHub response': '#59A14F', 'Downstream': '#E15759', 'Stability': '#76B7B2', 'Usability_overall': '#B279A2'}
    df = df_raw.copy()
    df = df.apply(pd.to_numeric, errors='coerce')
    ordered_cols = []
    col_to_cat = {}
    cat_spans = []
    x = 0
    for cat, cols in categories:
        cols_in = [c for c in cols if c in df.columns]
        if not cols_in:
            continue
        start = x
        for c in cols_in:
            ordered_cols.append(c)
            col_to_cat[c] = cat
            x += 1
        cat_spans.append((cat, start, x))
    df = df[ordered_cols]
    score = compute_normalized_scores(df)
    overall_raw = np.nanmean(score.values.astype(float), axis=1)
    n = len(score.index)
    if n <= 1:
        overall_rank01 = pd.Series(np.ones(n), index=score.index, dtype=float)
    else:
        r = pd.Series(overall_raw, index=score.index).rank(method='dense', ascending=True)
        rmax = r.max()
        if rmax <= 1:
            overall_rank01 = pd.Series(np.ones(n), index=score.index, dtype=float)
        else:
            overall_rank01 = (r - 1) / (rmax - 1)
    score['Usability_overall'] = overall_rank01
    ordered_cols_plus = ordered_cols + ['Usability_overall']
    col_to_cat['Usability_overall'] = 'Usability_overall'
    cat_spans.append(('Usability_overall', len(ordered_cols), len(ordered_cols) + 1))
    col_mode = {}
    for c in ordered_cols_plus:
        if c in BIN_METRICS:
            col_mode[c] = 'bin'
        elif c in TRI_METRICS:
            col_mode[c] = 'tri'
        elif c in CONT_METRICS:
            col_mode[c] = 'cont'
        elif c == 'Usability_overall':
            col_mode[c] = 'rank5'
        else:
            col_mode[c] = 'cont'
    GH_BASE_HEX = cat_color['GitHub response']
    gh_cmap = _make_soft_cmap(GH_BASE_HEX, 'cmap_github_response')

    def _activity_norm_params(xvals: np.ndarray):
        xvals = xvals[np.isfinite(xvals)]
        if xvals.size == 0:
            return (np.nan, np.nan, 1.0)
        vmin = float(np.nanmin(xvals))
        vmax = float(np.nanmax(xvals))
        if np.isclose(vmin, vmax):
            return (vmin, vmax, 1.0)
        span = vmax - vmin
        pad = max(0.08, span * 1.8)
        vmin_p = vmin - pad
        vmax_p = vmax + pad
        gamma = 1.6
        return (vmin_p, vmax_p, gamma)

    def _response_dist_params(xvals: np.ndarray, target: float=30.0):
        xvals = xvals[np.isfinite(xvals)]
        if xvals.size == 0:
            return (np.nan, np.nan, target)
        d = np.abs(xvals - target)
        dmax = float(np.nanmax(d))
        if np.isclose(dmax, 0):
            dmax = 1.0
        dmax = dmax * 1.15
        return (0.0, dmax, target)
    activity_raw = pd.to_numeric(df['Activity'], errors='coerce').to_numpy() if 'Activity' in df.columns else np.array([])
    response_raw = pd.to_numeric(df['Response'], errors='coerce').to_numpy() if 'Response' in df.columns else np.array([])
    act_vmin, act_vmax, act_gamma = _activity_norm_params(activity_raw)
    dist_vmin, dist_vmax, resp_target = _response_dist_params(response_raw, target=30.0)
    act_norm = None
    if np.isfinite(act_vmin) and np.isfinite(act_vmax) and (not np.isclose(act_vmin, act_vmax)):
        act_norm = mcolors.PowerNorm(gamma=act_gamma, vmin=act_vmin, vmax=act_vmax)
    dist_norm = None
    if np.isfinite(dist_vmax) and dist_vmax > 0:
        dist_norm = plt.Normalize(dist_vmin, dist_vmax)
    gh_cmap_response = gh_cmap.reversed()

    def _activity_raw_to_color(v):
        if not np.isfinite(v):
            return '#E6E6E6'
        if act_norm is None:
            return gh_cmap(0.55)
        return gh_cmap(act_norm(v))

    def _response_raw_to_color(v):
        if not np.isfinite(v):
            return '#E6E6E6'
        if dist_norm is None:
            return gh_cmap_response(0.55)
        d = abs(v - resp_target)
        return gh_cmap_response(dist_norm(d))
    fig = plt.figure(figsize=figsize)
    left, right, top, bottom = margins
    ax = fig.add_axes([left, bottom, right - left, top - bottom])
    n_rows = score.shape[0]
    n_cols = len(ordered_cols_plus)
    box_size = 1 - gap
    for cat, x0, x1 in cat_spans:
        ax.add_patch(patches.Rectangle((x0, 0), x1 - x0, n_rows, facecolor=_rgba(cat_color[cat], 0.1), edgecolor='none', zorder=0))
    for i, method in enumerate(score.index):
        for j, metric in enumerate(ordered_cols_plus):
            cat = col_to_cat[metric]
            if metric == 'Activity' and 'Activity' in df.columns:
                val_raw = float(df.loc[method, 'Activity'])
                face = _activity_raw_to_color(val_raw)
            elif metric == 'Response' and 'Response' in df.columns:
                val_raw = float(df.loc[method, 'Response'])
                face = _response_raw_to_color(val_raw)
            else:
                val = score.loc[method, metric]
                face = _value_to_color(val, cat_color[cat], col_mode[metric])
            ax.add_patch(patches.Rectangle((j + gap / 2, i + gap / 2), box_size, box_size, facecolor=face, edgecolor='none', zorder=2))
            if metric == 'Usability_overall':
                raw_v = overall_raw.loc[method] if isinstance(overall_raw, pd.Series) else overall_raw[i]
                if np.isfinite(raw_v):
                    ax.text(j + 0.5, i + 0.5, f'{raw_v:.2f}', ha='center', va='center', fontsize=5, color='black', zorder=3)
    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.tick_params(left=False, bottom=False)
    ax.set_xticks(np.arange(n_cols) + 0.5)
    display_xticks = [*ordered_cols, 'Usability overall']
    ax.set_xticklabels(display_xticks, rotation=90, ha='center', fontsize=6)
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels(score.index, fontsize=6, va='center')
    for spine in ax.spines.values():
        spine.set_visible(False)
    title_y = -0.55
    for cat, x0, x1 in cat_spans:
        mid = (x0 + x1) / 2
        ax.text(mid, title_y, cat if cat != 'Usability_overall' else 'Usability overall', ha='center', va='bottom', fontsize=6, fontweight='bold', color='black', clip_on=False)
    ax_pos = ax.get_position()
    cbar_h = ax_pos.height
    cbar_y = ax_pos.y0
    cbar_w = 0.012
    cbar_gap = 0.01
    cbar_x0 = ax_pos.x1 + 0.015
    appeared = [s[0] for s in cat_spans]
    seen = set()
    appeared = [c for c in appeared if not (c in seen or seen.add(c))]
    cat_to_cols = {cat: [] for cat in appeared}
    for col in ordered_cols_plus:
        cat_to_cols[col_to_cat[col]].append(col)

    def _is_discrete_col(col: str) -> bool:
        return col in BIN_METRICS or col in TRI_METRICS
    cat_kind = {}
    for cat in appeared:
        if cat == 'Usability_overall':
            cat_kind[cat] = 'rank5'
            continue
        if cat == 'GitHub response':
            cat_kind[cat] = 'skip'
            continue
        cols_here = cat_to_cols.get(cat, [])
        if len(cols_here) == 0:
            cat_kind[cat] = 'cont'
        elif all((_is_discrete_col(c) for c in cols_here)):
            any_tri = any((c in TRI_METRICS for c in cols_here))
            cat_kind[cat] = 'tri' if any_tri else 'bin'
        else:
            cat_kind[cat] = 'cont'
    cbar_items = []
    for cat in appeared:
        if cat == 'GitHub response':
            if 'Activity' in df.columns:
                cbar_items.append(('gh_activity_raw', 'Activity'))
            if 'Response' in df.columns:
                cbar_items.append(('gh_response_dist', 'Response'))
            continue
        cbar_items.append(('category', cat))
    for k, (kind, name) in enumerate(cbar_items):
        cax = fig.add_axes([cbar_x0 + k * (cbar_w + cbar_gap), cbar_y, cbar_w, cbar_h])
        if kind == 'gh_activity_raw':
            if act_norm is None:
                norm = plt.Normalize(0, 1)
                sm = plt.cm.ScalarMappable(cmap=gh_cmap, norm=norm)
                sm.set_array([])
                cb = fig.colorbar(sm, cax=cax)
                cb.set_ticks([])
            else:
                sm = plt.cm.ScalarMappable(cmap=gh_cmap, norm=act_norm)
                sm.set_array([])
                cb = fig.colorbar(sm, cax=cax)
                finite = activity_raw[np.isfinite(activity_raw)]
                if finite.size > 0:
                    v0 = float(np.nanmin(finite))
                    v1 = float(np.nanmedian(finite))
                    v2 = float(np.nanmax(finite))
                    cb.set_ticks([v0, v1, v2])
                    cb.set_ticklabels([f'{v0:.2f}', f'{v1:.2f}', f'{v2:.2f}'])
                else:
                    cb.set_ticks([])
            cb.set_label('Activity (raw)', rotation=270, labelpad=8, fontsize=6, color='black')
            cb.ax.tick_params(labelsize=5, length=2)
            continue
        if kind == 'gh_response_dist':
            if dist_norm is None:
                norm = plt.Normalize(0, 1)
                sm = plt.cm.ScalarMappable(cmap=gh_cmap_response, norm=norm)
                sm.set_array([])
                cb = fig.colorbar(sm, cax=cax)
                cb.set_ticks([])
            else:
                sm = plt.cm.ScalarMappable(cmap=gh_cmap_response, norm=dist_norm)
                sm.set_array([])
                cb = fig.colorbar(sm, cax=cax)
                d0 = 0.0
                d1 = dist_vmax * 0.5
                d2 = dist_vmax
                cb.set_ticks([d0, d1, d2])
                cb.set_ticklabels([f'{d0:.0f}', f'{d1:.0f}', f'{d2:.0f}'])
            cb.set_label('|Response-30| (days)', rotation=270, labelpad=8, fontsize=6, color='black')
            cb.ax.tick_params(labelsize=5, length=2)
            continue
        cat = name
        base_hex = cat_color[cat]
        lvl0, lvl05, base, lvl1 = _make_category_levels(base_hex)
        if cat_kind[cat] == 'bin':
            cmap = mcolors.ListedColormap([lvl0, lvl1])
            bounds = [0.0, 0.5, 1.0]
            norm = mcolors.BoundaryNorm(bounds, cmap.N)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cb = fig.colorbar(sm, cax=cax, boundaries=bounds, ticks=[0.25, 0.75])
            cb.set_ticklabels(['0', '1'])
        elif cat_kind[cat] == 'tri':
            cmap = mcolors.ListedColormap([lvl0, lvl05, lvl1])
            bounds = [0.0, 0.25, 0.75, 1.0]
            norm = mcolors.BoundaryNorm(bounds, cmap.N)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cb = fig.colorbar(sm, cax=cax, boundaries=bounds, ticks=[0.125, 0.5, 0.875])
            cb.set_ticklabels(['0', '0.5', '1'])
        elif cat_kind[cat] == 'rank5':
            cmap_cont = _make_soft_cmap(base_hex, f'cmap_{cat}_cont')
            pts = np.linspace(0.12, 0.95, 5)
            colors5 = [cmap_cont(p) for p in pts]
            cmap = mcolors.ListedColormap(colors5)
            bounds = np.linspace(0, 1, 6)
            norm = mcolors.BoundaryNorm(bounds, cmap.N)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            centers = (bounds[:-1] + bounds[1:]) / 2
            cb = fig.colorbar(sm, cax=cax, boundaries=bounds, ticks=centers)
            cb.set_ticklabels([str(i) for i in range(1, 6)])
        else:
            cmap = _make_soft_cmap(base_hex, f'cmap_{cat}')
            norm = plt.Normalize(0, 1)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cb = fig.colorbar(sm, cax=cax)
            cb.set_ticks([0, 0.5, 1.0])
            cb.set_ticklabels(['0', '0.5', '1'])
        cb.set_label(cat if cat != 'Usability_overall' else 'Usability overall', rotation=270, labelpad=8, fontsize=6, color='black')
        cb.ax.tick_params(labelsize=5, length=2)
    legend_handles = [Patch(facecolor=_rgba(cat_color[c], 0.22), edgecolor='black', linewidth=0.6, label=c if c != 'Usability_overall' else 'Usability overall') for c in appeared]
    legend_x = cbar_x0 + len(cbar_items) * (cbar_w + cbar_gap) + 0.015
    fig.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(legend_x, ax_pos.y1), frameon=False, fontsize=6, title='Category', title_fontsize=6)
    if save_pdf:
        fig.savefig(save_pdf, bbox_inches='tight', pad_inches=0.01)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return (fig, ax)

