#!/usr/bin/env python3
"""Quantify IF signal from multiplex OME-TIFF images.

This script implements the two IF quantifications used in the manuscript:
1. pSTAT3 binning in the pSTAT3 panel.
2. ECM compartment assignment of CD14+OPN+ events in the fibrin/collagen panel.

Expected inputs:
- One figure-parameter CSV for the pSTAT3 panel.
- One figure-parameter CSV for the fibrin/collagen panel.
- A directory containing the corresponding OME-TIFF files.

Main outputs:
- if_signal_quantification_summary.csv
- pstat3_seed_calls.csv
- ecm_seed_calls.csv
- qc/*.mask_qc.png

Method summary:
- DAPI is used to define nucleus-anchored events.
- CD14 and OPN are assigned in perinuclear neighborhoods.
- pSTAT3 is assigned in the nuclear neighborhood.
- Fibrin and CollagenI are assigned as local ECM-region signals.
- Optional panel-specific coarse-mask thresholds can be supplied to suppress background artifact.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import tifffile as tf
from scipy.ndimage import (
    binary_closing,
    binary_dilation,
    binary_fill_holes,
    binary_opening,
    distance_transform_edt,
    find_objects,
    gaussian_filter,
    label,
    maximum_filter,
    uniform_filter,
    zoom,
)


CHANNEL_INDEX = {
    'DAPI': 0,
    'Alexa488': 1,
    'Alexa594': 2,
    'Alexa647': 3,
    'CFP': 4,
}


@dataclass
class Seed:
    seed_id: int
    y: int
    x: int
    panel_id: str
    cd14: bool = False
    opn: bool = False
    pstat3: bool = False
    fibrin_local: bool = False
    collagen_local: bool = False


@dataclass
class PanelData:
    panel_id: str
    ome_file: str
    vmaps: dict[str, np.ndarray]
    masks: dict[str, np.ndarray]
    tissue_mask: np.ndarray
    dapi_binary: np.ndarray
    dapi_nucleus_mask: np.ndarray
    seeds: list[Seed]


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description='Quantify IF signal from OME-TIFFs using figure-parameter CSVs.')
    p.add_argument('--pstat3-csv', type=Path, required=True, help='Figure-parameter CSV for the pSTAT3 panel.')
    p.add_argument('--ecm-csv', type=Path, required=True, help='Figure-parameter CSV for the fibrin/collagen panel.')
    p.add_argument('--ome-search-root', type=Path, default=root / 'final_images' / 'raw_ome_data', help='Root directory searched recursively for matching OME-TIFFs.')
    p.add_argument('--out-dir', type=Path, default=root / 'if_signal_quantification_output', help='Output directory for quantification tables and QC images.')
    p.add_argument('--pstat3-mask-threshold', type=float, default=None, help='Optional coarse-mask threshold for the pSTAT3 panel.')
    p.add_argument('--ecm-mask-threshold', type=float, default=None, help='Optional coarse-mask threshold for the ECM panel.')
    p.add_argument('--seed-min-distance', type=int, default=3)
    p.add_argument('--nuclear-radius-px', type=int, default=6)
    p.add_argument('--perinuclear-radius-px', type=int, default=16)
    p.add_argument('--ecm-radius-px', type=int, default=22)
    return p.parse_args()


def canonical_marker_name(name: str) -> str:
    n = name.strip().lower().replace('_', '').replace('-', '').replace(' ', '')
    aliases = {
        'dapi': 'DAPI',
        'cd14': 'CD14',
        'opn': 'OPN',
        'spp1': 'OPN',
        'pstat3': 'pSTAT3',
        'stat3': 'pSTAT3',
        'fibrin': 'Fibrin',
        'collageni': 'CollagenI',
        'collagen1': 'CollagenI',
        'col1': 'CollagenI',
    }
    return aliases.get(n, name.strip())


def percentile_norm(arr: np.ndarray, lo: float = 1.0, hi: float = 99.7) -> np.ndarray:
    a = arr.astype(np.float32)
    p_lo = float(np.percentile(a, lo))
    p_hi = float(np.percentile(a, hi))
    if p_hi <= p_lo:
        p_hi = p_lo + 1.0
    out = (a - p_lo) / (p_hi - p_lo)
    return np.clip(out, 0.0, 1.0)


def transform(raw_arr: np.ndarray, mn: float, mx: float, gamma: float, gain: float) -> np.ndarray:
    mn = max(0.0, float(mn))
    mx = max(mn + 1.0, float(mx))
    gamma = max(0.3, min(3.0, float(gamma)))
    gain = max(0.0, min(3.0, float(gain)))
    v = (raw_arr.astype(np.float32) - mn) / (mx - mn)
    v = np.clip(v, 0.0, 1.0)
    v = np.power(v, gamma) * gain
    return np.clip(v, 0.0, 1.0)


def fill_small_holes(mask: np.ndarray, max_area: int) -> np.ndarray:
    filled = binary_fill_holes(mask)
    holes = filled & (~mask)
    lbl, n = label(holes)
    if n == 0:
        return mask
    areas = np.bincount(lbl.ravel())
    keep = np.zeros(n + 1, dtype=bool)
    for cid in range(1, n + 1):
        if int(areas[cid]) <= int(max_area):
            keep[cid] = True
    return mask | keep[lbl]


def resolve_ome_path(search_root: Path, ome_file: str) -> Path:
    matches = sorted(search_root.rglob(ome_file))
    if not matches:
        raise FileNotFoundError(f'Could not resolve OME-TIFF for {ome_file} under {search_root}')
    if len(matches) > 1:
        joined = '\n'.join(str(m) for m in matches)
        raise RuntimeError(f'Ambiguous OME-TIFF match for {ome_file} under {search_root}. Matches:\n{joined}')
    return matches[0]


def select_series_level(series: tf.TiffPageSeries, preferred_level: int) -> tuple[np.ndarray, int]:
    levels = list(series.levels) if getattr(series, 'levels', None) else []
    if not levels:
        return series.asarray(), 0
    level_index = min(max(0, int(preferred_level)), len(levels) - 1)
    return levels[level_index].asarray(), level_index


def infer_required_markers(panel_id: str, ome_file: str, csv_path: Path) -> set[str]:
    key = f'{panel_id} {ome_file} {csv_path.name}'.lower()
    if 'collagen' in key:
        return {'DAPI', 'CD14', 'OPN', 'Fibrin', 'CollagenI'}
    if 'pstat3' in key:
        return {'DAPI', 'CD14', 'OPN', 'pSTAT3'}
    raise ValueError(f'Could not infer required markers for {csv_path}; expected a pSTAT3 or CollagenI panel.')


def load_params(csv_path: Path) -> tuple[str, str, dict[str, dict]]:
    panel_id = ''
    ome_file = ''
    params: dict[str, dict] = {}
    with csv_path.open('r', encoding='utf-8', newline='') as f:
        for row in csv.DictReader(f):
            if row.get('item_type') != 'marker':
                continue
            panel_id = row['panel_id']
            ome_file = row['ome_file']
            marker = canonical_marker_name(row['marker'])
            params[marker] = {
                'channel': row['channel'],
                'min': float(row['min']),
                'max': float(row['max']),
                'gamma': float(row['gamma']),
                'gain': float(row['intensity']),
            }
    if not params:
        raise ValueError(f'No marker rows found in {csv_path}')
    return panel_id, ome_file, params


def compute_af_max(raw: dict[str, np.ndarray]) -> np.ndarray:
    af_stack = []
    for arr in raw.values():
        x = arr.astype(np.float32)
        p99 = float(np.percentile(x, 99.5))
        if p99 <= 1e-6:
            p99 = 1.0
        af_stack.append(np.clip(x / p99, 0.0, 1.0))
    return np.max(np.stack(af_stack, axis=0), axis=0)


def build_baseline_tissue_mask(raw: dict[str, np.ndarray]) -> np.ndarray:
    dapi = percentile_norm(raw['DAPI'], 1.0, 99.5)
    others = [percentile_norm(arr, 1.0, 99.5) for name, arr in raw.items() if name != 'DAPI']
    nondapi_max = np.max(np.stack(others, axis=0), axis=0)
    dapi_fine = gaussian_filter(dapi, sigma=0.45)
    dapi_bg = gaussian_filter(dapi, sigma=2.2)
    dapi_lc = np.clip(dapi_fine - dapi_bg, 0.0, 1.0)
    nondapi_smooth = uniform_filter(nondapi_max, size=5, mode='nearest')
    nondapi_coarse = uniform_filter(nondapi_max, size=17, mode='nearest')
    nondapi_broad = uniform_filter(nondapi_max, size=33, mode='nearest')
    dapi_lc_smooth = uniform_filter(dapi_lc, size=5, mode='nearest')
    strong = (nondapi_smooth >= 0.09) | (nondapi_coarse >= 0.06) | (dapi_lc >= 0.026)
    support = (nondapi_coarse >= 0.05) | (nondapi_broad >= 0.04) | (dapi_lc_smooth >= 0.013)
    seed = binary_closing(strong, structure=np.ones((3, 3), dtype=bool))
    seed = fill_small_holes(seed, max_area=200)
    grown = binary_dilation(seed, structure=np.ones((5, 5), dtype=bool))
    candidate = grown & support
    candidate = binary_opening(candidate, structure=np.ones((2, 2), dtype=bool))
    candidate = binary_closing(candidate, structure=np.ones((3, 3), dtype=bool))
    candidate = fill_small_holes(candidate, max_area=500)
    lbl, n = label(candidate)
    if n == 0:
        return candidate
    areas = np.bincount(lbl.ravel())
    keep = np.zeros(n + 1, dtype=bool)
    largest = float(np.max(areas[1:])) if n else 0.0
    for cid in range(1, n + 1):
        comp = lbl == cid
        area = int(areas[cid])
        strong_frac = float(strong[comp].mean()) if np.any(comp) else 0.0
        support_frac = float(support[comp].mean()) if np.any(comp) else 0.0
        border_touch = bool(comp[0, :].any() or comp[-1, :].any() or comp[:, 0].any() or comp[:, -1].any())
        if area < 60:
            continue
        if area < largest * 0.03 and strong_frac < 0.10:
            continue
        if border_touch and strong_frac < 0.16 and support_frac < 0.28:
            continue
        if support_frac < 0.18:
            continue
        keep[cid] = True
    return keep[lbl]


def build_special_tissue_mask(ome_path: Path, target_shape: tuple[int, int], threshold: float) -> np.ndarray:
    with tf.TiffFile(ome_path) as tif:
        lvl_arr, lvl_index = select_series_level(tif.series[0], preferred_level=6)
        lvl6 = lvl_arr.astype(np.float32)
    d6 = percentile_norm(lvl6[0], 1.0, 99.5)
    n6 = np.max(np.stack([percentile_norm(lvl6[i], 1.0, 99.5) for i in range(1, lvl6.shape[0])], axis=0), axis=0)
    mix6 = np.maximum(gaussian_filter(n6, 1.2), gaussian_filter(d6, 1.0) * 0.8)
    mask6 = mix6 >= float(threshold)
    mask6 = binary_closing(mask6, structure=np.ones((5, 5), dtype=bool))
    mask6 = binary_opening(mask6, structure=np.ones((3, 3), dtype=bool))
    mask6 = binary_fill_holes(mask6)
    scale = (target_shape[0] / mask6.shape[0], target_shape[1] / mask6.shape[1])
    mask4 = zoom(mask6.astype(np.float32), scale, order=0) > 0.5
    mask4 = binary_closing(mask4, structure=np.ones((3, 3), dtype=bool))
    mask4 = binary_fill_holes(mask4)
    return mask4


def build_dapi_local_contrast(dapi_vmap: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    fine = gaussian_filter(dapi_vmap, sigma=0.45)
    bg = gaussian_filter(dapi_vmap, sigma=2.2)
    local_contrast = np.clip(fine - bg, 0.0, 1.0)
    return fine, local_contrast


def build_nucleus_like_dapi(dapi_vmap: np.ndarray, tissue_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    fine, local_contrast = build_dapi_local_contrast(dapi_vmap)
    tissue_vals = local_contrast[tissue_mask]
    positive = tissue_vals[tissue_vals > 0]
    if positive.size:
        lc_floor = max(0.001, float(np.percentile(positive, 12)))
        seed_floor = max(0.004, float(np.percentile(positive, 48)))
    else:
        lc_floor = 0.001
        seed_floor = 0.004
    binary = (local_contrast >= lc_floor) & tissue_mask & (fine >= 0.02)
    binary = binary_opening(binary, structure=np.ones((2, 2), dtype=bool))
    binary = binary_closing(binary, structure=np.ones((2, 2), dtype=bool))
    binary = binary_fill_holes(binary)
    lbl, n = label(binary)
    keep = np.zeros(n + 1, dtype=bool)
    areas = np.bincount(lbl.ravel())
    for cid in range(1, n + 1):
        if int(areas[cid]) >= 3:
            keep[cid] = True
    nucleus_mask = keep[lbl]
    return binary, nucleus_mask, seed_floor


def detect_dapi_seeds(dapi_vmap: np.ndarray, nucleus_mask: np.ndarray, seed_floor: float, seed_min_distance: int) -> list[tuple[int, int]]:
    fine, local_contrast = build_dapi_local_contrast(dapi_vmap)
    dist = distance_transform_edt(nucleus_mask)
    lbl, n = label(nucleus_mask)
    areas = np.bincount(lbl.ravel())
    boxes = find_objects(lbl)
    small_areas = [int(a) for a in areas[1:] if 20 <= int(a) <= 180]
    typical_area = float(np.median(small_areas)) if small_areas else 55.0
    typical_area = max(30.0, min(120.0, typical_area))
    seeds: list[tuple[int, int]] = []
    min_dist2 = float(seed_min_distance ** 2)
    for cid, box in enumerate(boxes, start=1):
        if box is None:
            continue
        area = int(areas[cid])
        if area < 10:
            continue
        sub_lbl = lbl[box] == cid
        sub_dist = dist[box]
        score = sub_dist * (0.20 + fine[box] + 2.0 * local_contrast[box])
        size = max(3, int(seed_min_distance) * 2 + 1)
        local_max = score == maximum_filter(score, size=size, mode='nearest')
        peak_mask = local_max & sub_lbl & (local_contrast[box] >= seed_floor) & (sub_dist >= 0.8)
        ys, xs = np.where(peak_mask)
        order = np.argsort(score[ys, xs])[::-1] if len(ys) else np.array([], dtype=int)
        comp_seeds: list[tuple[int, int]] = []
        for idx in order:
            y = int(ys[idx])
            x = int(xs[idx])
            if any((y - py) ** 2 + (x - px) ** 2 < min_dist2 for py, px in comp_seeds):
                continue
            comp_seeds.append((y, x))
        desired = max(1, int(round(area / typical_area)))
        if len(comp_seeds) < desired:
            ys2, xs2 = np.where(sub_lbl)
            order2 = np.argsort(score[ys2, xs2])[::-1]
            for idx in order2:
                y = int(ys2[idx])
                x = int(xs2[idx])
                if sub_dist[y, x] < 0.6:
                    continue
                if any((y - py) ** 2 + (x - px) ** 2 < min_dist2 for py, px in comp_seeds):
                    continue
                comp_seeds.append((y, x))
                if len(comp_seeds) >= desired:
                    break
        for y, x in comp_seeds:
            seeds.append((y + box[0].start, x + box[1].start))
    return seeds


def neighborhood_mean(mask: np.ndarray, radius: int) -> np.ndarray:
    size = max(3, int(radius) * 2 + 1)
    return uniform_filter(mask.astype(np.float32), size=size, mode='nearest')


def ring_mean(inner_mean: np.ndarray, outer_mean: np.ndarray, inner_radius: int, outer_radius: int) -> np.ndarray:
    inner_area = float((2 * inner_radius + 1) ** 2)
    outer_area = float((2 * outer_radius + 1) ** 2)
    ring_area = max(1.0, outer_area - inner_area)
    return np.clip((outer_mean * outer_area - inner_mean * inner_area) / ring_area, 0.0, 1.0)


def classify_seeds(panel_id: str, masks: dict[str, np.ndarray], seeds_xy: list[tuple[int, int]], args: argparse.Namespace) -> list[Seed]:
    seeds: list[Seed] = []
    nuc_r = int(args.nuclear_radius_px)
    peri_r = int(args.perinuclear_radius_px)
    ecm_r = int(args.ecm_radius_px)
    cache: dict[str, np.ndarray] = {}
    for name, mask in masks.items():
        if name.startswith('DAPI'):
            continue
        cache[f'{name}_nuc'] = neighborhood_mean(mask, nuc_r)
        cache[f'{name}_peri_outer'] = neighborhood_mean(mask, peri_r)
        cache[f'{name}_ecm'] = neighborhood_mean(mask, ecm_r)
    peri_cd14 = ring_mean(cache['CD14_nuc'], cache['CD14_peri_outer'], nuc_r, peri_r) if 'CD14_peri_outer' in cache else None
    peri_opn = ring_mean(cache['OPN_nuc'], cache['OPN_peri_outer'], nuc_r, peri_r) if 'OPN_peri_outer' in cache else None
    for i, (y, x) in enumerate(seeds_xy, start=1):
        seed = Seed(seed_id=i, y=y, x=x, panel_id=panel_id)
        if peri_cd14 is not None:
            seed.cd14 = float(peri_cd14[y, x]) >= 0.08
        if peri_opn is not None:
            seed.opn = float(peri_opn[y, x]) >= 0.03
        if 'pSTAT3_nuc' in cache:
            seed.pstat3 = float(cache['pSTAT3_nuc'][y, x]) >= 0.45
        if 'Fibrin_ecm' in cache:
            seed.fibrin_local = float(cache['Fibrin_ecm'][y, x]) >= 0.25
        if 'CollagenI_ecm' in cache:
            seed.collagen_local = float(cache['CollagenI_ecm'][y, x]) >= 0.85
        seeds.append(seed)
    return seeds


def build_panel_data(csv_path: Path, ome_search_root: Path, special_mask_threshold: float | None, args: argparse.Namespace) -> PanelData:
    panel_id, ome_file, params = load_params(csv_path)
    required_markers = infer_required_markers(panel_id, ome_file, csv_path)
    ome_path = resolve_ome_path(ome_search_root, ome_file)
    with tf.TiffFile(ome_path) as tif:
        arr, level_index = select_series_level(tif.series[0], preferred_level=4)
    raw: dict[str, np.ndarray] = {}
    vmaps: dict[str, np.ndarray] = {}
    masks: dict[str, np.ndarray] = {}
    unknown_channels = sorted({p['channel'] for p in params.values() if p['channel'] not in CHANNEL_INDEX})
    if unknown_channels:
        raise RuntimeError(f'Unrecognized channel names in {csv_path}: {unknown_channels}')
    for marker, p in params.items():
        channel_name = p['channel']
        raw_marker = arr[CHANNEL_INDEX[channel_name]]
        raw[marker] = raw_marker
        if marker == 'DAPI':
            vmaps[marker] = percentile_norm(raw_marker, 1.0, 99.7)
        else:
            vmaps[marker] = transform(raw_marker, p['min'], p['max'], p['gamma'], p['gain'])
            masks[marker] = vmaps[marker] >= 0.25
    missing_markers = sorted(required_markers - set(params))
    if missing_markers:
        raise RuntimeError(f'Missing required markers in {csv_path}: {missing_markers}')
    if special_mask_threshold is None:
        tissue_mask = build_baseline_tissue_mask(raw)
    else:
        tissue_mask = build_special_tissue_mask(ome_path, raw['DAPI'].shape, special_mask_threshold)
    dapi_binary, dapi_nucleus_mask, seed_floor = build_nucleus_like_dapi(vmaps['DAPI'], tissue_mask)
    seeds_xy = detect_dapi_seeds(vmaps['DAPI'], dapi_nucleus_mask, seed_floor, args.seed_min_distance)
    masks = {k: (m & tissue_mask) for k, m in masks.items()}
    masks['DAPI'] = dapi_nucleus_mask & tissue_mask
    seeds = classify_seeds(panel_id, masks, seeds_xy, args)
    return PanelData(panel_id=panel_id, ome_file=ome_file, vmaps=vmaps, masks=masks, tissue_mask=tissue_mask, dapi_binary=dapi_binary, dapi_nucleus_mask=dapi_nucleus_mask, seeds=seeds)


def save_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text('', encoding='utf-8')
        return
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row.keys():
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def save_mask_qc(panel: PanelData, out_path: Path) -> None:
    fig, ax = plt.subplots(1, 4, figsize=(14, 4), dpi=160)
    ax[0].imshow(panel.vmaps['DAPI'], cmap='gray', vmin=0, vmax=1)
    ax[0].set_title('DAPI')
    ax[1].imshow(panel.dapi_binary, cmap='gray', vmin=0, vmax=1)
    ax[1].set_title('DAPI binary')
    ax[2].imshow(panel.tissue_mask, cmap='gray', vmin=0, vmax=1)
    ax[2].set_title('Tissue mask')
    overlay = np.stack([panel.vmaps['DAPI'] * 0.4] * 3, axis=-1)
    overlay[..., 1] = np.clip(overlay[..., 1] + panel.tissue_mask * 0.6, 0.0, 1.0)
    ax[3].imshow(overlay)
    ax[3].set_title(f'Overlay ({len(panel.seeds)} seeds)')
    for a in ax:
        a.axis('off')
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches='tight')
    plt.close(fig)


def summarize_pstat3_bins(panel: PanelData) -> list[dict]:
    seeds = panel.seeds
    pstat3_pos = [s for s in seeds if s.pstat3]
    cd14_opn = [s for s in seeds if s.cd14 and s.opn]
    cd14_not_opn = [s for s in seeds if s.cd14 and not s.opn]
    cd14_neg = [s for s in seeds if not s.cd14]
    rows = [
        {'analysis': 'pstat3_within_group', 'group': 'CD14+OPN+', 'n_group': len(cd14_opn), 'n_pstat3_pos': sum(s.pstat3 for s in cd14_opn), 'proportion': sum(s.pstat3 for s in cd14_opn) / max(1, len(cd14_opn))},
        {'analysis': 'pstat3_within_group', 'group': 'CD14+OPN-', 'n_group': len(cd14_not_opn), 'n_pstat3_pos': sum(s.pstat3 for s in cd14_not_opn), 'proportion': sum(s.pstat3 for s in cd14_not_opn) / max(1, len(cd14_not_opn))},
        {'analysis': 'pstat3_within_group', 'group': 'CD14-', 'n_group': len(cd14_neg), 'n_pstat3_pos': sum(s.pstat3 for s in cd14_neg), 'proportion': sum(s.pstat3 for s in cd14_neg) / max(1, len(cd14_neg))},
    ]
    composition = {
        'CD14+OPN+': sum(1 for s in pstat3_pos if s.cd14 and s.opn),
        'CD14+OPN-': sum(1 for s in pstat3_pos if s.cd14 and not s.opn),
        'CD14-OPN+': sum(1 for s in pstat3_pos if (not s.cd14) and s.opn),
        'CD14-OPN-': sum(1 for s in pstat3_pos if (not s.cd14) and (not s.opn)),
    }
    for group, n_group in composition.items():
        rows.append({'analysis': 'pstat3_composition', 'group': group, 'n_group': n_group, 'n_pstat3_pos': len(pstat3_pos), 'proportion': n_group / max(1, len(pstat3_pos))})
    rows.append({'analysis': 'pstat3_total', 'group': 'all_cells', 'n_group': len(seeds), 'n_pstat3_pos': len(pstat3_pos), 'proportion': len(pstat3_pos) / max(1, len(seeds))})
    return rows


def summarize_opnmac_ecm(panel: PanelData) -> list[dict]:
    opnmac = [s for s in panel.seeds if s.cd14 and s.opn]
    fibrin_only = sum(1 for s in opnmac if s.fibrin_local and not s.collagen_local)
    collagen_only = sum(1 for s in opnmac if s.collagen_local and not s.fibrin_local)
    both = sum(1 for s in opnmac if s.fibrin_local and s.collagen_local)
    neither = sum(1 for s in opnmac if not s.fibrin_local and not s.collagen_local)
    total = max(1, len(opnmac))
    return [
        {'analysis': 'opnmac_ecm_compartment', 'group': 'Fibrin only', 'n_group': fibrin_only, 'n_opnmac': len(opnmac), 'proportion': fibrin_only / total},
        {'analysis': 'opnmac_ecm_compartment', 'group': 'Collagen only', 'n_group': collagen_only, 'n_opnmac': len(opnmac), 'proportion': collagen_only / total},
        {'analysis': 'opnmac_ecm_compartment', 'group': 'Both', 'n_group': both, 'n_opnmac': len(opnmac), 'proportion': both / total},
        {'analysis': 'opnmac_ecm_compartment', 'group': 'Neither', 'n_group': neither, 'n_opnmac': len(opnmac), 'proportion': neither / total},
    ]


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    qc_dir = args.out_dir / 'qc'
    qc_dir.mkdir(exist_ok=True)

    pstat3_panel = build_panel_data(args.pstat3_csv, args.ome_search_root, args.pstat3_mask_threshold, args)
    ecm_panel = build_panel_data(args.ecm_csv, args.ome_search_root, args.ecm_mask_threshold, args)

    save_mask_qc(pstat3_panel, qc_dir / f'{pstat3_panel.panel_id}.mask_qc.png')
    save_mask_qc(ecm_panel, qc_dir / f'{ecm_panel.panel_id}.mask_qc.png')

    pstat3_seed_rows = []
    for s in pstat3_panel.seeds:
        pstat3_seed_rows.append({'panel_id': s.panel_id, 'seed_id': s.seed_id, 'x': s.x, 'y': s.y, 'cd14': int(s.cd14), 'opn': int(s.opn), 'pstat3': int(s.pstat3)})
    ecm_seed_rows = []
    for s in ecm_panel.seeds:
        ecm_seed_rows.append({'panel_id': s.panel_id, 'seed_id': s.seed_id, 'x': s.x, 'y': s.y, 'cd14': int(s.cd14), 'opn': int(s.opn), 'fibrin_local': int(s.fibrin_local), 'collagen_local': int(s.collagen_local)})

    save_csv(args.out_dir / 'pstat3_seed_calls.csv', pstat3_seed_rows)
    save_csv(args.out_dir / 'ecm_seed_calls.csv', ecm_seed_rows)
    save_csv(args.out_dir / 'if_signal_quantification_summary.csv', summarize_pstat3_bins(pstat3_panel) + summarize_opnmac_ecm(ecm_panel))

    (args.out_dir / 'README.txt').write_text(
        '\n'.join([
            'Quantification of IF signal',
            '',
            'Inputs:',
            f'- pSTAT3 panel CSV: {args.pstat3_csv}',
            f'- ECM panel CSV: {args.ecm_csv}',
            f'- OME search root: {args.ome_search_root}',
            '',
            'Outputs:',
            '- if_signal_quantification_summary.csv',
            '- pstat3_seed_calls.csv',
            '- ecm_seed_calls.csv',
            '- qc/*.mask_qc.png',
            '',
            'Notes:',
            '- DAPI is used to define nucleus-anchored events.',
            '- CD14 and OPN are assigned in perinuclear neighborhoods.',
            '- pSTAT3 is assigned in the nuclear neighborhood.',
            '- Fibrin and CollagenI are assigned as local ECM-region signals.',
            '- Optional coarse-mask thresholds can be provided per panel to suppress background artifact.',
        ]) + '\n', encoding='utf-8'
    )
    print(args.out_dir)


if __name__ == '__main__':
    main()
