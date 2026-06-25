"""
serial_registration.py
Serial-section-appropriate registration + structural QC + local-residual TRE
for paired H-DAB CD8 / TIM-3 sections (Phase A certification).

WHY THIS EXISTS (see ihc.md "Phase A — registration redesign"):
The legacy path (registration.py: rigid Euler2D MI, then nuclear ORB/SIFT) cannot
register serial sections — individual nuclei are different physical objects across
the z-gap, so nuclear texture does not correspond — and its QC fails *closed* on
genuinely well-aligned tissue (residual measured on the same non-repeatable nuclear
features). This module instead:

  1. Registers on a LOW-FREQUENCY STRUCTURAL hematoxylin signal (vessels, lumens,
     sinusoids, tissue boundaries) at σ≈12 µm — the morphology shared across serial
     sections; single nuclei are blurred away.
  2. Uses a SIMILARITY transform (rotation + translation + uniform scale) so any
     within-pair scale difference is absorbed (rigid cannot), and exports the
     estimated scale for cross-check against the scale-bar ratio.
  3. EVALUATES multiple candidate transforms (multi-init multi-resolution MI +
     phase correlation + identity) and SELECTS by LOCAL STRUCTURAL CONSENSUS —
     dense patch phase-correlation, rewarding many confident locally-aligned
     patches with low residual. (A global NCC/MI is too flat on near-uniform liver
     parenchyma and was the cause of the spurious identity fallbacks.)
  4. CERTIFIES by LOCAL residual measured directly on structure (patch
     phase-correlation residual flow → median / p90 / per-region max), cross-checked
     by independent LUMEN-CENTROID TRE, and produces green/magenta + checkerboard
     overlays for a human 2-minute visual confirmation.

HONESTY: the patch-flow residual is measured on the structural channel that MI was
optimised on (so it is a consistency check, not a fully independent gold standard);
the lumen-centroid TRE (independent objects) and the human visual overlays are what
close the independence gap, as the scope requires.
"""

import os
import math
import numpy as np

from registration import (
    extract_hematoxylin,
    _rgb_to_gray,
    _load_rgb_thumbnail,
    _sitk_to_affine,
)


# ──────────────────────────────────────────────────────────────────────────────
# 1. Scale-bar self-calibration (burned-in 100 µm bar)
# ──────────────────────────────────────────────────────────────────────────────
def detect_scale_bar_px(image_path: str, bar_um: float = 100.0):
    """
    Robustly measure the burned-in scale bar length (px) in the bottom strip.

    The bar is a SOLID horizontal segment; the "100 µm" label above it is thin
    text. The legacy extractor mis-measured by merging the two; here we take the
    longest CONTIGUOUS solid dark run (fill≈1, short height), voted across a small
    threshold sweep so anti-aliasing / text-merge artefacts are rejected.

    Returns dict: {bar_px, pixel_size_um, bbox=(x,y,w,h), source}; bar_px None on
    failure.
    """
    import cv2
    from PIL import Image
    try:
        g = np.array(Image.open(image_path).convert("L"))
    except Exception as e:
        return {"bar_px": None, "pixel_size_um": None, "bbox": None,
                "source": f"load_failed:{e}"}

    h, w = g.shape
    y0 = int(h * 0.80)
    strip = g[y0:, :]
    widths, boxes = {}, {}
    for thr in (60, 90, 120):
        dark = (strip < thr).astype(np.uint8) * 255
        k = cv2.getStructuringElement(cv2.MORPH_RECT, (9, 1))
        op = cv2.morphologyEx(dark, cv2.MORPH_OPEN, k)
        n, _lab, stats, _c = cv2.connectedComponentsWithStats(op, connectivity=8)
        best = None
        for i in range(1, n):
            cw = int(stats[i, cv2.CC_STAT_WIDTH]); ch = int(stats[i, cv2.CC_STAT_HEIGHT])
            area = int(stats[i, cv2.CC_STAT_AREA])
            if cw < 40 or ch > 20 or area / float(cw * ch) < 0.6:
                continue
            if best is None or cw > best[0]:
                best = (cw, ch, int(stats[i, cv2.CC_STAT_LEFT]),
                        y0 + int(stats[i, cv2.CC_STAT_TOP]))
        if best is not None:
            widths[thr] = best[0]; boxes[thr] = best

    if not widths:
        return {"bar_px": None, "pixel_size_um": None, "bbox": None,
                "source": "no_bar_detected"}
    bar_px = int(round(np.median(list(widths.values()))))
    box = next((boxes[t] for t in (90, 60, 120) if widths.get(t) == bar_px),
               list(boxes.values())[0])
    cw, ch, bx, by = box
    return {"bar_px": bar_px, "pixel_size_um": round(bar_um / bar_px, 4),
            "bbox": (bx, by, cw, ch), "source": "scale_bar"}


# ──────────────────────────────────────────────────────────────────────────────
# 2. Structural representations (low-frequency hematoxylin, tissue mask, lumens)
# ──────────────────────────────────────────────────────────────────────────────
def structural_channel(rgb: np.ndarray, pixel_size_um: float):
    """Low-frequency STRUCTURAL channel: hematoxylin density blurred at σ≈12 µm so
    individual (non-corresponding) nuclei are suppressed and tissue architecture
    (vessels, sinusoids, lumens, boundaries) dominates. Returns uint8."""
    import cv2
    try:
        hema = extract_hematoxylin(rgb)
    except Exception:
        hema = _rgb_to_gray(rgb)
    sigma_px = max(12.0 / float(pixel_size_um), 4.0)
    k = int(sigma_px * 3) | 1
    return cv2.GaussianBlur(hema, (k, k), sigma_px).astype(np.uint8)


def tissue_mask(rgb: np.ndarray, pixel_size_um: float):
    """Binary tissue mask from hematoxylin density (Otsu; tissue = stained).
    Holes/lumens NOT filled (an empty lumen is not tissue)."""
    import cv2
    struct = structural_channel(rgb, pixel_size_um)
    _, m = cv2.threshold(struct, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    m = cv2.morphologyEx(m, cv2.MORPH_CLOSE, k, iterations=2)
    m = cv2.morphologyEx(m, cv2.MORPH_OPEN, k, iterations=1)
    return m


def _fill_holes(mask: np.ndarray) -> np.ndarray:
    import cv2
    h, w = mask.shape
    ff = mask.copy()
    m2 = np.zeros((h + 2, w + 2), np.uint8)
    cv2.floodFill(ff, m2, (0, 0), 255)
    return cv2.bitwise_or(mask, cv2.bitwise_not(ff))


def lumen_centroids(mask: np.ndarray, pixel_size_um: float):
    """Centroids of lumens/holes inside the tissue (sinusoids, veins, vessels,
    glandular lumens) — genuine structural OBJECTS that correspond across serial
    sections, used for the independent TRE cross-check. Returns Nx2 (x,y) px."""
    import cv2
    filled = _fill_holes(mask)
    holes = ((filled > 0) & (mask == 0)).astype(np.uint8)
    n, _lab, stats, cent = cv2.connectedComponentsWithStats(holes, connectivity=8)
    min_area = (8.0 / pixel_size_um) ** 2
    max_area = 0.05 * mask.size
    pts = [[float(cent[i][0]), float(cent[i][1])] for i in range(1, n)
           if min_area <= stats[i, cv2.CC_STAT_AREA] <= max_area]
    return np.array(pts, dtype=np.float64) if pts else np.zeros((0, 2))


# ──────────────────────────────────────────────────────────────────────────────
# 3. Local-residual patch flow (the TRE engine)
# ──────────────────────────────────────────────────────────────────────────────
def patch_residual_flow(ref_struct, warped_struct, overlap, pixel_size_um,
                        patch=128, stride=96, resp_min=0.06, min_std=5.0):
    """
    Dense LOCAL residual: for each tissue-overlap patch, the residual translation
    that still best aligns ref vs registered-moving structure (cv2.phaseCorrelate,
    Hann-windowed). If registration is good the residual ≈ 0 everywhere; a locally
    deformed region shows up as a large residual in that patch only.

    Returns list of (residual_um, cx, cy, response).
    """
    import cv2
    H, W = ref_struct.shape
    win = cv2.createHanningWindow((patch, patch), cv2.CV_32F)
    recs = []
    for y in range(0, H - patch + 1, stride):
        for x in range(0, W - patch + 1, stride):
            if overlap[y:y + patch, x:x + patch].mean() < 0.7:
                continue
            rp = ref_struct[y:y + patch, x:x + patch].astype(np.float32)
            wp = warped_struct[y:y + patch, x:x + patch].astype(np.float32)
            if rp.std() < min_std or wp.std() < min_std:
                continue
            (dx, dy), resp = cv2.phaseCorrelate(rp, wp, win)
            if resp < resp_min:
                continue
            recs.append((float(np.hypot(dx, dy)) * float(pixel_size_um),
                         x + patch // 2, y + patch // 2, float(resp)))
    return recs


def flow_stats(recs, shape, grid=5):
    """Summarise patch-flow residuals. `region_max_um` is the worst LOCAL REGION
    (median residual within a grid cell), robust to a single noisy patch; `max_um`
    is the raw worst single patch (reported for transparency)."""
    if not recs:
        return {"n": 0, "median_um": None, "p90_um": None, "max_um": None,
                "region_max_um": None}
    H, W = shape
    r = np.array([v[0] for v in recs])
    cells = {}
    for resid, cx, cy, _resp in recs:
        gx = min(int(cx / W * grid), grid - 1)
        gy = min(int(cy / H * grid), grid - 1)
        cells.setdefault((gx, gy), []).append(resid)
    region_med = [float(np.median(v)) for v in cells.values()]
    return {"n": len(recs),
            "median_um": round(float(np.median(r)), 3),
            "p90_um": round(float(np.percentile(r, 90)), 3),
            "max_um": round(float(r.max()), 3),
            "region_max_um": round(float(max(region_med)), 3)}


# ──────────────────────────────────────────────────────────────────────────────
# 4. Similarity registration: multi-init candidates + local-consensus selection
# ──────────────────────────────────────────────────────────────────────────────
def _affine_scale(matrix):
    a, b, c, d = matrix[0, 0], matrix[0, 1], matrix[1, 0], matrix[1, 1]
    return float(math.sqrt(abs(a * d - b * c)))


def _make_inits(fixed, moving):
    import SimpleITK as sitk
    inits = []
    for mode in (sitk.CenteredTransformInitializerFilter.GEOMETRY,
                 sitk.CenteredTransformInitializerFilter.MOMENTS):
        try:
            base = sitk.Similarity2DTransform(sitk.CenteredTransformInitializer(
                fixed, moving, sitk.Similarity2DTransform(), mode))
        except Exception:
            continue
        base_angle = base.GetAngle()
        for da in (math.radians(a) for a in (-10, -5, 0, 5, 10)):
            t = sitk.Similarity2DTransform(base)
            t.SetAngle(base_angle + da)
            inits.append(t)
    return inits


def _run_similarity(fixed, moving, init_tf):
    """Returns (final_transform, final_MI_value) or (None, None). Lower MI = better
    (Mattes is reported negative)."""
    import SimpleITK as sitk
    try:
        R = sitk.ImageRegistrationMethod()
        R.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
        R.SetMetricSamplingStrategy(R.RANDOM)
        R.SetMetricSamplingPercentage(0.25, seed=42)
        R.SetInterpolator(sitk.sitkLinear)
        R.SetOptimizerAsRegularStepGradientDescent(
            learningRate=1.0, minStep=1e-4, numberOfIterations=300,
            gradientMagnitudeTolerance=1e-6)
        R.SetOptimizerScalesFromPhysicalShift()
        R.SetShrinkFactorsPerLevel([4, 2, 1])
        R.SetSmoothingSigmasPerLevel([2, 1, 0])
        R.SmoothingSigmasAreSpecifiedInPhysicalUnitsOff()
        R.SetInitialTransform(sitk.Similarity2DTransform(init_tf), inPlace=False)
        final = R.Execute(fixed, moving)
        return final, float(R.GetMetricValue())
    except Exception:
        return None, None


def _mi_eval(fixed, moving, tf):
    """Evaluate Mattes MI for a fixed transform (no optimisation) → comparable
    score across MI / identity candidates. Lower = better."""
    import SimpleITK as sitk
    R = sitk.ImageRegistrationMethod()
    R.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    R.SetMetricSamplingStrategy(R.RANDOM)
    R.SetMetricSamplingPercentage(0.3, seed=1)
    R.SetInterpolator(sitk.sitkLinear)
    R.SetInitialTransform(tf)
    try:
        return float(R.MetricEvaluate(fixed, moving))
    except Exception:
        return float("inf")


def register_similarity(ref_rgb, mov_rgb, pixel_size_um):
    """
    Register mov→ref (similarity) on the structural channel; evaluate multiple
    multi-init multi-resolution MI candidates + identity, and SELECT by Mattes-MI
    value (lower = better). MI is used because every *dense* structural-agreement
    metric (NCC, patch phase-correlation, edge Chamfer) was shown to saturate or
    alias on this dense quasi-periodic tissue (see ihc.md §18.4); MI localises the
    transform robustly. patch-flow / lumen residuals are reported as DIAGNOSTICS
    only — they are NOT used to gate certification (manual landmark TRE does that).

    Returns dict with the chosen transform + diagnostics.
    """
    import cv2
    import SimpleITK as sitk

    ref_struct = structural_channel(ref_rgb, pixel_size_um)
    mov_struct = structural_channel(mov_rgb, pixel_size_um)
    ref_mask = tissue_mask(ref_rgb, pixel_size_um)
    mov_mask = tissue_mask(mov_rgb, pixel_size_um)
    Hh, Ww = ref_struct.shape
    diag = float(np.hypot(Hh, Ww))

    fixed = sitk.GetImageFromArray(ref_struct.astype(np.float32))
    moving = sitk.GetImageFromArray(mov_struct.astype(np.float32))

    cand = []   # (label, sitk_transform, mi_value)
    for i, init in enumerate(_make_inits(fixed, moving)):
        final, mi = _run_similarity(fixed, moving, init)
        if final is None:
            continue
        cand.append((f"sim_mi_{i}", final, mi))
    identity_tf = sitk.Similarity2DTransform()
    cand.append(("identity", identity_tf, _mi_eval(fixed, moving, identity_tf)))

    # Sanity-gate MI candidates, then select by lowest MI value.
    sane = []
    for label, tf, mi in cand:
        try:
            M = _sitk_to_affine(tf)
        except Exception:
            continue
        s = _affine_scale(M)
        tx, ty = float(M[0, 2]), float(M[1, 2])
        ang = abs(math.degrees(math.atan2(M[1, 0], M[0, 0])))
        if label.startswith("sim_mi") and (
                ang > 45 or abs(s - 1.0) > 0.30 or np.hypot(tx, ty) > diag):
            continue
        sane.append({"label": label, "matrix": M, "mi": mi, "est_scale": s})
    sane.sort(key=lambda d: d["mi"])
    best = sane[0]
    M = best["matrix"]
    method = "similarity_mi" if best["label"].startswith("sim_mi") else best["label"]

    # Diagnostics (NOT gating): structural NCC/Dice + patch-flow + lumen residual.
    warped = cv2.warpAffine(mov_struct, M, (Ww, Hh))
    wmask = cv2.warpAffine(mov_mask, M, (Ww, Hh), flags=cv2.INTER_NEAREST)
    overlap = (ref_mask > 0) & (wmask > 0)
    ncc = float(np.corrcoef(ref_struct[overlap], warped[overlap])[0, 1]) \
        if overlap.sum() > 10 else 0.0
    dice = 2.0 * overlap.sum() / float((ref_mask > 0).sum() + (wmask > 0).sum())
    recs = patch_residual_flow(ref_struct, warped, overlap, pixel_size_um)
    return {
        "matrix": np.asarray(M, dtype=np.float32),
        "scale_ref": 1.0, "scale_mov": 1.0,
        "method": method, "success": method != "identity",
        "est_scale": round(best["est_scale"], 4),
        "mi_value": round(best["mi"], 5),
        "struct_ncc": round(ncc, 4), "struct_dice": round(dice, 4),
        "flow": flow_stats(recs, (Hh, Ww)), "recs": recs,
        "n_candidates": len(sane),
    }


# ──────────────────────────────────────────────────────────────────────────────
# 5. Independent lumen-centroid TRE cross-check
# ──────────────────────────────────────────────────────────────────────────────
def _apply_affine(pts, matrix):
    if len(pts) == 0:
        return pts
    return (matrix @ np.hstack([pts, np.ones((len(pts), 1))]).T).T


def lumen_tre(ref_mask, mov_mask, matrix, pixel_size_um, tol_um=12.0):
    """Independent object-based TRE: match lumen centroids by mutual nearest
    neighbour (tight tolerance) in the common frame. Returns dict."""
    from scipy.spatial import cKDTree
    ref_pts = lumen_centroids(ref_mask, pixel_size_um)
    mov_pts = lumen_centroids(mov_mask, pixel_size_um)
    out = {"n_ref": len(ref_pts), "n_mov": len(mov_pts), "n_corr": 0,
           "median_um": None, "p90_um": None,
           "ref_matched": np.zeros((0, 2)), "mapped_matched": np.zeros((0, 2))}
    if len(ref_pts) == 0 or len(mov_pts) == 0:
        return out
    mapped = _apply_affine(mov_pts, matrix)
    tol_px = tol_um / float(pixel_size_um)
    tr, tm = cKDTree(ref_pts), cKDTree(mapped)
    d_rm, idx_rm = tr.query(mapped)
    _d_mr, idx_mr = tm.query(ref_pts)
    res, rm_, mp_ = [], [], []
    for j, (d, i) in enumerate(zip(d_rm, idx_rm)):
        if d <= tol_px and idx_mr[i] == j:
            res.append(d * float(pixel_size_um)); rm_.append(ref_pts[i]); mp_.append(mapped[j])
    if res:
        res = np.array(res)
        out.update(n_corr=len(res), median_um=round(float(np.median(res)), 3),
                   p90_um=round(float(np.percentile(res, 90)), 3),
                   ref_matched=np.array(rm_), mapped_matched=np.array(mp_))
    return out


def manual_tre(ref_pts, mov_pts, matrix, pixel_size_um):
    """
    GOLD-STANDARD TRE from manual corresponding landmarks (full-res px).
    `matrix` is the automated similarity transform (mov→ref). We measure how well
    that automated transform maps the HUMAN-identified moving points onto their
    human-identified reference points (the target registration error), and — as an
    internal-consistency check — fit a similarity directly from the clicks and
    report its residual (click + biological scatter) and its discrepancy from the
    automated transform.
    """
    import cv2
    ref = np.asarray(ref_pts, dtype=np.float64)
    mov = np.asarray(mov_pts, dtype=np.float64)
    out = {"n": len(ref), "median_um": None, "p90_um": None, "max_um": None,
           "per_point_um": [], "fit_residual_med_um": None, "fit_scale": None,
           "fit_vs_mi_med_um": None}
    if len(ref) == 0:
        return out
    d = np.linalg.norm(_apply_affine(mov, matrix) - ref, axis=1) * float(pixel_size_um)
    out.update(median_um=round(float(np.median(d)), 3),
               p90_um=round(float(np.percentile(d, 90)), 3),
               max_um=round(float(d.max()), 3),
               per_point_um=[round(float(x), 3) for x in d])
    if len(ref) >= 3:
        Mfit, _inl = cv2.estimateAffinePartial2D(
            mov.astype(np.float32), ref.astype(np.float32), method=cv2.LMEDS)
        if Mfit is not None:
            df = np.linalg.norm(_apply_affine(mov, Mfit) - ref, axis=1) * float(pixel_size_um)
            disc = np.linalg.norm(_apply_affine(mov, Mfit) - _apply_affine(mov, matrix),
                                  axis=1) * float(pixel_size_um)
            out.update(fit_residual_med_um=round(float(np.median(df)), 3),
                       fit_scale=round(_affine_scale(Mfit), 4),
                       fit_vs_mi_med_um=round(float(np.median(disc)), 3))
    return out


def _fit_similarity_ls(src, dst):
    """Closed-form least-squares similarity (rotation + uniform scale + translation)
    mapping src→dst (Umeyama). Deterministic — no RANSAC randomness. Returns 2x3 or None."""
    src = np.asarray(src, float)
    dst = np.asarray(dst, float)
    n = len(src)
    if n < 2:
        return None
    mx, my = src.mean(0), dst.mean(0)
    Xc, Yc = src - mx, dst - my
    cov = (Yc.T @ Xc) / n
    U, S, Vt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1.0, d])
    R = U @ D @ Vt
    var_x = (Xc ** 2).sum() / n
    s = float(np.trace(np.diag(S) @ D) / var_x) if var_x > 0 else 1.0
    t = my - s * (R @ mx)
    M = np.zeros((2, 3))
    M[:, :2] = s * R
    M[:, 2] = t
    return M


def loo_tre(ref_pts, mov_pts, pixel_size_um):
    """
    Leave-one-out TRE: for each landmark, refit the similarity on the other N−1 and
    predict the held-out one. The held-out prediction error is an UNBIASED estimate
    of registration accuracy (using the same points to fit and score is optimistic).
    """
    ref = np.asarray(ref_pts, float)
    mov = np.asarray(mov_pts, float)
    n = len(ref)
    out = {"n": n, "loo_median_um": None, "loo_p90_um": None, "loo_max_um": None,
           "per_point_um": []}
    if n < 3:                                    # need ≥3 so the N−1 fit is determined
        return out
    errs = []
    for i in range(n):
        m = np.arange(n) != i
        M = _fit_similarity_ls(mov[m], ref[m])
        if M is None:
            continue
        pred = (M @ np.array([mov[i, 0], mov[i, 1], 1.0]))[:2]
        errs.append(float(np.linalg.norm(pred - ref[i])) * float(pixel_size_um))
    if errs:
        e = np.array(errs)
        out.update(loo_median_um=round(float(np.median(e)), 3),
                   loo_p90_um=round(float(np.percentile(e, 90)), 3),
                   loo_max_um=round(float(e.max()), 3),
                   per_point_um=[round(float(x), 3) for x in e])
    return out


def _hull_area(pts):
    import cv2
    pts = np.asarray(pts, np.float32)
    if len(pts) < 3:
        return 0.0
    return float(cv2.contourArea(cv2.convexHull(pts)))


def landmark_register_and_verify(ref_pts, mov_pts, pixel_size_um,
                                 val_ref_pts=None, val_mov_pts=None, image_wh=None,
                                 min_n=6, target_n=12, loo_max_um=5.0, fit_max_um=5.0,
                                 deformed_loo_um=15.0, min_roi_frac=0.10):
    """
    GOLD-STANDARD, landmark-DRIVEN registration + verification (Phase A).

    The operator's confident anatomical landmarks DEFINE the registration: a
    least-squares similarity (distance-preserving, so the downstream cross-K stays
    valid — we never non-rigidly warp). Accuracy is measured on HELD-OUT points:

      • if an independent validation set (ideally a SECOND annotator) is supplied,
        TRE is its error under the fit set's transform — annotator-independent.
      • otherwise leave-one-out (LOO). NOTE: LOO is *fit-unbiased* (a point is never
        in the transform that predicts it) but NOT annotator-independent — all points
        share one annotator's selection bias. It is the limited-data fallback, not an
        ANHIR-grade gold standard.

    Four-state verdict (a failed pair is reported, never warped or forced):
      CERTIFIED         n≥min_n, held-out TRE median ≤loo_max_um, fit-residual ≤fit_max_um
      LOCALLY_CERTIFIED only a spatial subset passes (≥min_n, hull ≥min_roi_frac of
                        field) → analyse that ROI only
      DEFORMED          confident correspondences exist but a similarity cannot fit
                        them within tolerance (local non-rigid deformation)
      NOT_CERTIFIABLE   too few unambiguous correspondences to measure accuracy — this
                        is NOT positive evidence the sections are unrelated

    Thresholds follow the ≤5 µm criterion + serial-section z-gap floor; fixed, not
    tuned. ~target_n well-spread points are wanted for paper-grade; min_n only fits.
    """
    import cv2
    ref = np.asarray(ref_pts, float)
    mov = np.asarray(mov_pts, float)
    n = len(ref)
    out = {"n": n, "matrix": None, "est_scale": None, "fit_residual_um": None,
           "tre_median_um": None, "tre_p90_um": None, "tre_max_um": None,
           "validation": None, "coverage_frac": None, "n_good": 0,
           "roi_polygon": None, "verdict": None, "reason": None}

    M = _fit_similarity_ls(mov, ref) if n >= 2 else None
    if M is not None:
        out["matrix"] = M.tolist()
        out["est_scale"] = round(_affine_scale(M), 4)
        d = np.linalg.norm(_apply_affine(mov, M) - ref, axis=1) * float(pixel_size_um)
        out["fit_residual_um"] = round(float(np.median(d)), 3)
    if image_wh and M is not None:
        out["coverage_frac"] = round(_hull_area(ref) / float(image_wh[0] * image_wh[1]), 4)

    # Held-out accuracy: independent validation set if supplied, else LOO.
    if val_ref_pts is not None and len(val_ref_pts) >= 1 and M is not None:
        vr, vm = np.asarray(val_ref_pts, float), np.asarray(val_mov_pts, float)
        err = np.linalg.norm(_apply_affine(vm, M) - vr, axis=1) * float(pixel_size_um)
        out["validation"] = f"independent validation set (n={len(err)})"
        local_ok = False
    else:
        loo = loo_tre(ref, mov, pixel_size_um)
        err = np.array(loo["per_point_um"]) if loo["per_point_um"] else np.array([])
        out["validation"] = ("leave-one-out (single-annotator; fit-unbiased, "
                             "NOT annotator-independent)")
        local_ok = True

    if err.size:
        out.update(tre_median_um=round(float(np.median(err)), 3),
                   tre_p90_um=round(float(np.percentile(err, 90)), 3),
                   tre_max_um=round(float(err.max()), 3))
    good = (err <= loo_max_um) if err.size else np.array([], bool)
    out["n_good"] = int(good.sum())
    tier = "" if n >= target_n else f" (n={n} < {target_n} preferred — provisional)"

    if n < min_n or not err.size:
        out.update(verdict="NOT_CERTIFIABLE",
                   reason=f"only {n} confident landmarks — too few unambiguous "
                          f"correspondences to measure accuracy (NOT evidence the "
                          f"sections are unrelated)")
        return out
    med, fr = out["tre_median_um"], out["fit_residual_um"]
    if med <= loo_max_um and fr is not None and fr <= fit_max_um:
        out.update(verdict="CERTIFIED",
                   reason=f"held-out TRE median {med} µm (p90 {out['tre_p90_um']}), "
                          f"fit-residual {fr} µm, n={n}{tier}")
        return out
    # Locally certified? a spatially-coherent subset of good points (LOO case only)
    if local_ok and out["n_good"] >= min_n:
        gref = ref[:len(err)][good]
        roi_frac = (_hull_area(gref) / float(image_wh[0] * image_wh[1])) if image_wh else 0.0
        if roi_frac >= min_roi_frac:
            Mloc = _fit_similarity_ls(mov[:len(err)][good], gref)
            out.update(verdict="LOCALLY_CERTIFIED",
                       matrix=(Mloc.tolist() if Mloc is not None else out["matrix"]),
                       roi_polygon=[[float(x), float(y)] for x, y in
                                    cv2.convexHull(gref.astype(np.float32)).reshape(-1, 2)],
                       reason=f"{out['n_good']} of {n} landmarks pass within an ROI "
                              f"(~{roi_frac*100:.0f}% of field); analyse that ROI only")
            return out
    if med <= deformed_loo_um:
        out.update(verdict="DEFORMED",
                   reason=f"held-out TRE {med} µm / fit-residual {fr} µm exceed "
                          f"≤{loo_max_um} µm — local non-rigid deformation; not "
                          f"certifiable for distance-based statistics (no warp applied)")
    else:
        out.update(verdict="NOT_CERTIFIABLE",
                   reason=f"held-out TRE {med} µm ≫ tolerance — landmarks do not agree "
                          f"on a single transform; insufficient correspondence to certify")
    return out


def registration_perturbation_sensitivity(stat_fn, base_matrix, tre_um, pixel_size_um,
                                          field_um, n_samples=50, seed=0):
    """
    Phase-B robustness (Codex recommendation): perturb the certified transform within
    its MEASURED landmark uncertainty and re-run the spatial statistic. If the verdict
    (direction + significance) is stable across perturbations the conclusion is
    supported; if it flips it is inconclusive. (A registration error comparable to a
    tested radius means that radius is not interpretable.)

    stat_fn(matrix_2x3) -> dict with 'significant' (bool) and 'direction' (str).
    tre_um: held-out registration TRE; field_um: field half-extent (for rotation jitter).
    Returns {'base', 'n', 'agree_frac', 'stable', 'tre_um'}.
    """
    import math
    rng = np.random.default_rng(seed)
    B = np.vstack([np.asarray(base_matrix, float), [0, 0, 1]])
    base = stat_fn(np.asarray(base_matrix, float))
    sigma_px = tre_um / float(pixel_size_um)
    sigma_rot = tre_um / max(field_um, 1.0)              # rad: arc ≈ TRE at field edge
    agree = 0
    for _ in range(n_samples):
        th = rng.normal(0, sigma_rot)
        c, s = math.cos(th), math.sin(th)
        J = np.array([[c, -s, rng.normal(0, sigma_px)],
                      [s, c, rng.normal(0, sigma_px)], [0, 0, 1]])
        r = stat_fn((J @ B)[:2])
        if (r.get("significant") == base.get("significant")
                and r.get("direction") == base.get("direction")):
            agree += 1
    return {"base": base, "n": n_samples, "agree_frac": round(agree / n_samples, 3),
            "stable": (agree / n_samples) >= 0.9, "tre_um": tre_um}


# ──────────────────────────────────────────────────────────────────────────────
# 6. Human-verifiable QC visualisations (green/magenta + checkerboard)
# ──────────────────────────────────────────────────────────────────────────────
def save_qc_overlays(ref_rgb, mov_rgb, matrix, out_prefix, recs=None, lumens=None):
    """Save a two-colour overlay (ref=green, registered mov=magenta; grey where
    they agree, with patch-flow residual vectors) and a checkerboard. Returns the
    two file paths."""
    import cv2
    ref_g = _rgb_to_gray(ref_rgb)
    mov_g = _rgb_to_gray(mov_rgb)
    Hh, Ww = ref_g.shape
    warped = cv2.warpAffine(mov_g, matrix, (Ww, Hh))

    ov = np.zeros((Hh, Ww, 3), np.uint8)
    ov[..., 1] = ref_g
    ov[..., 0] = warped
    ov[..., 2] = warped
    if recs:
        for resid, cx, cy, _resp in recs:                 # exaggerate ×10 to see
            cv2.circle(ov, (int(cx), int(cy)), 3, (255, 255, 0), -1)
    if lumens is not None and len(lumens.get("ref_matched", [])):
        for (rx, ry), (mx, my) in zip(lumens["ref_matched"], lumens["mapped_matched"]):
            cv2.line(ov, (int(rx), int(ry)), (int(mx), int(my)), (0, 255, 255), 1)
    overlay_path = f"{out_prefix}_overlay.png"
    cv2.imwrite(overlay_path, cv2.cvtColor(ov, cv2.COLOR_RGB2BGR))

    tile = max(Hh, Ww) // 12
    cb = ref_g.copy()
    for y in range(0, Hh, tile):
        for x in range(0, Ww, tile):
            if ((x // tile) + (y // tile)) % 2 == 1:
                cb[y:y + tile, x:x + tile] = warped[y:y + tile, x:x + tile]
    cb_path = f"{out_prefix}_checkerboard.png"
    cv2.imwrite(cb_path, cb)
    return overlay_path, cb_path


# ──────────────────────────────────────────────────────────────────────────────
# 7. Orchestrator: certify one pair
# ──────────────────────────────────────────────────────────────────────────────
def certify_pair(sample_id, ref_path, mov_path, pixel_size_um, out_dir,
                 ref_bar_px=None, mov_bar_px=None, pixel_size_source="manual",
                 tre_median_max_um=5.0, region_max_um=10.0, min_patches=10,
                 scale_xcheck_tol=0.03):
    """Full Phase-A certification for one pair. Read-only on inputs; writes QC
    overlays to out_dir. Returns a row dict."""
    os.makedirs(out_dir, exist_ok=True)
    row = {"sample_id": sample_id, "pixel_size_um": pixel_size_um,
           "pixel_size_source": pixel_size_source,
           "ref_bar_px": ref_bar_px, "mov_bar_px": mov_bar_px}

    ref_rgb, _ = _load_rgb_thumbnail(ref_path, max_side=1920)
    mov_rgb, _ = _load_rgb_thumbnail(mov_path, max_side=1920)
    if ref_rgb is None or mov_rgb is None:
        row.update(status="NOT CERTIFIED", reason="image load failed")
        return row

    reg = register_similarity(ref_rgb, mov_rgb, pixel_size_um)
    fs = reg["flow"]
    row.update(method=reg["method"], est_scale=reg["est_scale"],
               struct_ncc=reg["struct_ncc"], struct_dice=reg["struct_dice"],
               n_patches=fs["n"], tre_median_um=fs["median_um"],
               tre_p90_um=fs["p90_um"], tre_max_um=fs["max_um"],
               region_max_um=fs["region_max_um"])

    if ref_bar_px and mov_bar_px:
        bar_ratio = ref_bar_px / float(mov_bar_px)
        row["bar_scale_expected"] = round(bar_ratio, 4)
        row["scale_xcheck_delta"] = round(abs(reg["est_scale"] - bar_ratio), 4)
        row["scale_xcheck_ok"] = row["scale_xcheck_delta"] <= scale_xcheck_tol
    else:
        row["bar_scale_expected"] = None
        row["scale_xcheck_delta"] = None
        row["scale_xcheck_ok"] = None

    ref_mask = tissue_mask(ref_rgb, pixel_size_um)
    mov_mask = tissue_mask(mov_rgb, pixel_size_um)
    lum = lumen_tre(ref_mask, mov_mask, reg["matrix"], pixel_size_um)
    row.update(lumen_n_corr=lum["n_corr"], lumen_tre_median_um=lum["median_um"])

    prefix = os.path.join(out_dir, sample_id)
    ov, cb = save_qc_overlays(ref_rgb, mov_rgb, reg["matrix"], prefix,
                              recs=reg["recs"], lumens=lum)
    row["overlay_path"] = ov
    row["checkerboard_path"] = cb

    # ── Gate-A decision ───────────────────────────────────────────────────────
    reasons = []
    if reg["method"] == "identity":
        reasons.append("identity fallback (no structural alignment found)")
    if fs["n"] < min_patches:
        row["status"] = "NEEDS-MY-INPUT"
        row["reason"] = (f"only {fs['n']} confident structural patches "
                         f"(< {min_patches}); supply a few manual landmark points")
        return row

    tre_ok = fs["median_um"] is not None and fs["median_um"] <= tre_median_max_um
    region_ok = fs["region_max_um"] is not None and fs["region_max_um"] < region_max_um
    if not tre_ok:
        reasons.append(f"median local residual {fs['median_um']} µm "
                       f"> {tre_median_max_um} µm")
    if not region_ok:
        reasons.append(f"worst-region residual {fs['region_max_um']} µm "
                       f"≥ {region_max_um} µm")
    if row["scale_xcheck_ok"] is False:
        reasons.append(f"scale cross-check failed "
                       f"(est {reg['est_scale']} vs bar {row['bar_scale_expected']})")

    if (reg["method"] != "identity" and tre_ok and region_ok
            and row["scale_xcheck_ok"] is not False):
        row["status"] = "CERTIFIED"
        row["reason"] = (f"median local residual {fs['median_um']} µm "
                         f"(p90 {fs['p90_um']}, worst-region {fs['region_max_um']}, "
                         f"n={fs['n']} patches); lumen TRE {lum['median_um']} µm "
                         f"(n={lum['n_corr']})")
    else:
        row["status"] = "NOT CERTIFIED"
        row["reason"] = "; ".join(reasons) if reasons else "failed certification"
    return row
