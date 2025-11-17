import re
import os
import io
import base64
import skimage
import warnings
import numpy as np
import pandas as pd
from PIL import Image
from pathlib import Path
from itertools import product
from cellpose import models, transforms
from skimage.morphology import dilation, disk
from skimage.segmentation import find_boundaries
warnings.filterwarnings("ignore", message=".*weights_only=False.*")
warnings.filterwarnings("ignore", message=".*low contrast image.*")


def expand_range_token(tok):
    if '-' in tok:
        start, end = tok.split('-', 1)
        width = max(len(start), len(end))
        return [f"{i:0{width}d}" for i in range(int(start), int(end) + 1)]
    return [tok]

_TILE_RE = re.compile(r"^L(?P<L>\d+)R(?P<R>[\d\-]+)C(?P<C>[\d\-]+)S(?P<S>\d+)$")

def expand_tile_spec(spec):
    m = _TILE_RE.match(spec)
    if not m:
        raise ValueError(f"Invalid tile spec format: {spec}")
    L = m.group("L")
    R_list = expand_range_token(m.group("R"))
    C_list = expand_range_token(m.group("C"))
    S = m.group("S")
    return [f"L{L}R{r}C{c}S{S}" for r, c in product(R_list, C_list)]

# Well definitions and tile mapping
TILE_MAPS_SPEC = {
    "TwelveWellStandard": {
        "A1": "L2R01-03C01-04S1", "A2": "L1R01-03C01-04S1",
        "B1": "L2R04-06C01-04S1", "B2": "L1R04-06C01-04S1",
        "C1": "L2R07-09C01-04S1", "C2": "L1R07-09C01-04S1",
        "D1": "L2R10-12C01-04S1", "D2": "L1R10-12C01-04S1",
        "E1": "L2R13-15C01-04S1", "E2": "L1R13-15C01-04S1",
        "F1": "L2R16-18C01-04S1", "F2": "L1R16-18C01-04S1"
    },
    "OneWell": {"A1": "L1R01-21C01-11S1"},
    "FortyEightWell": {
        f"{row}{col}": f"L1R{str(i).zfill(2)}C{str((col-1)*2+1).zfill(2)}-{str(col*2).zfill(2)}S1"
        for i, row in enumerate("ABCDEFGHIJKL", 1)
        for col in range(1, 5)
    }
}

WELL_LAYOUT_MAP = {
    "OneWell": ["A1"],
    "TwelveWellStandard": ["A1","A2","B1","B2","C1","C2","D1","D2","E1","E2","F1","F2"],
    "FortyEightWell": [f"{r}{c}" for r in "ABCDEFGHIJKL" for c in range(1, 5)],
}

TILE_MAPS = {layout: {w: expand_tile_spec(spec) for w, spec in d.items()} for layout, d in TILE_MAPS_SPEC.items()}


OVERLAY_MIN_CROP_SIZE = 256
OVERLAY_MAX_CROP_SIZE = 512
OVERLAY_CROP_MARGIN = 40


def is_two_channel_model(model_name: str) -> bool:
    return model_name.lower().endswith("2ch")


def is_three_channel_model(model_name: str) -> bool:
    return model_name.lower().endswith("3ch")


def is_general_model(model_name: str) -> bool:
    return "general" in model_name.lower()


def _compute_crop_bounds(center: int, size: int, max_len: int):
    half = size // 2
    start = max(center - half, 0)
    end = min(start + size, max_len)
    # Adjust start in case we hit the boundary
    start = max(0, end - size)
    return start, end

def _find_dense_crop(cell_presence: np.ndarray, crop_h: int, crop_w: int):
    """Return top-left coordinates for the crop window that captures the most cells."""
    h, w = cell_presence.shape
    if crop_h >= h and crop_w >= w:
        return 0, 0

    prefix = np.pad(cell_presence.astype(np.int32), ((1, 0), (1, 0)), mode="constant")
    prefix = prefix.cumsum(axis=0).cumsum(axis=1)

    best_sum = -1
    best_top = 0
    best_left = 0

    for top in range(0, h - crop_h + 1):
        bottom = top + crop_h
        window_sums = (
            prefix[bottom, crop_w:] - prefix[bottom, : w - crop_w + 1]
            - prefix[top, crop_w:] + prefix[top, : w - crop_w + 1]
        )
        if window_sums.size == 0:
            continue
        col_idx = int(np.argmax(window_sums))
        current_sum = int(window_sums[col_idx])
        if current_sum > best_sum:
            best_sum = current_sum
            best_top = top
            best_left = col_idx

    return best_top, best_left

def overlay_thumbnail(cell_mask, binary_nuclei, cell_image=None, nuclei_color=(1, 0, 1)):
    """
    Create an RGB overlay of cells, nuclei, and background image.

    Args:
        cell_mask (ndarray): integer mask of cells.
        binary_nuclei (ndarray): binary mask of nuclei.
        cell_image (ndarray, optional): HxWxC image; all channels are normalized and composited.
        nuclei_color (tuple): RGB tuple for nuclei overlay (default = magenta).

    Returns:
        ndarray: RGB image in float32 [0, 1].
    """
    cell_presence = cell_mask > 0
    h, w = cell_mask.shape

    # --- Find crop around dense cell region ---
    if cell_presence.any():
        coords = np.argwhere(cell_presence)
        min_r, min_c = coords.min(axis=0)
        max_r, max_c = coords.max(axis=0)
        bbox_h = max_r - min_r + 1
        bbox_w = max_c - min_c + 1

        crop_h = min(OVERLAY_MAX_CROP_SIZE, max(bbox_h + 2 * OVERLAY_CROP_MARGIN, OVERLAY_MIN_CROP_SIZE, 1))
        crop_w = min(OVERLAY_MAX_CROP_SIZE, max(bbox_w + 2 * OVERLAY_CROP_MARGIN, OVERLAY_MIN_CROP_SIZE, 1))
        crop_h = min(crop_h, h)
        crop_w = min(crop_w, w)

        top, left = _find_dense_crop(cell_presence, crop_h, crop_w)
        r_start, r_end = top, top + crop_h
        c_start, c_end = left, left + crop_w
    else:
        crop_h = min(OVERLAY_MAX_CROP_SIZE, h)
        crop_w = min(OVERLAY_MAX_CROP_SIZE, w)
        center_r = h // 2
        center_c = w // 2
        r_start, r_end = _compute_crop_bounds(center_r, crop_h, h)
        c_start, c_end = _compute_crop_bounds(center_c, crop_w, w)

    # --- Crop data ---
    cell_mask = cell_mask[r_start:r_end, c_start:c_end]
    binary_nuclei = binary_nuclei[r_start:r_end, c_start:c_end]

    # --- Normalize all channels in cell_image ---
    if cell_image is not None:
        cell_image = np.asarray(cell_image, dtype=np.float32)
        if cell_image.ndim == 2:
            cell_image = cell_image[..., np.newaxis]
        cell_background = cell_image[r_start:r_end, c_start:c_end, :]

        # Robust percentile normalization per channel
        norm_channels = []
        for ch in range(cell_background.shape[2]):
            channel = cell_background[..., ch]
            v_min, v_max = np.percentile(channel, (2, 98))
            if v_max <= v_min:
                v_min, v_max = channel.min(), channel.max()
            if v_max > v_min:
                norm = np.clip((channel - v_min) / (v_max - v_min), 0.0, 1.0)
            else:
                norm = np.zeros_like(channel)
            norm_channels.append(norm)

        norm_stack = np.stack(norm_channels, axis=-1)
        if norm_stack.shape[2] > 3:
            # Take first 3 channels (assumes biologically relevant order)
            rgb = norm_stack[..., :3]
        elif norm_stack.shape[2] == 3:
            rgb = norm_stack
        else:
            repeats = -(-3 // norm_stack.shape[2])
            rgb = np.repeat(norm_stack, repeats, axis=2)[..., :3]
    else:
        rgb = np.ones((r_end - r_start, c_end - c_start, 3), dtype=np.float32)

    # --- Overlays ---
    from skimage.morphology import disk, dilation
    from skimage.segmentation import find_boundaries

    cell_borders = find_boundaries(cell_mask, mode="inner")
    thick_borders = dilation(cell_borders.astype(bool), disk(1))

    # Cell borders overlay (white)
    if np.any(thick_borders):
        white = np.array([1, 1, 1], dtype=float)
        rgb[thick_borders] = 0.4 * rgb[thick_borders] + 0.6 * white

    return np.clip(rgb, 0, 1)

def image_to_thumbnail_html(rgb, max_size=256):
    """Convert uint8 RGB numpy array to base64 PNG <img/> HTML string for DataFrame embedding."""
    if rgb.dtype!= np.uint8:
        rgb = (np.clip(rgb, 0, 1) * 255).astype(np.uint8)
        
    img = Image.fromarray(rgb, mode='RGB')
    
    img.thumbnail((max_size, max_size))
    with io.BytesIO() as buf:
        img.save(buf, format="PNG",optimize=True)
        #b64 = base64.b64encode(buf.getvalue()).decode()
        encoded = base64.b64encode(buf.getvalue()).decode("utf-8")
        html = (
            f'<div style="background-color:white;display:flex;justify-content:center;'
            f'align-items:center;width:{max_size}px;height:{max_size}px;'
            f'border:1px solid #ddd;border-radius:8px;box-shadow:0 0 3px rgba(0,0,0,0.2);">'
            f'<img src="data:image/png;base64,{encoded}" style="max-width:100%;max-height:100%;"/></div>'
        )
    return html

#Function to Normalize images from the Projection folder
def normalize_image(image, region_size=1824):
    """
    Normalize an image within a defined region size.
    """
    image_norm = np.zeros_like(image, np.single)

    for xi in range(int(image.shape[1] / region_size)):
        for yi in range(int(image.shape[0] / region_size)):
            cropped = image[
                            yi * region_size: (yi + 1) * region_size, xi * region_size: (xi + 1) * region_size
                        ]
            cropped = transforms.normalize_img(cropped.reshape(
                            cropped.shape[0], cropped.shape[1], 1)).reshape(cropped.shape[0], cropped.shape[1])
            image_norm[
                            yi * region_size: (yi + 1) * region_size, xi * region_size: (xi + 1) * region_size
                        ] = cropped

    return image_norm

def segment_cells(cell_image, nuclear_image, cell_model_path, nuclear_model_path, cell_diameter=None, normalization=True, use_gpu=False, actin_image=None):
    """Run segmentation with given models + diameter, return metrics + overlay thumbnail."""
    cell_image_original = np.asarray(cell_image)
    nuclear_image_original = np.asarray(nuclear_image)
    actin_image_original = np.asarray(actin_image) if actin_image is not None else None
    composite_original = None

    def maybe_normalize(img, do_norm):
        return normalize_image(img) if do_norm else np.asarray(img, dtype=np.float32)

    # Normalize if requested
    cell_image = maybe_normalize(cell_image_original, normalization)
    nuclear_image = maybe_normalize(nuclear_image_original, normalization)
    if actin_image_original is not None:
        actin_image = maybe_normalize(actin_image_original, normalization)

    model_name = Path(cell_model_path).name
    three_channel_model = is_three_channel_model(model_name)
    general_model = is_general_model(model_name)
    if three_channel_model and actin_image is None:
        raise ValueError(f"Actin image is required for 3-channel model: {cell_model_path}")
    if three_channel_model and actin_image.shape != cell_image.shape:
        raise ValueError(
            f"Actin image shape {actin_image.shape} does not match cell image shape {cell_image.shape} "
            f"for model {cell_model_path}"
        )

    if three_channel_model:
        composite = np.stack([cell_image, nuclear_image, actin_image], axis=-1)
        composite_original = np.stack([cell_image_original, nuclear_image_original, actin_image_original], axis=-1)
        nchan = 3
    else:
        composite = np.stack([cell_image, nuclear_image], axis=-1)
        composite_original = np.stack([cell_image_original, nuclear_image_original], axis=-1)
        nchan = 2

    # --- Cell segmentation ---
    
    cell_model = models.CellposeModel(gpu=use_gpu, pretrained_model=False, model_type=cell_model_path, nchan=nchan)
    eval_kwargs = dict(channels=None, normalize=False, resample=False)
    if general_model:
        if cell_diameter is None:
            raise ValueError(f"Cell diameter is required for general model: {cell_model_path}")
        eval_kwargs["diameter"] = float(cell_diameter)
    cell_mask, _, _ = cell_model.eval(composite, **eval_kwargs)
    cell_mask = cell_mask.astype(np.uint16)

    # --- Nuclear segmentation ---
    nuc_model = models.CellposeModel(gpu=use_gpu, pretrained_model=False, model_type=nuclear_model_path)
    nuclear_mask, _, _ = nuc_model.eval(nuclear_image, resample=False)
    binary_nuclei = nuclear_mask.copy()
    binary_nuclei[nuclear_mask > 0] = 1

    # --- Metrics ---
    num_cells = len(np.unique(cell_mask)) - 1
    num_nuclei = len(np.unique(nuclear_mask)) - 1

    # Relationship checks
    #nuc_labels, _ = label(nuclear_mask > 0)
    cells_without_nuc = 0
    nuclei_outside_cells = 0
    
    for c in np.unique(cell_mask):
        if c == 0: 
            continue
        mask_c = (cell_mask == c)
        if not np.any(np.logical_and(mask_c, nuclear_mask > 0)):
            cells_without_nuc += 1

    for n in np.unique(nuclear_mask):
        if n == 0: 
            continue
        mask_n = (nuclear_mask == n)
        if not np.any(np.logical_and(mask_n, cell_mask > 0)):
            nuclei_outside_cells += 1

    pct_cells_wo_nuc = 100 * cells_without_nuc / max(num_cells, 1)
    pct_nuc_outside = 100 * nuclei_outside_cells / max(num_nuclei, 1)

    binary_cell = cell_mask > 0

    # --- Overlay thumbnail ---
    thumb = overlay_thumbnail(cell_mask, binary_nuclei, composite_original)
    thumb_html = image_to_thumbnail_html(thumb, max_size=256)
    metrics = {
        "Number of Cells": num_cells,
        "Number of Nuclei": num_nuclei,
        "Mean Cell Area": f"{binary_cell.sum() / max(num_cells, 1):.2f}",
        "Percent Cells without Nucleus": f"{pct_cells_wo_nuc:.2f}",
        "Percent Nuclei outside cells": f"{pct_nuc_outside:.2f}",
        "Overlay": thumb_html
    } 

    return  cell_mask, binary_nuclei, metrics

def cropImg(im, r):
    h, w = im.shape

    crop_img = im[int(h/2 - r):int(h/2 + r), int(w/2 - r):int(w/2 + r)]

    return crop_img
def run_segmentation(selection_well_tiles,selected_models,model_dir, run_path, output_location,cell_diameters=None):
    # Extract selection info
    layout = selection_well_tiles.get("layout")  # stored in the widget output
    selected_wells = selection_well_tiles.get("wells", [])
    selected_tiles = selection_well_tiles.get("tiles", [])
    cell_models = selected_models.get("cell_models", [])
    nuclear_model = selected_models.get("nuclear_model", [])
    selection_diameters = selection_well_tiles.get("cell_diameters")
    optimize_crop = bool(selection_well_tiles.get("optimize_crop", True))

    if selection_diameters is not None:
        provided_diameters = selection_diameters
    else:
        provided_diameters = cell_diameters

    if provided_diameters is None:
        diameter_values = []
    elif isinstance(provided_diameters, (list, tuple)):
        diameter_values = [float(d) for d in provided_diameters if d is not None]
    else:
        diameter_values = [float(provided_diameters)]

    # Prepare results container
    results = []

    # Build a proper mapping: well -> list of tiles (filtered if user chose specific tiles)
    well_tile_map = {}
    for well in selected_wells:
        # Get all possible tiles for this well from TILE_MAPS
        all_tiles = TILE_MAPS[layout].get(well, [])
        
        # If user manually selected tiles, intersect them
        tiles = [t for t in all_tiles if not selected_tiles or t in selected_tiles]
        
        if tiles:
            well_tile_map[well] = tiles

    # Main segmentation loop
    for well, tiles in well_tile_map.items():
        for tile in tiles:

            # Load images from the run folders
            cell_path = os.path.join(run_path, "Projection", f"Well{well}", f"CP01_{tile}_Cell-Membrane.tif")
            nuc_path  = os.path.join(run_path, "Projection", f"Well{well}", f"CP01_{tile}_Nucleus.tif")
            actin_path = os.path.join(run_path, "Projection", f"Well{well}", f"CP01_{tile}_Actin.tif")
        
            # Check file existence
            if not os.path.exists(cell_path):
                print(f"[WARN] Missing cell image: {cell_path} - skipping tile {tile} in well {well}")
                continue
            if not os.path.exists(nuc_path):
                print(f"[WARN] Missing nuclear image: {nuc_path} - skipping tile {tile} in well {well}")
                continue

            # Read images
            try:
                cell_image = skimage.io.imread(cell_path)
            except OSError as exc:
                print(f"[WARN] Failed to read cell image {cell_path}: {exc} - skipping tile {tile} in well {well}")
                continue

            try:
                nuclear_image = skimage.io.imread(nuc_path)
            except OSError as exc:
                print(f"[WARN] Failed to read nuclear image {nuc_path}: {exc} - skipping tile {tile} in well {well}")
                continue
            actin_image_cache = None

            # Run segmentation for each cell model
            for cell_model in cell_models:
                cell_model_path = os.path.join(model_dir, cell_model)
                nuclear_model_path = os.path.join(model_dir, nuclear_model)
                print(f"Running segmentation for {well} - {tile} using {cell_model} / {nuclear_model}")

                actin_image = None
                if is_three_channel_model(cell_model):
                    if not os.path.exists(actin_path):
                        print(f"Missing actin image: {actin_path} -> skipping model {cell_model}")
                        continue
                    if actin_image_cache is None:
                        try:
                            actin_image_cache = skimage.io.imread(actin_path)
                        except OSError as exc:
                            print(f"[WARN] Failed to read actin image {actin_path}: {exc} -> skipping model {cell_model}")
                            continue
                    actin_image = actin_image_cache
                
                is_general = is_general_model(cell_model)
                diameter_options = diameter_values if is_general else [None]
                if is_general and not diameter_options:
                    raise ValueError(
                        f"General model {cell_model} requires at least one diameter. "
                        "Please provide comma-separated diameters in the selection step."
                    )

                cropped_cell = cropImg(cell_image,912) if optimize_crop else cell_image
                cropped_nuc = cropImg(nuclear_image,912) if optimize_crop else nuclear_image
                cropped_actin = cropImg(actin_image,912) if (optimize_crop and actin_image is not None) else actin_image

                for diameter in diameter_options:
                    cell_input = cropped_cell
                    nuclear_input = cropped_nuc
                    actin_input = cropped_actin

                    # Run segmentation
                    cell_mask, binary_nuclei, seg_result = segment_cells(
                        cell_input,
                        nuclear_input,
                        cell_model_path,
                        nuclear_model_path,
                        cell_diameter=diameter,
                        actin_image=actin_input
                    )

                    # Save outputs
                    out_dir = os.path.join(output_location, cell_model, f"Well{well}")
                    os.makedirs(out_dir, exist_ok=True)

                    cell_outfile = os.path.join(out_dir, f"CP01_{tile}_Cell-Membrane.tif")
                    nuc_outfile  = os.path.join(out_dir, f"CP01_{tile}_Nucleus.tif")

                    skimage.io.imsave(cell_outfile, cell_mask.astype("uint16"))
                    skimage.io.imsave(nuc_outfile, binary_nuclei.astype("uint8"))

                    # Record result
                    seg_entry = {
                        "Cell Membrane Model": cell_model,
                        "Nuclear Model": nuclear_model,
                        "Well": well,
                        "Tile": tile,
                        "Diameter (px)": f"{diameter:.2f}" if diameter is not None else "",
                        **seg_result
                    }
                    results.append(seg_entry)

    # Build dataframe for visualization
    results_df = pd.DataFrame(results)
    return results_df
