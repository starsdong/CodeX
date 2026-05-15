#!/usr/bin/env python3
"""Merge selected sPHENIX pointing-resolution curves and data points.

Inputs:
- two curves from pixel_intt_single_track_pointing_resolution.png
- red data points from sPHENIX_simu.png
- black iterative-Gaussian points from the local CSV matching the ROOT plot
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
from PIL import Image, ImageDraw, ImageFont


OUT_DIR = Path(__file__).resolve().parent
FIRST_FIGURE = Path("/Users/starsdong/Work/CodeX/FastSim/pixel_intt_single_track_pointing_resolution_ymax002.png")
SECOND_FIGURE = Path("/Users/starsdong/Desktop/sPHENIX_simu.png")
THIRD_CSV = OUT_DIR / "resolution_vs_x_methods_final/output_PlottingMacro_data/hdcaxy_vs_preco_allparticles_resolution_comparison_vs_x.csv"

FONT_REGULAR = "/System/Library/Fonts/Supplemental/Arial.ttf"

WIDTH, HEIGHT = 1800, 1800
LEFT, TOP, RIGHT, BOTTOM = 300, 110, 1740, 1550
X_MIN, X_MAX = 0.0, 5.5
Y_MIN, Y_MAX = 0.0, 0.010

BLUE = (0, 114, 178)
GRAY = (86, 86, 86)
RED = (255, 0, 0)
BLACK = (0, 0, 0)


def load_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    try:
        return ImageFont.truetype(FONT_REGULAR, size)
    except OSError:
        return ImageFont.load_default()


def x_to_px(x: float) -> float:
    return LEFT + (x - X_MIN) / (X_MAX - X_MIN) * (RIGHT - LEFT)


def y_to_px(y: float) -> float:
    return BOTTOM - (y - Y_MIN) / (Y_MAX - Y_MIN) * (BOTTOM - TOP)


def points_to_pixels(points: np.ndarray) -> list[tuple[float, float]]:
    return [(x_to_px(float(x)), y_to_px(float(y))) for x, y in points[:, :2]]


def visible_points(points: np.ndarray) -> np.ndarray:
    return points[
        (points[:, 0] >= X_MIN)
        & (points[:, 0] <= X_MAX)
        & (points[:, 1] >= Y_MIN)
        & (points[:, 1] <= Y_MAX)
    ]


def draw_text_center(draw: ImageDraw.ImageDraw, xy: tuple[float, float], text: str, font, fill=BLACK) -> None:
    bbox = draw.textbbox((0, 0), text, font=font)
    draw.text((xy[0] - (bbox[2] - bbox[0]) / 2, xy[1] - (bbox[3] - bbox[1]) / 2), text, font=font, fill=fill)


def draw_rotated_center(base: Image.Image, center: tuple[float, float], text: str, font, angle: int = 90) -> None:
    bbox = ImageDraw.Draw(Image.new("RGBA", (1, 1))).textbbox((0, 0), text, font=font)
    layer = Image.new("RGBA", (bbox[2] - bbox[0] + 24, bbox[3] - bbox[1] + 24), (255, 255, 255, 0))
    ImageDraw.Draw(layer).text((12, 12), text, font=font, fill=BLACK)
    layer = layer.rotate(angle, expand=True, resample=Image.Resampling.BICUBIC)
    base.alpha_composite(layer, (int(center[0] - layer.width / 2), int(center[1] - layer.height / 2)))


def split_contiguous(indices: np.ndarray, gap: int = 4) -> list[np.ndarray]:
    if indices.size == 0:
        return []
    groups: list[np.ndarray] = []
    start = prev = int(indices[0])
    for value in indices[1:]:
        value = int(value)
        if value - prev > gap:
            groups.append(np.arange(start, prev + 1))
            start = value
        prev = value
    groups.append(np.arange(start, prev + 1))
    return groups


def bin_curve(points: np.ndarray, bin_width: float = 0.015) -> np.ndarray:
    if points.size == 0:
        return points
    bins = np.arange(points[:, 0].min(), points[:, 0].max() + bin_width, bin_width)
    binned = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (points[:, 0] >= lo) & (points[:, 0] < hi)
        if np.count_nonzero(mask) >= 2:
            binned.append(((lo + hi) / 2, float(np.median(points[mask, 1]))))
    return np.array(binned, dtype=float)


def digitize_first_figure_curves() -> tuple[np.ndarray, np.ndarray]:
    image = Image.open(FIRST_FIGURE).convert("RGB")
    arr = np.array(image)

    x_left, x_right = 172, 1276
    y_top, y_bottom = 60, 827
    y_data_max = 0.020
    crop = arr[y_top : y_bottom + 1, x_left : x_right + 1]

    blue_mask = (
        (crop[:, :, 2] > 130)
        & (crop[:, :, 1] > 70)
        & (crop[:, :, 1] < 170)
        & (crop[:, :, 0] < 60)
    )
    gray_mask = (
        (crop[:, :, 0] > 60)
        & (crop[:, :, 0] < 130)
        & (np.abs(crop[:, :, 0].astype(int) - crop[:, :, 1].astype(int)) < 12)
        & (np.abs(crop[:, :, 1].astype(int) - crop[:, :, 2].astype(int)) < 12)
    )

    def pixel_curve(mask: np.ndarray, gray_filter: bool = False) -> np.ndarray:
        out = []
        height, width = mask.shape
        for xx in range(width):
            ys = np.flatnonzero(mask[:, xx])
            if ys.size == 0:
                continue
            groups = split_contiguous(ys)
            group = max(groups, key=lambda values: float(values.mean()))
            x_pix = xx + x_left
            y_pix = float(group.mean() + y_top)
            p = (x_pix - x_left) / (x_right - x_left) * 5.0
            sigma = (y_bottom - y_pix) / (y_bottom - y_top) * y_data_max
            if 0.06 < p <= 5.02 and 0.0008 < sigma < 0.0205:
                if not gray_filter or p < 0.75 or y_pix > 220:
                    out.append((p, sigma))
        return np.array(out, dtype=float)

    blue_curve = bin_curve(pixel_curve(blue_mask), 0.012)
    gray_ys, gray_xs = np.where(gray_mask)
    gray_x = gray_xs + x_left
    gray_y = gray_ys + y_top
    gray_p = (gray_x - x_left) / (x_right - x_left) * 5.0
    gray_sigma = (y_bottom - gray_y) / (y_bottom - y_top) * y_data_max
    gray_keep = (
        (gray_p > 0.06)
        & (gray_p < 5.02)
        & (gray_sigma > 0.0022)
        & (gray_sigma < 0.0205)
        & (
            (gray_p < 0.75)
            | ((gray_y > 220) & (gray_sigma < 0.006))
        )
    )
    gray_raw = np.column_stack([gray_p[gray_keep], gray_sigma[gray_keep]])
    gray_curve = bin_curve(gray_raw, 0.02)
    return blue_curve, gray_curve


def digitize_red_points() -> np.ndarray:
    image = Image.open(SECOND_FIGURE).convert("RGB")
    arr = np.array(image)
    mask = (arr[:, :, 0] > 180) & (arr[:, :, 1] < 90) & (arr[:, :, 2] < 90)

    x_left, x_right = 204, 869
    y_top, y_bottom = 50, 955
    x_min, x_max = 0.2, 5.5
    crop = mask[y_top : y_bottom + 1, x_left : x_right + 1]
    visited = np.zeros_like(crop, dtype=bool)
    comps = []

    height, width = crop.shape
    for yy in range(height):
        for xx in np.flatnonzero(crop[yy] & ~visited[yy]):
            if visited[yy, xx] or not crop[yy, xx]:
                continue
            stack = [(yy, xx)]
            visited[yy, xx] = True
            pixels = []
            while stack:
                y, x = stack.pop()
                pixels.append((y, x))
                for ny in range(max(0, y - 1), min(height, y + 2)):
                    for nx in range(max(0, x - 1), min(width, x + 2)):
                        if (ny != y or nx != x) and crop[ny, nx] and not visited[ny, nx]:
                            visited[ny, nx] = True
                            stack.append((ny, nx))
            if len(pixels) < 8:
                continue
            ys = np.array([p[0] for p in pixels], dtype=float)
            xs = np.array([p[1] for p in pixels], dtype=float)
            cx = xs.mean() + x_left
            cy = ys.mean() + y_top
            p_val = x_min + (cx - x_left) / (x_right - x_left) * (x_max - x_min)
            sigma = (y_bottom - cy) / (y_bottom - y_top) * 0.005
            comps.append((p_val, sigma, len(pixels)))

    # Some red ROOT markers are split into left/right visible halves by overlaid
    # colored markers. Merge close components that share the same y level.
    components = sorted(comps)
    merged = []
    i = 0
    while i < len(components):
        group = [components[i]]
        j = i + 1
        while j < len(components) and len(group) < 2:
            close_x = components[j][0] - group[-1][0] < 0.22
            close_y = abs(components[j][1] - group[-1][1]) < 0.00005
            if not (close_x and close_y):
                break
            group.append(components[j])
            j += 1
        weight = np.array([item[2] for item in group], dtype=float)
        xs = np.array([item[0] for item in group], dtype=float)
        ys = np.array([item[1] for item in group], dtype=float)
        merged.append((float(np.average(xs, weights=weight)), float(np.average(ys, weights=weight))))
        i = j
    return np.array(merged, dtype=float)


def load_black_points() -> np.ndarray:
    rows = []
    with THIRD_CSV.open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                (
                    float(row["x_center"]),
                    float(row["iterative_sigma"]),
                    float(row["iterative_sigma_error"]),
                )
            )
    return np.array(rows, dtype=float)


def draw_polyline(draw: ImageDraw.ImageDraw, points: list[tuple[float, float]], fill, width: int) -> None:
    if len(points) >= 2:
        draw.line([(round(x), round(y)) for x, y in points], fill=fill, width=width, joint="curve")


def draw_dashed_polyline(draw: ImageDraw.ImageDraw, points: list[tuple[float, float]], fill, width: int, dash: int = 34, gap: int = 22) -> None:
    if len(points) < 2:
        return
    draw_dash = True
    remain = dash
    for start, end in zip(points[:-1], points[1:]):
        x0, y0 = start
        x1, y1 = end
        length = math.hypot(x1 - x0, y1 - y0)
        if length == 0:
            continue
        consumed = 0.0
        while consumed < length:
            step = min(remain, length - consumed)
            t0 = consumed / length
            t1 = (consumed + step) / length
            a = (x0 + (x1 - x0) * t0, y0 + (y1 - y0) * t0)
            b = (x0 + (x1 - x0) * t1, y0 + (y1 - y0) * t1)
            if draw_dash:
                draw.line((a, b), fill=fill, width=width)
            consumed += step
            remain -= step
            if remain <= 1e-9:
                draw_dash = not draw_dash
                remain = dash if draw_dash else gap


def y_tick_label(value: float) -> str:
    if value == 0:
        return "0"
    return f"{value:.3f}".rstrip("0").rstrip(".")


def draw_axes(draw: ImageDraw.ImageDraw, tick_font, label_font) -> None:
    axis_w = 8
    draw.rectangle((LEFT, TOP, RIGHT, BOTTOM), outline=BLACK, width=axis_w)
    major_len, minor_len = 38, 20

    for tick in np.arange(0.5, 5.51, 0.5):
        x = x_to_px(float(tick))
        draw.line((x, BOTTOM, x, BOTTOM - minor_len), fill=BLACK, width=4)
        draw.line((x, TOP, x, TOP + minor_len), fill=BLACK, width=4)

    for tick in range(0, 6):
        x = x_to_px(float(tick))
        draw.line((x, BOTTOM, x, BOTTOM - major_len), fill=BLACK, width=5)
        draw.line((x, TOP, x, TOP + major_len), fill=BLACK, width=5)
        draw_text_center(draw, (x, BOTTOM + 76), str(tick), tick_font)

    y_major_ticks = [0.0, 0.005, 0.010]
    for tick in np.arange(0.0025, Y_MAX, 0.0025):
        if any(abs(float(tick) - major) < 1e-12 for major in y_major_ticks):
            continue
        y = y_to_px(float(tick))
        draw.line((LEFT, y, LEFT + minor_len, y), fill=BLACK, width=4)
        draw.line((RIGHT, y, RIGHT - minor_len, y), fill=BLACK, width=4)

    for tick in y_major_ticks:
        label = y_tick_label(tick)
        y = y_to_px(tick)
        draw.line((LEFT, y, LEFT + major_len, y), fill=BLACK, width=5)
        draw.line((RIGHT, y, RIGHT - major_len, y), fill=BLACK, width=5)
        bbox = draw.textbbox((0, 0), label, font=tick_font)
        draw.text((LEFT - 24 - (bbox[2] - bbox[0]), y - (bbox[3] - bbox[1]) / 2), label, font=tick_font, fill=BLACK)

    draw_text_center(draw, ((LEFT + RIGHT) / 2, HEIGHT - 72), "Track momentum p (GeV/c)", label_font)


def draw_title_box(draw: ImageDraw.ImageDraw, title_font) -> None:
    title = "sPHENIX pointing-resolution comparison"
    title_bbox = draw.textbbox((0, 0), title, font=title_font)
    title_w = title_bbox[2] - title_bbox[0]
    title_h = title_bbox[3] - title_bbox[1]
    title_x = (LEFT + RIGHT) / 2
    title_y = TOP + 74
    pad_x, pad_y = 46, 24
    draw.rectangle(
        (
            title_x - title_w / 2 - pad_x,
            title_y - title_h / 2 - pad_y,
            title_x + title_w / 2 + pad_x,
            title_y + title_h / 2 + pad_y,
        ),
        fill="white",
        outline=BLACK,
        width=8,
    )
    draw_text_center(draw, (title_x, title_y), title, title_font)


def draw_png_pdf(blue_curve: np.ndarray, gray_curve: np.ndarray, red_points: np.ndarray, black_points: np.ndarray) -> None:
    image = Image.new("RGBA", (WIDTH, HEIGHT), "white")
    draw = ImageDraw.Draw(image)

    tick_font = load_font(58)
    label_font = load_font(70)
    legend_font = load_font(48)
    title_font = load_font(59)

    draw_polyline(draw, points_to_pixels(visible_points(blue_curve)), BLUE, 8)
    draw_dashed_polyline(draw, points_to_pixels(visible_points(gray_curve)), GRAY, 7)

    for x, y in red_points:
        xpix, ypix = x_to_px(float(x)), y_to_px(float(y))
        r = 11
        draw.rectangle((xpix - r, ypix - r, xpix + r, ypix + r), fill=RED, outline=RED)

    for x, y, err in black_points:
        xpix, ypix = x_to_px(float(x)), y_to_px(float(y))
        ylo, yhi = y_to_px(float(y - err)), y_to_px(float(y + err))
        draw.line((xpix, yhi, xpix, ylo), fill=BLACK, width=4)
        draw.line((xpix - 16, yhi, xpix + 16, yhi), fill=BLACK, width=4)
        draw.line((xpix - 16, ylo, xpix + 16, ylo), fill=BLACK, width=4)
        r = 8
        draw.ellipse((xpix - r, ypix - r, xpix + r, ypix + r), fill=BLACK, outline=BLACK)

    draw_axes(draw, tick_font, label_font)
    draw_rotated_center(image, (64, (TOP + BOTTOM) / 2), "\u03c3(dcaxy) (cm)", label_font)
    draw_title_box(draw, title_font)

    # Legend in the same direct, unboxed style as the STAR overlay.
    lx, ly = 950, 355
    draw.line((lx, ly, lx + 115, ly), fill=BLUE, width=8)
    draw.text((lx + 145, ly - 31), "PXL + INTT fit", font=legend_font, fill=BLACK)
    ly += 76
    draw_dashed_polyline(draw, [(lx, ly), (lx + 115, ly)], GRAY, 7)
    draw.text((lx + 145, ly - 31), "3 PXL layers only", font=legend_font, fill=BLACK)
    ly += 76
    draw.rectangle((lx + 47, ly - 11, lx + 69, ly + 11), fill=RED, outline=RED)
    draw.text((lx + 145, ly - 31), "sPHENIX simulation", font=legend_font, fill=BLACK)
    ly += 76
    draw.ellipse((lx + 49, ly - 9, lx + 67, ly + 9), fill=BLACK, outline=BLACK)
    draw.text((lx + 145, ly - 31), "data: iterative Gaussian", font=legend_font, fill=BLACK)

    png_path = OUT_DIR / "sphenix_merged_resolution_overlay.png"
    pdf_path = OUT_DIR / "sphenix_merged_resolution_overlay.pdf"
    image.convert("RGB").save(png_path)
    image.convert("RGB").save(pdf_path, "PDF", resolution=300.0)


def svg_text(x: float, y: float, text: str, size: int, anchor: str = "start", extra: str = "") -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" text-anchor="{anchor}" fill="black" {extra}>{escape(text)}</text>'
    )


def draw_svg(blue_curve: np.ndarray, gray_curve: np.ndarray, red_points: np.ndarray, black_points: np.ndarray) -> None:
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="white"/>',
    ]
    blue = " ".join(f"{x_to_px(x):.1f},{y_to_px(y):.1f}" for x, y in visible_points(blue_curve))
    gray = " ".join(f"{x_to_px(x):.1f},{y_to_px(y):.1f}" for x, y in visible_points(gray_curve))
    parts.append(f'<polyline points="{blue}" fill="none" stroke="rgb{BLUE}" stroke-width="8" stroke-linejoin="round"/>')
    parts.append(f'<polyline points="{gray}" fill="none" stroke="rgb{GRAY}" stroke-width="7" stroke-dasharray="34 22" stroke-linejoin="round"/>')
    for x, y in red_points:
        parts.append(f'<rect x="{x_to_px(x) - 11:.1f}" y="{y_to_px(y) - 11:.1f}" width="22" height="22" fill="red"/>')
    for x, y, err in black_points:
        xpix, ypix = x_to_px(float(x)), y_to_px(float(y))
        yhi, ylo = y_to_px(float(y + err)), y_to_px(float(y - err))
        parts.append(f'<line x1="{xpix:.1f}" y1="{yhi:.1f}" x2="{xpix:.1f}" y2="{ylo:.1f}" stroke="black" stroke-width="4"/>')
        parts.append(f'<line x1="{xpix - 16:.1f}" y1="{yhi:.1f}" x2="{xpix + 16:.1f}" y2="{yhi:.1f}" stroke="black" stroke-width="4"/>')
        parts.append(f'<line x1="{xpix - 16:.1f}" y1="{ylo:.1f}" x2="{xpix + 16:.1f}" y2="{ylo:.1f}" stroke="black" stroke-width="4"/>')
        parts.append(f'<circle cx="{xpix:.1f}" cy="{ypix:.1f}" r="8" fill="black"/>')

    parts.append(f'<rect x="{LEFT}" y="{TOP}" width="{RIGHT - LEFT}" height="{BOTTOM - TOP}" fill="none" stroke="black" stroke-width="8"/>')
    for tick in np.arange(0.5, 5.51, 0.5):
        x = x_to_px(float(tick))
        parts.append(f'<line x1="{x:.1f}" y1="{BOTTOM}" x2="{x:.1f}" y2="{BOTTOM - 20}" stroke="black" stroke-width="4"/>')
        parts.append(f'<line x1="{x:.1f}" y1="{TOP}" x2="{x:.1f}" y2="{TOP + 20}" stroke="black" stroke-width="4"/>')
    for tick in range(0, 6):
        x = x_to_px(float(tick))
        parts.append(f'<line x1="{x:.1f}" y1="{BOTTOM}" x2="{x:.1f}" y2="{BOTTOM - 38}" stroke="black" stroke-width="5"/>')
        parts.append(f'<line x1="{x:.1f}" y1="{TOP}" x2="{x:.1f}" y2="{TOP + 38}" stroke="black" stroke-width="5"/>')
        parts.append(svg_text(x, BOTTOM + 96, str(tick), 58, "middle"))
    y_major_ticks = [0.0, 0.005, 0.010]
    for tick in np.arange(0.0025, Y_MAX, 0.0025):
        if any(abs(float(tick) - major) < 1e-12 for major in y_major_ticks):
            continue
        y = y_to_px(float(tick))
        parts.append(f'<line x1="{LEFT}" y1="{y:.1f}" x2="{LEFT + 20}" y2="{y:.1f}" stroke="black" stroke-width="4"/>')
        parts.append(f'<line x1="{RIGHT}" y1="{y:.1f}" x2="{RIGHT - 20}" y2="{y:.1f}" stroke="black" stroke-width="4"/>')
    for tick in y_major_ticks:
        label = y_tick_label(tick)
        y = y_to_px(tick)
        parts.append(f'<line x1="{LEFT}" y1="{y:.1f}" x2="{LEFT + 38}" y2="{y:.1f}" stroke="black" stroke-width="5"/>')
        parts.append(f'<line x1="{RIGHT}" y1="{y:.1f}" x2="{RIGHT - 38}" y2="{y:.1f}" stroke="black" stroke-width="5"/>')
        parts.append(svg_text(LEFT - 24, y + 20, label, 58, "end"))
    parts.append(svg_text((LEFT + RIGHT) / 2, HEIGHT - 52, "Track momentum p (GeV/c)", 70, "middle"))
    parts.append(svg_text(0, 0, "\u03c3(dcaxy) (cm)", 70, "middle", f'transform="translate(64 {(TOP + BOTTOM) / 2:.1f}) rotate(-90)"'))

    title = "sPHENIX pointing-resolution comparison"
    title_x = (LEFT + RIGHT) / 2
    title_y = TOP + 74
    rect_w, rect_h = 1180, 108
    parts.append(f'<rect x="{title_x - rect_w / 2:.1f}" y="{title_y - rect_h / 2:.1f}" width="{rect_w}" height="{rect_h}" fill="white" stroke="black" stroke-width="8"/>')
    parts.append(svg_text(title_x, title_y + 20, title, 59, "middle"))

    lx, ly = 950, 355
    parts.append(f'<line x1="{lx}" y1="{ly}" x2="{lx + 115}" y2="{ly}" stroke="rgb{BLUE}" stroke-width="8"/>')
    parts.append(svg_text(lx + 145, ly + 18, "PXL + INTT fit", 48))
    ly += 76
    parts.append(f'<line x1="{lx}" y1="{ly}" x2="{lx + 115}" y2="{ly}" stroke="rgb{GRAY}" stroke-width="7" stroke-dasharray="34 22"/>')
    parts.append(svg_text(lx + 145, ly + 18, "3 PXL layers only", 48))
    ly += 76
    parts.append(f'<rect x="{lx + 47}" y="{ly - 11}" width="22" height="22" fill="red"/>')
    parts.append(svg_text(lx + 145, ly + 18, "sPHENIX simulation", 48))
    ly += 76
    parts.append(f'<circle cx="{lx + 58}" cy="{ly}" r="9" fill="black"/>')
    parts.append(svg_text(lx + 145, ly + 18, "data: iterative Gaussian", 48))

    parts.append("</svg>\n")
    (OUT_DIR / "sphenix_merged_resolution_overlay.svg").write_text("\n".join(parts), encoding="utf-8")


def save_csv(blue_curve: np.ndarray, gray_curve: np.ndarray, red_points: np.ndarray, black_points: np.ndarray) -> None:
    out = OUT_DIR / "sphenix_merged_resolution_data.csv"
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["series", "p_GeV_c", "sigma_dcaxy_cm", "sigma_error_cm"])
        for x, y in blue_curve:
            writer.writerow(["PXL+INTT fit", f"{x:.5f}", f"{y:.8f}", ""])
        for x, y in gray_curve:
            writer.writerow(["3 PXL layers only", f"{x:.5f}", f"{y:.8f}", ""])
        for x, y in red_points:
            writer.writerow(["sPHENIX simulation red", f"{x:.5f}", f"{y:.8f}", ""])
        for x, y, err in black_points:
            writer.writerow(["data iterative Gaussian black", f"{x:.5f}", f"{y:.8f}", f"{err:.8f}"])


def main() -> None:
    blue_curve, gray_curve = digitize_first_figure_curves()
    red_points = digitize_red_points()
    black_points = load_black_points()
    save_csv(blue_curve, gray_curve, red_points, black_points)
    draw_png_pdf(blue_curve, gray_curve, red_points, black_points)
    draw_svg(blue_curve, gray_curve, red_points, black_points)


if __name__ == "__main__":
    main()
