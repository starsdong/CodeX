#!/usr/bin/env python3
"""Overlay digitized STAR pion DCAxy points with the PXL estimate curve.

This version draws the plot directly with Pillow so it does not depend on
Matplotlib's font cache.
"""

from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
from PIL import Image, ImageDraw, ImageFont


OUT_DIR = Path(__file__).resolve().parent
FONT_REGULAR = "/System/Library/Fonts/Supplemental/Arial.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"

# Digitized from /Users/starsdong/Desktop/STAR_data.png, left DCA_xy panel,
# keeping only the filled black pion markers.
STAR_PION_DCAXY = np.array(
    [
        [0.213, 160.2],
        [0.239, 147.6],
        [0.264, 131.1],
        [0.289, 120.5],
        [0.315, 108.7],
        [0.340, 100.3],
        [0.365, 91.6],
        [0.390, 86.1],
        [0.415, 79.6],
        [0.440, 75.3],
        [0.465, 70.0],
        [0.484, 68.0],
        [0.515, 63.1],
        [0.541, 60.2],
        [0.576, 57.0],
        [0.627, 52.5],
        [0.676, 48.4],
        [0.726, 45.2],
        [0.766, 43.1],
        [0.805, 41.7],
        [0.849, 39.3],
        [0.892, 38.0],
        [0.941, 36.6],
        [1.002, 33.8],
        [1.091, 32.3],
        [1.161, 30.7],
        [1.214, 28.9],
        [1.287, 28.1],
        [1.358, 26.9],
        [1.439, 26.2],
        [1.518, 24.8],
        [1.616, 23.8],
        [1.713, 22.6],
        [1.959, 20.7],
        [2.030, 20.3],
        [2.171, 19.3],
        [2.301, 18.5],
        [2.505, 17.9],
        [2.703, 17.3],
        [2.942, 16.7],
        [3.160, 16.5],
        [3.487, 16.1],
        [3.679, 15.9],
        [3.882, 15.5],
        [4.077, 15.3],
        [4.321, 15.2],
        [4.642, 14.8],
    ]
)

X_MIN, X_MAX = 0.18, 5.0
Y_MIN, Y_MAX = 0.0, 170.0
WIDTH, HEIGHT = 1800, 1800
LEFT, TOP, RIGHT, BOTTOM = 190, 75, 1665, 1550


def pxl_estimate(p):
    """STAR PXL estimate from the reference plot: 9.74 oplus 24.25 / p."""

    return np.sqrt(9.74**2 + (24.25 / p) ** 2)


def load_font(path, size):
    try:
        return ImageFont.truetype(path, size)
    except OSError:
        return ImageFont.load_default()


def x_to_px(p):
    frac = (np.log10(p) - np.log10(X_MIN)) / (np.log10(X_MAX) - np.log10(X_MIN))
    return LEFT + frac * (RIGHT - LEFT)


def y_to_px(y):
    frac = (y - Y_MIN) / (Y_MAX - Y_MIN)
    return BOTTOM - frac * (BOTTOM - TOP)


def draw_text_center(draw, xy, text, font, fill="black"):
    bbox = draw.textbbox((0, 0), text, font=font)
    width = bbox[2] - bbox[0]
    height = bbox[3] - bbox[1]
    draw.text((xy[0] - width / 2, xy[1] - height / 2), text, font=font, fill=fill)


def draw_rotated_center(base, center, text, font, angle=90):
    text_bbox = ImageDraw.Draw(Image.new("RGBA", (1, 1))).textbbox((0, 0), text, font=font)
    text_w = text_bbox[2] - text_bbox[0]
    text_h = text_bbox[3] - text_bbox[1]
    layer = Image.new("RGBA", (text_w + 20, text_h + 20), (255, 255, 255, 0))
    layer_draw = ImageDraw.Draw(layer)
    layer_draw.text((10, 10), text, font=font, fill="black")
    layer = layer.rotate(angle, expand=True, resample=Image.Resampling.BICUBIC)
    base.alpha_composite(layer, (int(center[0] - layer.width / 2), int(center[1] - layer.height / 2)))


def draw_polyline(draw, points, fill, width):
    int_points = [(int(round(x)), int(round(y))) for x, y in points]
    if len(int_points) > 1:
        draw.line(int_points, fill=fill, width=width, joint="curve")


def draw_png_pdf():
    image = Image.new("RGBA", (WIDTH, HEIGHT), "white")
    draw = ImageDraw.Draw(image)

    tick_font = load_font(FONT_REGULAR, 58)
    label_font = load_font(FONT_REGULAR, 70)
    text_font = load_font(FONT_REGULAR, 67)
    legend_font = load_font(FONT_REGULAR, 52)
    title_font = load_font(FONT_REGULAR, 66)

    # Estimate curve.
    p_curve = np.logspace(np.log10(X_MIN), np.log10(X_MAX), 600)
    curve_points = [(x_to_px(p), y_to_px(y)) for p, y in zip(p_curve, pxl_estimate(p_curve))]
    draw_polyline(draw, curve_points, "black", 8)

    # Pion data markers.
    for p, y in STAR_PION_DCAXY:
        x_pix, y_pix = x_to_px(p), y_to_px(y)
        r = 11
        draw.ellipse((x_pix - r, y_pix - r, x_pix + r, y_pix + r), fill="black", outline="black")

    # Frame and ticks.
    axis_w = 8
    draw.rectangle((LEFT, TOP, RIGHT, BOTTOM), outline="black", width=axis_w)

    major_tick = 36
    minor_tick = 19
    x_major = [0.2, 0.5, 1.0, 2.0, 3.0]
    for decade in [0.1, 1.0]:
        for n in range(2, 10):
            tick = n * decade
            if X_MIN <= tick <= X_MAX and all(abs(tick - m) > 1e-6 for m in x_major):
                x = x_to_px(tick)
                draw.line((x, BOTTOM, x, BOTTOM - minor_tick), fill="black", width=4)
                draw.line((x, TOP, x, TOP + minor_tick), fill="black", width=4)

    for tick, label in zip(x_major, ["0.2", "0.5", "1", "2", "3"]):
        x = x_to_px(tick)
        draw.line((x, BOTTOM, x, BOTTOM - major_tick), fill="black", width=5)
        draw.line((x, TOP, x, TOP + major_tick), fill="black", width=5)
        draw_text_center(draw, (x, BOTTOM + 68), label, tick_font)

    y_major = [0, 50, 100, 150]
    for tick in range(10, 170, 10):
        if tick not in y_major:
            y = y_to_px(tick)
            draw.line((LEFT, y, LEFT + minor_tick, y), fill="black", width=4)
            draw.line((RIGHT, y, RIGHT - minor_tick, y), fill="black", width=4)

    for tick in y_major:
        y = y_to_px(tick)
        draw.line((LEFT, y, LEFT + major_tick, y), fill="black", width=5)
        draw.line((RIGHT, y, RIGHT - major_tick, y), fill="black", width=5)
        bbox = draw.textbbox((0, 0), str(tick), font=tick_font)
        draw.text((LEFT - 25 - (bbox[2] - bbox[0]), y - (bbox[3] - bbox[1]) / 2), str(tick), font=tick_font, fill="black")

    draw_text_center(draw, ((LEFT + RIGHT) / 2, HEIGHT - 70), "Momentum p (GeV/c)", label_font)
    draw_rotated_center(image, (48, (TOP + BOTTOM) / 2), "\u03c3xy (\u03bcm)", label_font, angle=90)

    # STAR-style title box.
    title = "Au+Au \u221asNN = 200 GeV, 0-80%"
    title_bbox = draw.textbbox((0, 0), title, font=title_font)
    title_w = title_bbox[2] - title_bbox[0]
    title_h = title_bbox[3] - title_bbox[1]
    title_x = (LEFT + RIGHT) / 2
    title_y = TOP + 70
    pad_x, pad_y = 52, 24
    title_rect = (
        title_x - title_w / 2 - pad_x,
        title_y - title_h / 2 - pad_y,
        title_x + title_w / 2 + pad_x,
        title_y + title_h / 2 + pad_y,
    )
    draw.rectangle(title_rect, fill="white", outline="black", width=axis_w)
    draw_text_center(draw, (title_x, title_y), title, title_font)

    # Legend in the original figure's direct, unboxed style.
    leg_x, leg_y = 930, 270
    draw.line((leg_x, leg_y, leg_x + 115, leg_y), fill="black", width=8)
    draw.text((leg_x + 150, leg_y - 33), "STAR PXL estimate", font=legend_font, fill="black")
    marker_y = leg_y + 80
    draw.ellipse((leg_x + 45, marker_y - 11, leg_x + 67, marker_y + 11), fill="black", outline="black")
    draw.text((leg_x + 150, marker_y - 33), "STAR \u03c0\u00b1", font=legend_font, fill="black")

    draw.text((LEFT + 72, BOTTOM - 118), "DCAxy", font=text_font, fill="black")

    png_path = OUT_DIR / "star_pion_pxl_overlay.png"
    pdf_path = OUT_DIR / "star_pion_pxl_overlay.pdf"
    image.convert("RGB").save(png_path)
    image.convert("RGB").save(pdf_path, "PDF", resolution=300.0)


def svg_text(x, y, text, size, anchor="start", extra=""):
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" text-anchor="{anchor}" fill="black" {extra}>{escape(text)}</text>'
    )


def draw_svg():
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="white"/>',
    ]
    p_curve = np.logspace(np.log10(X_MIN), np.log10(X_MAX), 600)
    curve = " ".join(f"{x_to_px(p):.1f},{y_to_px(y):.1f}" for p, y in zip(p_curve, pxl_estimate(p_curve)))
    parts.append(f'<polyline points="{curve}" fill="none" stroke="black" stroke-width="8" stroke-linejoin="round"/>')
    for p, y in STAR_PION_DCAXY:
        parts.append(f'<circle cx="{x_to_px(p):.1f}" cy="{y_to_px(y):.1f}" r="11" fill="black"/>')

    parts.append(f'<rect x="{LEFT}" y="{TOP}" width="{RIGHT - LEFT}" height="{BOTTOM - TOP}" fill="none" stroke="black" stroke-width="8"/>')
    x_major = [0.2, 0.5, 1.0, 2.0, 3.0]
    for decade in [0.1, 1.0]:
        for n in range(2, 10):
            tick = n * decade
            if X_MIN <= tick <= X_MAX and all(abs(tick - m) > 1e-6 for m in x_major):
                x = x_to_px(tick)
                parts.append(f'<line x1="{x:.1f}" y1="{BOTTOM}" x2="{x:.1f}" y2="{BOTTOM - 19}" stroke="black" stroke-width="4"/>')
                parts.append(f'<line x1="{x:.1f}" y1="{TOP}" x2="{x:.1f}" y2="{TOP + 19}" stroke="black" stroke-width="4"/>')
    for tick, label in zip(x_major, ["0.2", "0.5", "1", "2", "3"]):
        x = x_to_px(tick)
        parts.append(f'<line x1="{x:.1f}" y1="{BOTTOM}" x2="{x:.1f}" y2="{BOTTOM - 36}" stroke="black" stroke-width="5"/>')
        parts.append(f'<line x1="{x:.1f}" y1="{TOP}" x2="{x:.1f}" y2="{TOP + 36}" stroke="black" stroke-width="5"/>')
        parts.append(svg_text(x, BOTTOM + 90, label, 58, anchor="middle"))

    y_major = [0, 50, 100, 150]
    for tick in range(10, 170, 10):
        if tick not in y_major:
            y = y_to_px(tick)
            parts.append(f'<line x1="{LEFT}" y1="{y:.1f}" x2="{LEFT + 19}" y2="{y:.1f}" stroke="black" stroke-width="4"/>')
            parts.append(f'<line x1="{RIGHT}" y1="{y:.1f}" x2="{RIGHT - 19}" y2="{y:.1f}" stroke="black" stroke-width="4"/>')
    for tick in y_major:
        y = y_to_px(tick)
        parts.append(f'<line x1="{LEFT}" y1="{y:.1f}" x2="{LEFT + 36}" y2="{y:.1f}" stroke="black" stroke-width="5"/>')
        parts.append(f'<line x1="{RIGHT}" y1="{y:.1f}" x2="{RIGHT - 36}" y2="{y:.1f}" stroke="black" stroke-width="5"/>')
        parts.append(svg_text(LEFT - 28, y + 20, str(tick), 58, anchor="end"))

    parts.append(svg_text((LEFT + RIGHT) / 2, HEIGHT - 50, "Momentum p (GeV/c)", 70, anchor="middle"))
    parts.append(svg_text(0, 0, "\u03c3xy (\u03bcm)", 70, anchor="middle", extra=f'transform="translate(48 {(TOP + BOTTOM) / 2:.1f}) rotate(-90)"'))

    title = "Au+Au \u221asNN = 200 GeV, 0-80%"
    title_x, title_y = (LEFT + RIGHT) / 2, TOP + 70
    rect_w, rect_h = 935, 105
    parts.append(f'<rect x="{title_x - rect_w / 2:.1f}" y="{title_y - rect_h / 2:.1f}" width="{rect_w}" height="{rect_h}" fill="white" stroke="black" stroke-width="8"/>')
    parts.append(svg_text(title_x, title_y + 22, title, 66, anchor="middle"))

    leg_x, leg_y = 930, 270
    parts.append(f'<line x1="{leg_x}" y1="{leg_y}" x2="{leg_x + 115}" y2="{leg_y}" stroke="black" stroke-width="8"/>')
    parts.append(svg_text(leg_x + 150, leg_y + 18, "STAR PXL estimate", 52))
    marker_y = leg_y + 80
    parts.append(f'<circle cx="{leg_x + 56}" cy="{marker_y}" r="11" fill="black"/>')
    parts.append(svg_text(leg_x + 150, marker_y + 18, "STAR \u03c0\u00b1", 52))
    parts.append(svg_text(LEFT + 72, BOTTOM - 70, "DCAxy", 67))
    parts.append("</svg>\n")

    (OUT_DIR / "star_pion_pxl_overlay.svg").write_text("\n".join(parts), encoding="utf-8")


def main():
    np.savetxt(
        OUT_DIR / "star_pion_dcaxy_digitized.csv",
        STAR_PION_DCAXY,
        delimiter=",",
        header="p_GeV_c,sigma_xy_um",
        comments="",
        fmt=["%.4f", "%.2f"],
    )
    draw_png_pdf()
    draw_svg()


if __name__ == "__main__":
    main()
