#!/usr/bin/env python3
"""Build a labelled contact sheet for the nine submission figures."""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, ImageOps


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    figure_dir = Path(args.figure_dir)
    paths = sorted(figure_dir.glob("Figure*.png"))
    if len(paths) != 9:
        raise SystemExit(f"expected 9 PNG figures, found {len(paths)}")

    columns = 3
    tile_width = 520
    tile_height = 430
    margin = 20
    label_height = 34
    rows = (len(paths) + columns - 1) // columns
    sheet = Image.new(
        "RGB",
        (
            columns * tile_width + (columns + 1) * margin,
            rows * (tile_height + label_height)
            + (rows + 1) * margin,
        ),
        "white",
    )
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.load_default()

    for index, path in enumerate(paths):
        row, column = divmod(index, columns)
        x = margin + column * (tile_width + margin)
        y = margin + row * (tile_height + label_height + margin)
        image = Image.open(path).convert("RGB")
        image.thumbnail((tile_width, tile_height), Image.Resampling.LANCZOS)
        framed = ImageOps.expand(image, border=1, fill="#D1D5DB")
        offset_x = x + (tile_width - framed.width) // 2
        offset_y = y + (tile_height - framed.height) // 2
        sheet.paste(framed, (offset_x, offset_y))
        draw.text(
            (x, y + tile_height + 8),
            path.stem,
            fill="#111827",
            font=font,
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(output, dpi=(200, 200))
    print(output)


if __name__ == "__main__":
    main()
