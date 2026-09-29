"""Deterministic figure-writing helpers."""

from __future__ import annotations

from pathlib import Path


def save_figure(fig, path, dpi: int = 300, bbox_inches=None) -> None:
    """Write a figure without timestamp metadata in PDF output."""
    path = Path(path)
    kwargs = {"dpi": dpi}
    if bbox_inches is not None:
        kwargs["bbox_inches"] = bbox_inches
    if path.suffix.lower() == ".pdf":
        kwargs["metadata"] = {"CreationDate": None, "ModDate": None}
    fig.savefig(path, **kwargs)


def save_figure_pair(fig, output_dir: Path, stem: str) -> None:
    """Write one figure as a PNG/PDF pair and release its pyplot state."""
    import matplotlib.pyplot as plt

    for suffix in ("png", "pdf"):
        save_figure(
            fig,
            Path(output_dir) / f"{stem}.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)
