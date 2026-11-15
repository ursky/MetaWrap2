"""Every PNG under docs/figures carries a resolution tag of at least 300 dpi.

The tag is read from the PNG pHYs chunk with struct, so the check needs nothing beyond the
standard library. A figure with no pHYs chunk fails too: an untagged PNG is shown at whatever
resolution the viewer assumes, which is how 150 dpi figures went unnoticed.
"""

import glob
import os
import struct

import pytest

FIGURES = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs", "figures"
)
MIN_DPI = 300
#: matplotlib and Pillow store 300 dpi as 11811 pixels per metre, which reads back as 299.9994 dpi.
TOLERANCE = 0.5
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def png_dpi(path):
    """(x_dpi, y_dpi) from the pHYs chunk, or None when the file has no metre-based pHYs chunk."""
    with open(path, "rb") as handle:
        if handle.read(8) != PNG_SIGNATURE:
            raise ValueError("%s is not a PNG" % path)
        while True:
            header = handle.read(8)
            if len(header) < 8:
                return None
            length, kind = struct.unpack(">I4s", header)
            if kind == b"pHYs":
                x_ppu, y_ppu, unit = struct.unpack(">IIB", handle.read(9))
                return (x_ppu * 0.0254, y_ppu * 0.0254) if unit == 1 else None
            if kind in (b"IDAT", b"IEND"):  # pHYs must come before the image data
                return None
            handle.seek(length + 4, os.SEEK_CUR)  # chunk data and CRC


FIGURE_FILES = sorted(glob.glob(os.path.join(FIGURES, "*.png")))


def test_docs_figures_exist():
    assert FIGURE_FILES, "no PNGs found under docs/figures"


@pytest.mark.parametrize("path", FIGURE_FILES, ids=os.path.basename)
def test_docs_figure_is_at_least_300_dpi(path):
    dpi = png_dpi(path)
    assert dpi is not None, "%s has no pHYs resolution tag" % os.path.basename(path)
    assert min(dpi) >= MIN_DPI - TOLERANCE, "%s is tagged %.1f x %.1f dpi, below %d" % (
        os.path.basename(path),
        dpi[0],
        dpi[1],
        MIN_DPI,
    )
