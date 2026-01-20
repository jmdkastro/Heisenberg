"""
DS9 region file I/O for peak visualization.

Thin wrapper around the astropy-regions package for DS9 region
file export. Used for displaying peaks in DS9.

Translates IDL ds9_display_peaks.pro and mask_ds9_box_vertices.pro functionality.
"""

from pathlib import Path
from typing import List, Optional, Union, Tuple, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from heisenberg.peaks.detection import DetectedPeak


def write_peak_regions(
    peaks: List['DetectedPeak'],
    output_path: Union[str, Path],
    color: str = 'green',
    symbol: str = 'cross',
    width: int = 1,
    size: int = 10,
    coordinate_system: str = 'image',
    include_labels: bool = False,
) -> None:
    """
    Write peak positions to DS9 region file.

    Uses the astropy-regions package for proper DS9 format output.

    Args:
        peaks: List of DetectedPeak objects with x, y positions
        output_path: Path for output .reg file
        color: Region color (green, red, blue, cyan, magenta, yellow, white)
        symbol: Point symbol (cross, x, circle, box, diamond)
        width: Line width
        size: Point size in pixels
        coordinate_system: Coordinate system (image, physical, fk5, galactic)
        include_labels: If True, add peak ID labels

    Example:
        >>> from heisenberg.peaks import find_peaks
        >>> peaks = find_peaks(image, config)
        >>> write_peak_regions(peaks, "star_peaks.reg", color="green")
    """
    try:
        from regions import PointPixelRegion, PixCoord, Regions
        from regions import write_ds9
    except ImportError:
        # Fallback to manual DS9 format writing
        _write_peak_regions_manual(
            peaks, output_path, color, symbol, width, size,
            coordinate_system, include_labels
        )
        return

    output_path = Path(output_path)

    # Create region list
    region_list = []
    for i, peak in enumerate(peaks):
        center = PixCoord(x=peak.x, y=peak.y)
        region = PointPixelRegion(center=center)

        # Set visual properties
        region.visual = {
            'color': color,
            'width': width,
            'point': symbol,
        }

        if include_labels:
            region.meta = {'text': str(i + 1)}

        region_list.append(region)

    # Write using regions package
    regions = Regions(region_list)
    regions.write(str(output_path), format='ds9', overwrite=True)


def _write_peak_regions_manual(
    peaks: List['DetectedPeak'],
    output_path: Union[str, Path],
    color: str = 'green',
    symbol: str = 'cross',
    width: int = 1,
    size: int = 10,
    coordinate_system: str = 'image',
    include_labels: bool = False,
) -> None:
    """Fallback manual DS9 region file writer."""
    output_path = Path(output_path)

    # Map symbol names
    symbol_map = {
        'cross': 'cross',
        'x': 'x',
        'circle': 'circle',
        'box': 'box',
        'diamond': 'diamond',
    }
    ds9_symbol = symbol_map.get(symbol, 'cross')

    lines = [
        '# Region file format: DS9 version 4.1',
        f'global color={color} width={width} point={ds9_symbol} {size}',
        f'{coordinate_system}',
    ]

    for i, peak in enumerate(peaks):
        # DS9 uses 1-indexed coordinates for image coordinate system
        x = peak.x + 1 if coordinate_system == 'image' else peak.x
        y = peak.y + 1 if coordinate_system == 'image' else peak.y

        if include_labels:
            lines.append(f'point({x:.2f},{y:.2f}) # text={{{i+1}}}')
        else:
            lines.append(f'point({x:.2f},{y:.2f})')

    output_path.write_text('\n'.join(lines) + '\n')


def write_circle_regions(
    centers: List[tuple],
    radii: List[float],
    output_path: Union[str, Path],
    color: str = 'cyan',
    width: int = 1,
    coordinate_system: str = 'image',
    labels: Optional[List[str]] = None,
) -> None:
    """
    Write circular regions to DS9 file.

    Args:
        centers: List of (x, y) center coordinates
        radii: List of radii (same length as centers)
        output_path: Path for output .reg file
        color: Region color
        width: Line width
        coordinate_system: Coordinate system
        labels: Optional list of labels
    """
    try:
        from regions import CirclePixelRegion, PixCoord, Regions
    except ImportError:
        _write_circle_regions_manual(
            centers, radii, output_path, color, width,
            coordinate_system, labels
        )
        return

    output_path = Path(output_path)
    region_list = []

    for i, ((x, y), radius) in enumerate(zip(centers, radii)):
        center = PixCoord(x=x, y=y)
        region = CirclePixelRegion(center=center, radius=radius)
        region.visual = {'color': color, 'width': width}

        if labels and i < len(labels):
            region.meta = {'text': labels[i]}

        region_list.append(region)

    regions = Regions(region_list)
    regions.write(str(output_path), format='ds9', overwrite=True)


def _write_circle_regions_manual(
    centers: List[tuple],
    radii: List[float],
    output_path: Union[str, Path],
    color: str = 'cyan',
    width: int = 1,
    coordinate_system: str = 'image',
    labels: Optional[List[str]] = None,
) -> None:
    """Fallback manual circle region writer."""
    output_path = Path(output_path)

    lines = [
        '# Region file format: DS9 version 4.1',
        f'global color={color} width={width}',
        f'{coordinate_system}',
    ]

    for i, ((x, y), radius) in enumerate(zip(centers, radii)):
        # DS9 uses 1-indexed coordinates
        x_ds9 = x + 1 if coordinate_system == 'image' else x
        y_ds9 = y + 1 if coordinate_system == 'image' else y

        if labels and i < len(labels):
            lines.append(f'circle({x_ds9:.2f},{y_ds9:.2f},{radius:.2f}) # text={{{labels[i]}}}')
        else:
            lines.append(f'circle({x_ds9:.2f},{y_ds9:.2f},{radius:.2f})')

    output_path.write_text('\n'.join(lines) + '\n')


def read_peak_regions(
    input_path: Union[str, Path],
) -> List[tuple]:
    """
    Read peak positions from DS9 region file.

    Args:
        input_path: Path to .reg file

    Returns:
        List of (x, y) tuples (0-indexed pixel coordinates)
    """
    try:
        from regions import Regions
    except ImportError:
        return _read_peak_regions_manual(input_path)

    input_path = Path(input_path)
    regions = Regions.read(str(input_path), format='ds9')

    positions = []
    for region in regions:
        if hasattr(region, 'center'):
            # Convert to 0-indexed
            positions.append((region.center.x, region.center.y))

    return positions


def _read_peak_regions_manual(
    input_path: Union[str, Path],
) -> List[tuple]:
    """Fallback manual region reader."""
    import re

    input_path = Path(input_path)
    text = input_path.read_text()

    positions = []
    # Match point(x,y) or circle(x,y,r)
    pattern = r'(?:point|circle)\s*\(\s*([\d.]+)\s*,\s*([\d.]+)'
    for match in re.finditer(pattern, text):
        x = float(match.group(1)) - 1  # Convert to 0-indexed
        y = float(match.group(2)) - 1
        positions.append((x, y))

    return positions


def box_vertices(
    x: float,
    y: float,
    width: float,
    height: float,
    angle: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate the vertices of a rotated box.

    This is a Python translation of IDL mask_ds9_box_vertices.pro.

    Args:
        x: X coordinate of box center
        y: Y coordinate of box center
        width: Width of the box
        height: Height of the box
        angle: Rotation angle in radians (counterclockwise)

    Returns:
        Tuple of (x_vertices, y_vertices) arrays, each with 4 elements
        ordered: top-left, top-right, bottom-right, bottom-left
    """
    # Half dimensions
    wr = width / 2.0
    hr = height / 2.0

    cos_a = np.cos(angle)
    sin_a = np.sin(angle)

    # Calculate vertices
    # Top-left from box frame of reference
    xx_tl = x + (-hr * sin_a) + (-wr * cos_a)
    yy_tl = y + (hr * cos_a) + (-wr * sin_a)

    # Top-right from box frame of reference
    xx_tr = x + (-hr * sin_a) + (wr * cos_a)
    yy_tr = y + (hr * cos_a) + (wr * sin_a)

    # Bottom-right from box frame of reference
    xx_br = x + (hr * sin_a) + (wr * cos_a)
    yy_br = y + (-hr * cos_a) + (wr * sin_a)

    # Bottom-left from box frame of reference
    xx_bl = x + (hr * sin_a) + (-wr * cos_a)
    yy_bl = y + (-hr * cos_a) + (-wr * sin_a)

    x_vertices = np.array([xx_tl, xx_tr, xx_br, xx_bl])
    y_vertices = np.array([yy_tl, yy_tr, yy_br, yy_bl])

    return x_vertices, y_vertices


def write_box_regions(
    centers: List[Tuple[float, float]],
    widths: List[float],
    heights: List[float],
    output_path: Union[str, Path],
    angles: Optional[List[float]] = None,
    color: str = 'green',
    width: int = 1,
    coordinate_system: str = 'image',
    labels: Optional[List[str]] = None,
) -> None:
    """
    Write box regions to DS9 file.

    Args:
        centers: List of (x, y) center coordinates
        widths: List of box widths
        heights: List of box heights
        output_path: Path for output .reg file
        angles: List of rotation angles in degrees (0 if not specified)
        color: Region color
        width: Line width
        coordinate_system: Coordinate system
        labels: Optional list of labels
    """
    try:
        from regions import RectanglePixelRegion, PixCoord, Regions
    except ImportError:
        _write_box_regions_manual(
            centers, widths, heights, output_path,
            angles, color, width, coordinate_system, labels
        )
        return

    output_path = Path(output_path)
    region_list = []

    if angles is None:
        angles = [0.0] * len(centers)

    for i, ((x, y), w, h, ang) in enumerate(zip(centers, widths, heights, angles)):
        from astropy import units as u
        center = PixCoord(x=x, y=y)
        region = RectanglePixelRegion(
            center=center,
            width=w,
            height=h,
            angle=ang * u.deg,
        )
        region.visual = {'color': color, 'width': width}

        if labels and i < len(labels):
            region.meta = {'text': labels[i]}

        region_list.append(region)

    regions = Regions(region_list)
    regions.write(str(output_path), format='ds9', overwrite=True)


def _write_box_regions_manual(
    centers: List[Tuple[float, float]],
    widths: List[float],
    heights: List[float],
    output_path: Union[str, Path],
    angles: Optional[List[float]] = None,
    color: str = 'green',
    width: int = 1,
    coordinate_system: str = 'image',
    labels: Optional[List[str]] = None,
) -> None:
    """Fallback manual box region writer."""
    output_path = Path(output_path)

    if angles is None:
        angles = [0.0] * len(centers)

    lines = [
        '# Region file format: DS9 version 4.1',
        f'global color={color} width={width}',
        f'{coordinate_system}',
    ]

    for i, ((x, y), w, h, ang) in enumerate(zip(centers, widths, heights, angles)):
        # DS9 uses 1-indexed coordinates
        x_ds9 = x + 1 if coordinate_system == 'image' else x
        y_ds9 = y + 1 if coordinate_system == 'image' else y

        if labels and i < len(labels):
            lines.append(f'box({x_ds9:.2f},{y_ds9:.2f},{w:.2f},{h:.2f},{ang:.2f}) # text={{{labels[i]}}}')
        else:
            lines.append(f'box({x_ds9:.2f},{y_ds9:.2f},{w:.2f},{h:.2f},{ang:.2f})')

    output_path.write_text('\n'.join(lines) + '\n')


def read_box_regions(
    input_path: Union[str, Path],
) -> List[Tuple[float, float, float, float, float]]:
    """
    Read box regions from DS9 region file.

    Args:
        input_path: Path to .reg file

    Returns:
        List of (x, y, width, height, angle) tuples (0-indexed pixel coordinates)
    """
    try:
        from regions import Regions
    except ImportError:
        return _read_box_regions_manual(input_path)

    input_path = Path(input_path)
    regions = Regions.read(str(input_path), format='ds9')

    boxes = []
    for region in regions:
        if hasattr(region, 'width') and hasattr(region, 'height'):
            angle = float(region.angle.value) if hasattr(region, 'angle') else 0.0
            boxes.append((
                region.center.x,
                region.center.y,
                float(region.width),
                float(region.height),
                angle,
            ))

    return boxes


def _read_box_regions_manual(
    input_path: Union[str, Path],
) -> List[Tuple[float, float, float, float, float]]:
    """Fallback manual box region reader."""
    import re

    input_path = Path(input_path)
    text = input_path.read_text()

    boxes = []
    # Match box(x,y,width,height,angle)
    pattern = r'box\s*\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*(?:,\s*([\d.-]+))?\s*\)'
    for match in re.finditer(pattern, text, re.IGNORECASE):
        x = float(match.group(1)) - 1  # Convert to 0-indexed
        y = float(match.group(2)) - 1
        w = float(match.group(3))
        h = float(match.group(4))
        ang = float(match.group(5)) if match.group(5) else 0.0
        boxes.append((x, y, w, h, ang))

    return boxes
