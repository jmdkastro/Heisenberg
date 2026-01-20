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


# =============================================================================
# Full DS9 Region File Parsing (mask_ds9_file_mask.pro translation)
# =============================================================================

from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict


class RegionType(Enum):
    """DS9 region types."""
    CIRCLE = 'circle'
    ELLIPSE = 'ellipse'
    BOX = 'box'
    POLYGON = 'polygon'
    POINT = 'point'
    LINE = 'line'
    ANNULUS = 'annulus'


class CoordinateSystem(Enum):
    """DS9 coordinate systems."""
    IMAGE = 'image'
    PHYSICAL = 'physical'
    FK4 = 'fk4'
    FK5 = 'fk5'
    J2000 = 'j2000'
    B1950 = 'b1950'
    ICRS = 'icrs'
    GALACTIC = 'galactic'
    ECLIPTIC = 'ecliptic'


@dataclass
class DS9Region:
    """
    Parsed DS9 region.

    Attributes:
        region_type: Type of region (circle, ellipse, box, polygon)
        params: Dictionary of region parameters
        coordinate_system: Coordinate system of the region
        is_exclude: If True, region is excluded (prefixed with -)
    """
    region_type: RegionType
    params: Dict[str, Any]
    coordinate_system: CoordinateSystem = CoordinateSystem.IMAGE
    is_exclude: bool = False


def parse_ds9_region_file(
    filepath: Union[str, Path],
    require_image_coords: bool = True,
) -> Tuple[List[DS9Region], CoordinateSystem]:
    """
    Parse a DS9 region file.

    This is a Python translation of mask_ds9_file_mask.pro parsing logic.

    Args:
        filepath: Path to DS9 region file
        require_image_coords: If True, raise error for non-image coordinates

    Returns:
        Tuple of (list of DS9Region objects, coordinate system)

    Raises:
        ValueError: If non-image coordinates detected and require_image_coords=True
    """
    import re

    filepath = Path(filepath)
    text = filepath.read_text()

    regions = []
    current_coord_system = CoordinateSystem.IMAGE

    for line in text.splitlines():
        line = line.strip()

        # Skip empty lines
        if not line:
            continue

        # Handle semicolon-separated definitions (DS9 allows ; as line break)
        for part in line.split(';'):
            part = part.strip()
            if not part:
                continue

            # Skip comments
            if part.startswith('#'):
                # Check for composite region warning
                if 'composite' in part.lower():
                    continue
                continue

            # Remove trailing comments
            if '#' in part:
                part = part.split('#')[0].strip()

            # Remove composite markers (||)
            if '||' in part:
                part = part.split('||')[0].strip()

            part_lower = part.lower()

            # Handle coordinate system declarations
            if part_lower == 'image':
                current_coord_system = CoordinateSystem.IMAGE
                continue
            elif part_lower in ('fk4', 'b1950'):
                current_coord_system = CoordinateSystem.FK4
                if require_image_coords:
                    raise ValueError(
                        f"DS9 file contains {part} coordinates. "
                        "Only IMAGE coordinates are supported. "
                        "Convert using ds9_convert_to_image()."
                    )
                continue
            elif part_lower in ('fk5', 'j2000'):
                current_coord_system = CoordinateSystem.FK5
                if require_image_coords:
                    raise ValueError(
                        f"DS9 file contains {part} coordinates. "
                        "Only IMAGE coordinates are supported. "
                        "Convert using ds9_convert_to_image()."
                    )
                continue
            elif part_lower == 'icrs':
                current_coord_system = CoordinateSystem.ICRS
                if require_image_coords:
                    raise ValueError(
                        "DS9 file contains ICRS coordinates. "
                        "Only IMAGE coordinates are supported. "
                        "Convert using ds9_convert_to_image()."
                    )
                continue
            elif part_lower == 'galactic':
                current_coord_system = CoordinateSystem.GALACTIC
                if require_image_coords:
                    raise ValueError(
                        "DS9 file contains GALACTIC coordinates. "
                        "Only IMAGE coordinates are supported. "
                        "Convert using ds9_convert_to_image()."
                    )
                continue
            elif part_lower == 'physical':
                current_coord_system = CoordinateSystem.PHYSICAL
                if require_image_coords:
                    raise ValueError(
                        "DS9 file contains PHYSICAL coordinates. "
                        "Only IMAGE coordinates are supported. "
                        "Convert using ds9_convert_to_image()."
                    )
                continue
            elif part_lower.startswith('global'):
                continue

            # Check for exclusion prefix
            is_exclude = part.startswith('-')
            if is_exclude:
                part = part[1:].strip()

            # Parse region definitions
            region = _parse_region_definition(part, current_coord_system, is_exclude)
            if region is not None:
                regions.append(region)

    return regions, current_coord_system


def _parse_region_definition(
    definition: str,
    coord_system: CoordinateSystem,
    is_exclude: bool,
) -> Optional[DS9Region]:
    """Parse a single region definition."""
    import re

    definition = definition.strip()
    if not definition:
        return None

    # Extract region type and parameters
    match = re.match(r'(\w+)\s*\((.*)\)', definition, re.IGNORECASE)
    if not match:
        return None

    region_type_str = match.group(1).lower()
    params_str = match.group(2)

    # Parse parameters (comma-separated, may have units)
    params_raw = [p.strip() for p in params_str.split(',')]

    # Remove units and convert to float
    def parse_value(val: str) -> float:
        """Parse a value, removing any unit suffix."""
        val = val.strip()
        # Remove common units: ", ', d, r, p, i
        val = re.sub(r'["\'\s]*$', '', val)
        val = re.sub(r'[dpri]$', '', val, flags=re.IGNORECASE)
        try:
            return float(val)
        except ValueError:
            # May be sexagesimal (HH:MM:SS or DD:MM:SS)
            if ':' in val:
                parts = val.split(':')
                if len(parts) == 3:
                    deg = float(parts[0])
                    sign = -1 if deg < 0 or val.startswith('-') else 1
                    return sign * (abs(deg) + float(parts[1]) / 60 + float(parts[2]) / 3600)
            return 0.0

    params = [parse_value(p) for p in params_raw]

    # Handle different region types
    if region_type_str == 'circle':
        if len(params) >= 3:
            return DS9Region(
                region_type=RegionType.CIRCLE,
                params={'x': params[0], 'y': params[1], 'radius': params[2]},
                coordinate_system=coord_system,
                is_exclude=is_exclude,
            )

    elif region_type_str == 'ellipse':
        if len(params) >= 5:
            return DS9Region(
                region_type=RegionType.ELLIPSE,
                params={
                    'x': params[0],
                    'y': params[1],
                    'semi_major': params[2],
                    'semi_minor': params[3],
                    'angle': params[4],
                },
                coordinate_system=coord_system,
                is_exclude=is_exclude,
            )

    elif region_type_str == 'box':
        if len(params) >= 5:
            return DS9Region(
                region_type=RegionType.BOX,
                params={
                    'x': params[0],
                    'y': params[1],
                    'width': params[2],
                    'height': params[3],
                    'angle': params[4],
                },
                coordinate_system=coord_system,
                is_exclude=is_exclude,
            )
        elif len(params) >= 4:
            return DS9Region(
                region_type=RegionType.BOX,
                params={
                    'x': params[0],
                    'y': params[1],
                    'width': params[2],
                    'height': params[3],
                    'angle': 0.0,
                },
                coordinate_system=coord_system,
                is_exclude=is_exclude,
            )

    elif region_type_str == 'polygon':
        if len(params) >= 6:  # At least 3 vertices
            x_coords = params[0::2]  # Odd indices (0, 2, 4, ...)
            y_coords = params[1::2]  # Even indices (1, 3, 5, ...)
            return DS9Region(
                region_type=RegionType.POLYGON,
                params={'x': x_coords, 'y': y_coords},
                coordinate_system=coord_system,
                is_exclude=is_exclude,
            )

    elif region_type_str == 'point':
        if len(params) >= 2:
            return DS9Region(
                region_type=RegionType.POINT,
                params={'x': params[0], 'y': params[1]},
                coordinate_system=coord_system,
                is_exclude=is_exclude,
            )

    elif region_type_str in ('line', 'vector', 'text', 'ruler', 'compass', 'projection'):
        # Non-masking regions - skip
        return None

    elif region_type_str == 'annulus':
        raise ValueError("Annulus regions are not yet supported")

    elif region_type_str in ('panda', 'epanda', 'bpanda'):
        raise ValueError(f"{region_type_str} regions are not yet supported")

    return None


def ds9_convert_to_image(
    image_path: Union[str, Path],
    region_inpath: Union[str, Path],
    region_outpath: Union[str, Path],
) -> bool:
    """
    Convert DS9 region file to image coordinates using DS9.

    This is a Python translation of mask_ds9_file_convert.pro.

    Args:
        image_path: Path to FITS image file
        region_inpath: Path to input DS9 region file
        region_outpath: Path for output region file in image coordinates

    Returns:
        True if conversion successful, False otherwise

    Raises:
        FileNotFoundError: If DS9 is not installed
    """
    import subprocess
    import shutil

    # Check if DS9 is available
    if shutil.which('ds9') is None:
        raise FileNotFoundError(
            "DS9 not found in PATH. Install DS9 or convert regions manually."
        )

    image_path = Path(image_path)
    region_inpath = Path(region_inpath)
    region_outpath = Path(region_outpath)

    # Build DS9 command
    cmd = [
        'ds9',
        str(image_path),
        '-regions', 'load', str(region_inpath),
        '-regions', 'system', 'image',
        '-scale', 'log',
        '-regions', 'save', str(region_outpath),
        '-exit',
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60,
        )
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        return False
    except Exception:
        return False


def write_ellipse_regions(
    centers: List[Tuple[float, float]],
    semi_majors: List[float],
    semi_minors: List[float],
    output_path: Union[str, Path],
    angles: Optional[List[float]] = None,
    color: str = 'magenta',
    width: int = 1,
    coordinate_system: str = 'image',
    labels: Optional[List[str]] = None,
) -> None:
    """
    Write ellipse regions to DS9 file.

    Args:
        centers: List of (x, y) center coordinates
        semi_majors: List of semi-major axes
        semi_minors: List of semi-minor axes
        output_path: Path for output .reg file
        angles: List of rotation angles in degrees (0 if not specified)
        color: Region color
        width: Line width
        coordinate_system: Coordinate system
        labels: Optional list of labels
    """
    output_path = Path(output_path)

    if angles is None:
        angles = [0.0] * len(centers)

    lines = [
        '# Region file format: DS9 version 4.1',
        f'global color={color} width={width}',
        f'{coordinate_system}',
    ]

    for i, ((x, y), a, b, ang) in enumerate(zip(centers, semi_majors, semi_minors, angles)):
        # DS9 uses 1-indexed coordinates
        x_ds9 = x + 1 if coordinate_system == 'image' else x
        y_ds9 = y + 1 if coordinate_system == 'image' else y

        if labels and i < len(labels):
            lines.append(f'ellipse({x_ds9:.2f},{y_ds9:.2f},{a:.2f},{b:.2f},{ang:.2f}) # text={{{labels[i]}}}')
        else:
            lines.append(f'ellipse({x_ds9:.2f},{y_ds9:.2f},{a:.2f},{b:.2f},{ang:.2f})')

    output_path.write_text('\n'.join(lines) + '\n')


def write_polygon_regions(
    vertices_list: List[List[Tuple[float, float]]],
    output_path: Union[str, Path],
    color: str = 'yellow',
    width: int = 1,
    coordinate_system: str = 'image',
    labels: Optional[List[str]] = None,
) -> None:
    """
    Write polygon regions to DS9 file.

    Args:
        vertices_list: List of polygons, each a list of (x, y) vertices
        output_path: Path for output .reg file
        color: Region color
        width: Line width
        coordinate_system: Coordinate system
        labels: Optional list of labels
    """
    output_path = Path(output_path)

    lines = [
        '# Region file format: DS9 version 4.1',
        f'global color={color} width={width}',
        f'{coordinate_system}',
    ]

    for i, vertices in enumerate(vertices_list):
        # DS9 uses 1-indexed coordinates
        if coordinate_system == 'image':
            coords = ','.join(f'{x+1:.2f},{y+1:.2f}' for x, y in vertices)
        else:
            coords = ','.join(f'{x:.2f},{y:.2f}' for x, y in vertices)

        if labels and i < len(labels):
            lines.append(f'polygon({coords}) # text={{{labels[i]}}}')
        else:
            lines.append(f'polygon({coords})')

    output_path.write_text('\n'.join(lines) + '\n')


def read_ellipse_regions(
    input_path: Union[str, Path],
) -> List[Tuple[float, float, float, float, float]]:
    """
    Read ellipse regions from DS9 region file.

    Args:
        input_path: Path to .reg file

    Returns:
        List of (x, y, semi_major, semi_minor, angle) tuples (0-indexed)
    """
    import re

    input_path = Path(input_path)
    text = input_path.read_text()

    ellipses = []
    pattern = r'ellipse\s*\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*(?:,\s*([\d.-]+))?\s*\)'
    for match in re.finditer(pattern, text, re.IGNORECASE):
        x = float(match.group(1)) - 1
        y = float(match.group(2)) - 1
        a = float(match.group(3))
        b = float(match.group(4))
        ang = float(match.group(5)) if match.group(5) else 0.0
        ellipses.append((x, y, a, b, ang))

    return ellipses


def read_polygon_regions(
    input_path: Union[str, Path],
) -> List[List[Tuple[float, float]]]:
    """
    Read polygon regions from DS9 region file.

    Args:
        input_path: Path to .reg file

    Returns:
        List of polygons, each a list of (x, y) vertices (0-indexed)
    """
    import re

    input_path = Path(input_path)
    text = input_path.read_text()

    polygons = []
    pattern = r'polygon\s*\(([^)]+)\)'
    for match in re.finditer(pattern, text, re.IGNORECASE):
        coords_str = match.group(1)
        coords = [float(c.strip()) for c in coords_str.split(',')]
        # Convert to list of (x, y) tuples, 0-indexed
        vertices = [(coords[i] - 1, coords[i+1] - 1) for i in range(0, len(coords), 2)]
        polygons.append(vertices)

    return polygons


def read_circle_regions(
    input_path: Union[str, Path],
) -> List[Tuple[float, float, float]]:
    """
    Read circle regions from DS9 region file.

    Args:
        input_path: Path to .reg file

    Returns:
        List of (x, y, radius) tuples (0-indexed)
    """
    import re

    input_path = Path(input_path)
    text = input_path.read_text()

    circles = []
    pattern = r'circle\s*\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*\)'
    for match in re.finditer(pattern, text, re.IGNORECASE):
        x = float(match.group(1)) - 1
        y = float(match.group(2)) - 1
        r = float(match.group(3))
        circles.append((x, y, r))

    return circles
