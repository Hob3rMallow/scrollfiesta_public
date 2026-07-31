use std::fmt;

pub const ARENA_MAGIC: u64 = 0x4D4C_5343_5542_4531; // "MLSCUBE1"
pub const ARENA_VERSION: u32 = 1;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u32)]
pub enum Backend {
    Auto = 0,
    RustCpu = 1,
    CubeCpu = 2,
    CubeHip = 3,
    CubeCuda = 4,
}

impl TryFrom<u32> for Backend {
    type Error = ArenaError;

    fn try_from(value: u32) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::Auto),
            1 => Ok(Self::RustCpu),
            2 => Ok(Self::CubeCpu),
            3 => Ok(Self::CubeHip),
            4 => Ok(Self::CubeCuda),
            _ => Err(ArenaError::InvalidBackend(value)),
        }
    }
}

#[derive(Debug)]
pub enum ArenaError {
    EmptyPoints,
    NonFinitePoint { index: usize },
    InvalidCellSize(f32),
    GridTooLarge,
    InvalidGrid(&'static str),
    InvalidBackend(u32),
}

impl fmt::Display for ArenaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyPoints => write!(f, "point cloud is empty"),
            Self::NonFinitePoint { index } => {
                write!(f, "point {index} contains a non-finite coordinate")
            }
            Self::InvalidCellSize(v) => write!(f, "cell size must be finite and > 0, got {v}"),
            Self::GridTooLarge => write!(f, "cell grid is too large for u32 indexing"),
            Self::InvalidGrid(msg) => write!(f, "invalid cell grid: {msg}"),
            Self::InvalidBackend(v) => {
                write!(f, "unknown backend id {v}; expected 0, 1, 2, 3 or 4")
            }
        }
    }
}

impl std::error::Error for ArenaError {}

/// Rust-owned arena consumed by the explicit static-support projection API.
///
/// Points are stored as flat xyz triples. `cell_offsets` is a counting-sort
/// prefix array of length `dims product + 1`, and `point_indices` stores point
/// indices grouped by cell.
#[derive(Debug)]
pub struct MlsArena {
    pub magic: u64,
    pub version: u32,
    pub backend: Backend,
    pub points: Vec<f32>,
    pub cell_offsets: Vec<u32>,
    pub point_indices: Vec<u32>,
    pub dims: [u32; 3],
    pub origin: [f32; 3],
    pub cell_size: f32,
}

impl MlsArena {
    pub fn from_points(
        points: &[[f32; 3]],
        cell_size: f32,
        origin: Option<[f32; 3]>,
    ) -> Result<Self, ArenaError> {
        if points.is_empty() {
            return Err(ArenaError::EmptyPoints);
        }
        if !cell_size.is_finite() || cell_size <= 0.0 {
            return Err(ArenaError::InvalidCellSize(cell_size));
        }
        for (index, p) in points.iter().enumerate() {
            if !p.iter().all(|x| x.is_finite()) {
                return Err(ArenaError::NonFinitePoint { index });
            }
        }

        let mut min = points[0];
        let mut max = points[0];
        for p in &points[1..] {
            for axis in 0..3 {
                min[axis] = min[axis].min(p[axis]);
                max[axis] = max[axis].max(p[axis]);
            }
        }

        // A one-cell halo keeps boundary points away from the outer edge and
        // makes the neighbour-cell walk deterministic.
        let origin = origin.unwrap_or([min[0] - cell_size, min[1] - cell_size, min[2] - cell_size]);
        if !origin.iter().all(|x| x.is_finite()) {
            return Err(ArenaError::InvalidGrid(
                "origin contains a non-finite value",
            ));
        }
        if (0..3).any(|axis| origin[axis] > min[axis]) {
            return Err(ArenaError::InvalidGrid(
                "origin must not lie above the minimum point coordinate",
            ));
        }

        let mut dims = [0_u32; 3];
        for axis in 0..3 {
            let span = ((max[axis] - origin[axis]) / cell_size).ceil();
            if !span.is_finite() || span < 0.0 || span > (u32::MAX - 2) as f32 {
                return Err(ArenaError::GridTooLarge);
            }
            dims[axis] = span as u32 + 2;
        }

        let ncells_u64 = u64::from(dims[0])
            .checked_mul(u64::from(dims[1]))
            .and_then(|value| value.checked_mul(u64::from(dims[2])))
            .ok_or(ArenaError::GridTooLarge)?;
        if ncells_u64 == 0 || ncells_u64 > (u32::MAX - 1) as u64 {
            return Err(ArenaError::GridTooLarge);
        }
        let ncells = ncells_u64 as usize;

        let mut counts = vec![0_u32; ncells];
        let mut point_cells = Vec::with_capacity(points.len());
        for p in points {
            let cell = point_cell(*p, origin, cell_size, dims);
            counts[cell] = counts[cell]
                .checked_add(1)
                .ok_or(ArenaError::GridTooLarge)?;
            point_cells.push(cell);
        }

        let mut offsets = vec![0_u32; ncells + 1];
        for i in 0..ncells {
            offsets[i + 1] = offsets[i]
                .checked_add(counts[i])
                .ok_or(ArenaError::GridTooLarge)?;
        }
        if offsets[ncells] as usize != points.len() {
            return Err(ArenaError::GridTooLarge);
        }

        let mut cursors = offsets[..ncells].to_vec();
        let mut point_indices = vec![0_u32; points.len()];
        for (point_index, &cell) in point_cells.iter().enumerate() {
            let slot = cursors[cell] as usize;
            point_indices[slot] =
                u32::try_from(point_index).map_err(|_| ArenaError::GridTooLarge)?;
            cursors[cell] += 1;
        }

        // Match ScrollFiesta's deterministic within-cell support order:
        // position z, then y, then x, then original index. This matters for
        // reproducible FP32 accumulation near tied eigenspaces.
        for cell in 0..ncells {
            let begin = offsets[cell] as usize;
            let end = offsets[cell + 1] as usize;
            point_indices[begin..end].sort_unstable_by(|left, right| {
                let a = points[*left as usize];
                let b = points[*right as usize];
                a[2].total_cmp(&b[2])
                    .then_with(|| a[1].total_cmp(&b[1]))
                    .then_with(|| a[0].total_cmp(&b[0]))
                    .then_with(|| left.cmp(right))
            });
        }

        let mut flat = Vec::with_capacity(points.len() * 3);
        for p in points {
            flat.extend_from_slice(p);
        }

        Ok(Self {
            magic: ARENA_MAGIC,
            version: ARENA_VERSION,
            backend: Backend::Auto,
            points: flat,
            cell_offsets: offsets,
            point_indices,
            dims,
            origin,
            cell_size,
        })
    }

    pub fn from_grid(
        points: Vec<f32>,
        cell_offsets: Vec<u32>,
        point_indices: Vec<u32>,
        dims: [u32; 3],
        origin: [f32; 3],
        cell_size: f32,
    ) -> Result<Self, ArenaError> {
        if points.is_empty() || !points.len().is_multiple_of(3) {
            return Err(ArenaError::InvalidGrid(
                "points must contain one or more xyz triples",
            ));
        }
        if !points.iter().all(|x| x.is_finite()) {
            return Err(ArenaError::InvalidGrid("points contain a non-finite value"));
        }
        if !cell_size.is_finite() || cell_size <= 0.0 {
            return Err(ArenaError::InvalidCellSize(cell_size));
        }
        if !origin.iter().all(|x| x.is_finite()) || dims.contains(&0) {
            return Err(ArenaError::InvalidGrid("origin/dimensions are invalid"));
        }
        let ncells_u64 = u64::from(dims[0])
            .checked_mul(u64::from(dims[1]))
            .and_then(|value| value.checked_mul(u64::from(dims[2])))
            .ok_or(ArenaError::GridTooLarge)?;
        let ncells = usize::try_from(ncells_u64).map_err(|_| ArenaError::GridTooLarge)?;
        let offsets_len = ncells.checked_add(1).ok_or(ArenaError::GridTooLarge)?;
        if cell_offsets.len() != offsets_len {
            return Err(ArenaError::InvalidGrid(
                "cell_offsets length must equal cell count + 1",
            ));
        }
        if cell_offsets.first() != Some(&0) {
            return Err(ArenaError::InvalidGrid("cell_offsets must start at zero"));
        }
        if !cell_offsets.windows(2).all(|w| w[0] <= w[1]) {
            return Err(ArenaError::InvalidGrid("cell_offsets must be monotonic"));
        }
        if cell_offsets.last().copied().unwrap_or_default() as usize != point_indices.len() {
            return Err(ArenaError::InvalidGrid(
                "last cell offset must equal point_indices length",
            ));
        }
        let npoints = points.len() / 3;
        if point_indices.len() != npoints {
            return Err(ArenaError::InvalidGrid(
                "point_indices must contain every point exactly once",
            ));
        }
        if point_indices.iter().any(|&i| i as usize >= npoints) {
            return Err(ArenaError::InvalidGrid(
                "point_indices contains an out-of-range index",
            ));
        }
        let mut seen = vec![false; npoints];
        for &index in &point_indices {
            let slot = index as usize;
            if seen[slot] {
                return Err(ArenaError::InvalidGrid(
                    "point_indices contains a duplicate index",
                ));
            }
            seen[slot] = true;
        }

        Ok(Self {
            magic: ARENA_MAGIC,
            version: ARENA_VERSION,
            backend: Backend::Auto,
            points,
            cell_offsets,
            point_indices,
            dims,
            origin,
            cell_size,
        })
    }

    #[inline]
    pub fn npoints(&self) -> usize {
        self.points.len() / 3
    }

    #[inline]
    pub fn validate_tag(&self) -> bool {
        self.magic == ARENA_MAGIC && self.version == ARENA_VERSION
    }
}

#[inline]
pub fn linear_cell(x: u32, y: u32, z: u32, dims: [u32; 3]) -> usize {
    (x as usize) + (dims[0] as usize) * ((y as usize) + (dims[1] as usize) * (z as usize))
}

#[inline]
pub fn point_cell_coords(
    p: [f32; 3],
    origin: [f32; 3],
    cell_size: f32,
    dims: [u32; 3],
) -> [u32; 3] {
    let mut out = [0_u32; 3];
    for axis in 0..3 {
        let raw = ((p[axis] - origin[axis]) / cell_size).floor();
        out[axis] = if raw <= 0.0 {
            0
        } else if raw >= (dims[axis] - 1) as f32 {
            dims[axis] - 1
        } else {
            raw as u32
        };
    }
    out
}

#[inline]
pub fn point_cell(p: [f32; 3], origin: [f32; 3], cell_size: f32, dims: [u32; 3]) -> usize {
    let c = point_cell_coords(p, origin, cell_size, dims);
    linear_cell(c[0], c[1], c[2], dims)
}
