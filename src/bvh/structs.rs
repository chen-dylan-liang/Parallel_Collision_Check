use apollo_rust_spatial::vectors::V3;

#[derive(Debug, Clone, Copy)]
pub struct AABB{
    pub min_coords: V3,
    pub max_coords: V3,
    pub center: V3,
}

impl AABB {
    pub fn new(min_coords: V3, max_coords:V3)->Self {
        debug_assert!(
            (0..3).all(|i| min_coords[i] <= max_coords[i]),
            "Invalid AABB: min > max"
        );
        Self { min_coords, max_coords, center:0.5*(min_coords + max_coords) }
    }

    pub fn intersects(&self, other: &AABB) -> bool {
        // Separate on X axis
        if self.max_coords.x < other.min_coords.x || other.max_coords.x < self.min_coords.x {
            return false;
        }
        // Separate on Y axis
        if self.max_coords.y < other.min_coords.y || other.max_coords.y < self.min_coords.y {
            return false;
        }
        // Separate on Z axis
        if self.max_coords.z < other.min_coords.z || other.max_coords.z < self.min_coords.z {
            return false;
        }
        true
    }

    pub fn union(&self, other: &AABB) -> AABB {
        let min = self.min_coords.inf(&other.min_coords);
        let max = self.max_coords.sup(&other.max_coords);
        AABB::new(min, max)
    }
}

pub trait BVHNode: Send + Sync {
    fn is_leaf(&self) -> bool;

    fn children(&self) -> (Option<&dyn BVHNode>, Option<&dyn BVHNode>);

    fn leaf_indices(&self) -> Option<&[usize]>;

    fn aabb_ref(&self)->&AABB;

    fn union_aabb(&self, other: &dyn BVHNode) -> AABB{
        self.aabb_ref().union(other.aabb_ref())
    }

    fn intersects(&self, other: &dyn BVHNode) -> bool {
        self.aabb_ref().intersects(other.aabb_ref())
    }
}

pub struct BVHInternalNode{
    pub aabb: AABB,
    pub left: Box<dyn BVHNode>,
    pub right: Box<dyn BVHNode>,
}

impl BVHInternalNode{
    pub fn new(aabb: AABB,  left: Box<dyn BVHNode>, right:  Box<dyn BVHNode>)->Self{
        Self{
            aabb,
            left,
            right
        }
    }
}

pub struct BVHLeafNode{
    pub aabb: AABB,
    pub shape_indices: Vec<usize>,
}

impl BVHLeafNode{
    pub fn new(shape_indices: Vec<usize>, all_aabbs:&[AABB])->Self{
        let mut min_coords = V3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut max_coords = V3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for i in &shape_indices{
            min_coords = min_coords.inf(&all_aabbs[i].min_coords);
            max_coords = max_coords.sup(&all_aabbs[i].max_coords);
        }
        Self{
            shape_indices,
            aabb:  AABB::new(min_coords, max_coords),
        }
    }
}

impl BVHNode for BVHInternalNode{
    fn is_leaf(&self) -> bool{false}

    fn children(&self) -> (Option<&dyn BVHNode>, Option<&dyn BVHNode>) {
        (Some(&self.left), Some(&self.right))
    }

    fn leaf_indices(&self) -> Option<&[usize]> {
        None
    }

    fn aabb_ref(&self) -> &AABB {
        &self.aabb
    }
}

impl BVHNode for BVHLeafNode{
    fn is_leaf(&self) -> bool{true}

    fn children(&self) -> (Option<&dyn BVHNode>, Option<&dyn BVHNode>) {
        (None, None)
    }

    fn leaf_indices(&self) -> Option<&[usize]> {
        Some(self.shape_indices.as_slice())
    }

    fn aabb_ref(&self) -> &AABB {
        &self.aabb
    }
}
