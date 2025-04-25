use std::cmp::min;
use std::sync::Mutex;
use apollo_rust_spatial::vectors::V3;
use rayon::prelude::*;
use crate::shape::shape::ShapeTrait;

#[derive(Debug, Clone, Copy)]
struct AABB{
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

trait BVHNode: Send + Sync {
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

struct BVHInternalNode{
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

struct BVHLeafNode{
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

// essentially divide-and-conquer in parallel
fn longest_extent_axis(aabb_indices: &[usize], all_aabbs:&[AABB])->(usize, f64){
    let (min_v, max_v) = aabb_indices
        .par_iter()
        .copied()
        .map(|i| {
            let bb = &all_aabbs[i];
            (bb.min_coords, bb.max_coords)
        })
        .reduce(
            // identity element
            || (
                V3::new(f64::INFINITY,  f64::INFINITY,  f64::INFINITY),
                V3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
            ),
            |(min_a, max_a), (min_b, max_b)| (min_a.inf(&min_b), max_a.sup(&max_b)),
        );

    let extent = [
        max_v[0] - min_v[0],
        max_v[1] - min_v[1],
        max_v[2] - min_v[2],
    ];
    let mut axis = 0;
    if extent[1]>=extent[0]&&extent[1]>=extent[2] {axis=1}
    else if extent[2]>=extent[0]&&extent[2]>=extent[1] {axis=2}
    (axis, 0.5*(max_v[axis] + min_v[axis]))
}

// essentially a parallel filter
fn split_at_axis<'a>(aabb_indices: &'a mut [usize], all_aabbs:&[AABB], axis:usize, midpoint: f64)->(&'a mut [usize], &'a mut [usize]){
    let mid = aabb_indices.par_partition_in_place(|&idx| {
        all_aabbs[idx].center[axis] < midpoint
    });
    aabb_indices.split_at_mut(mid)
}

pub fn parallel_build_bvh(aabb_indices: &mut [usize], all_aabbs:&[AABB], cut_off_size:usize)->Box<dyn BVHNode>{
    // S=O(1)
    if aabb_indices.len() <= cut_off_size{
        return Box::new(BVHLeafNode::new( aabb_indices.to_vec(), all_aabbs));
    }
    // S=O(logn)
    let (axis,midpoint) = longest_extent_axis(aabb_indices, all_aabbs);
    // S=O(logn)
    let (indices_left, indices_right) = split_at_axis(aabb_indices, all_aabbs, axis, midpoint);
    // safeguard for degenerate cases where all bounding boxes equal
    if indices_left.is_empty() || indices_right.is_empty() {
        return Box::new(BVHLeafNode::new(aabb_indices.to_vec(), all_aabbs));
    }
    // do recursive calls
    let (left_tree, right_tree)=rayon::join(
        || parallel_build_bvh(indices_left, all_aabbs, cut_off_size),
        || parallel_build_bvh(indices_right, all_aabbs, cut_off_size));
    // S=O(1)
    let aabb = left_tree.union_aabb(&right_tree);
    Box::new(BVHInternalNode::new(aabb, left_tree, right_tree))
}

pub fn parallel_broad_phase_check(
    s1: &dyn BVHNode,
    s2: &dyn BVHNode,
    contacts: &Mutex<Vec<(usize, usize)>>,
) {
    rayon::scope(|scope| {
        scope.spawn(|_| dfs_bvh_pairs(s1, s2, contacts, scope));
    });
}

fn dfs_bvh_pairs<'scope>(
    s1: &dyn BVHNode,
    s2: &dyn BVHNode,
    contacts: &Mutex<Vec<(usize, usize)>>,
    scope: &rayon::Scope<'scope>) {
    if !s1.intersects(s2) {
        return;
    }

    match (s1.is_leaf(), s2.is_leaf()) {
        (true, true) => {
            let is1 = s1.leaf_indices().unwrap();
            let is2 = s2.leaf_indices().unwrap();
            let mut buf = contacts.lock().unwrap();
            for &i in is1 {
                for &j in is2 {
                    if i < j { buf.push((i, j)); }       // skip self/dup
                }
            }
        }

        (true, false) => {
            let (left, right) = s2.children();
            dfs_bvh_pairs(s1, left.unwrap(), contacts, scope);
            dfs_bvh_pairs(s1, right.unwrap(), contacts, scope);
        }
        (false, true) => {
            let (left, right) = s1.children();
            dfs_bvh_pairs(left.unwrap(), s2, contacts, scope);
            dfs_bvh_pairs(right.unwrap(), s2, contacts, scope);
        }

        (false, false) => {
            let (s1_left, s1_right) = s1.children();
            let (s2_left, s2_right) = s2.children();

            // Spawn a new task
            scope.spawn(move |_| {
                dfs_bvh_pairs(s1_left.unwrap(), s2_left.unwrap(), contacts, scope);
                dfs_bvh_pairs(s1_left.unwrap(), s2_right.unwrap(), contacts, scope);
            });
            // Current thread does the other two pairs
            dfs_bvh_pairs(s1_right.unwrap(), s2_left.unwrap(), contacts, scope);
            dfs_bvh_pairs(s1_right.unwrap(), s2_right.unwrap(), contacts, scope);
        }
    }
}