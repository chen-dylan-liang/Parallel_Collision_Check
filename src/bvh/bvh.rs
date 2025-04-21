use std::cmp::min;
use std::sync::Mutex;
use apollo_rust_spatial::vectors::V3;
use rayon::prelude::*;
use crate::shape::shape::ShapeTrait;

#[derive(Clone)]
struct AABB{
    pub min_coords: V3,
    pub max_coords: V3,
    pub center: V3,
}

impl AABB {
    pub fn new(min_coords: V3, max_coords:V3)->Self {
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
        let aabb = self.clone();
        for i in 0..3{
            if aabb.min_coords[i]<other.min_coords[i] {aabb.min_coords[i]=other.min_coords[i]}
            if aabb.max_coords[i]>other.max_coords[i] {aabb.max_coords[i]=other.max_coords[i]}
        }
        aabb
    }
}

trait BVHNode{
    fn is_leaf(&self) -> bool;

    fn children(&self) -> (Option<&Box<dyn BVHNode>>, Option<&Box<dyn BVHNode>>);

    fn leaf_indices(&self) -> Option<&[usize]>;

    fn aabb_ref(&self)->&AABB;

    fn union_aabb(&self, other: &Box<dyn BVHNode>) -> AABB{
        self.aabb_ref().union(other.aabb_ref())
    }

    fn intersects(&self, other: &Box<dyn BVHNode>) -> bool {
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
    pub fn new(aabb:AABB, shape_indices: Vec<usize>)->Self{
        Self{
            aabb,
            shape_indices,
        }
    }
}

impl BVHNode for BVHInternalNode{
    fn is_leaf(&self) -> bool{false}

    fn children(&self) -> (Option<&Box<dyn BVHNode>>, Option<&Box<dyn BVHNode>>) {
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

    fn children(&self) -> (Option<&Box<dyn BVHNode>>, Option<&Box<dyn BVHNode>>) {
        (None, None)
    }

    fn leaf_indices(&self) -> Option<&[usize]> {
        Some(self.shape_indices.as_slice())
    }

    fn aabb_ref(&self) -> &AABB {
        &self.aabb
    }
}

fn bvh_leaf_aabb(aabb_indices: &[usize], all_aabbs:&Vec<AABB>)->AABB{
   let mut min_coords = V3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max_coords = V3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
   for aabb_index in aabb_indices{
      for i in 0..3{
          if all_aabbs[aabb_index].min_coords[i]<min_coords[i] {min_coords[i]=all_aabbs[aabb_index].min_coords[i]}
          if all_aabbs[aabb_index].max_coords[i]>max_coords[i] {max_coords[i]=all_aabbs[aabb_index].max_coords[i]}
      }
   }
    AABB::new(min_coords, max_coords)
}

// essentially divide-and-conquer in parallel
fn longest_extent_axis(aabb_indices: &[usize], all_aabbs:&Vec<AABB>)->(usize, f64){
    let (min_v, max_v) = aabb_indices
        .par_iter()                       // iterate indices in parallel
        .map(|&i| {
            let bb = &all_aabbs[i];
            (bb.min_coords, bb.max_coords)
        })
        .reduce(
            || (V3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
                V3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY)),
            |(min_a, max_a), (min_b, max_b)| {
                (
                    V3::new( min_a[0].min(min_b[0]),
                        min_a[1].min(min_b[1]),
                        min_a[2].min(min_b[2])),
                    V3::new( max_a[0].max(max_b[0]),
                        max_a[1].max(max_b[1]),
                        max_a[2].max(max_b[2])),
                )
            },
        );

    let extent = [
        max_v[0] - min_v[0],
        max_v[1] - min_v[1],
        max_v[2] - min_v[2],
    ];
    let mut axis = 0;
    if extent[1]>extent[0]&&extent[1]>extent[2] {axis=1}
    else if extent[2]>extent[0]&&extent[2]>extent[1] {axis=2}
    (axis, 0.5*extent[axis])
}

// essentially a parallel filter
fn split_at_axis<'a>(aabb_indices: &'a mut [usize], all_aabbs:&Vec<AABB>, axis:usize, midpoint: f64)->(&'a mut [usize], &'a mut [usize]){
    let (left, right): (Vec<usize>, Vec<usize>) =
        aabb_indices
            .par_iter()
            .cloned()
            .partition(|&idx| all_aabbs[idx].center[axis] < midpoint);

    let left_len = left.len();
    aabb_indices[..left_len].copy_from_slice(&left);
    aabb_indices[left_len..].copy_from_slice(&right);
    let (left_slice, right_slice) = aabb_indices.split_at_mut(left_len);
    (left_slice, right_slice)
}

pub fn parallel_build_bvh(aabb_indices: &mut [usize], all_aabbs:&Vec<AABB>, cut_off_size:usize)->Box<dyn BVHNode>{
    // S=O(1)
    if aabb_indices.len() <= cut_off_size{
        let aabb = bvh_leaf_aabb(aabb_indices, all_aabbs);
        return Box::new(BVHLeafNode::new(aabb, aabb_indices.to_vec()));
    }
    // S=O(logn)
    let (axis,midpoint) = longest_extent_axis(aabb_indices, all_aabbs);
    // S=O(logn)
    let (indices_left, indices_right) = split_at_axis(aabb_indices, all_aabbs, axis, midpoint);
    // do recursive calls
    let (left_tree, right_tree)=rayon::join(
        || parallel_build_bvh(indices_left, all_aabbs, cut_off_size),
        || parallel_build_bvh(indices_right, all_aabbs, cut_off_size));
    // S=O(1)
    let aabb = left_tree.union_aabb(&right_tree);
    Box::new(BVHInternalNode::new(aabb, left_tree, right_tree))
}

pub fn parallel_broad_phase_check(
    s1: &Box<dyn BVHNode>,
    s2: &Box<dyn BVHNode>,
    contacts: &Mutex<Vec<(usize, usize)>>,
) {
    rayon::scope(|scope| {
        scope.spawn(|_| dfs_bvh_pairs(s1, s2, contacts, scope));
    });
}

fn dfs_bvh_pairs<'scope>(
    s1: &Box<dyn BVHNode>,
    s2: &Box<dyn BVHNode>,
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