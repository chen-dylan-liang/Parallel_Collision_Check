use apollo_rust_spatial::vectors::V3;
use super::structs::{AABB, BVHNode, BVHInternalNode, BVHLeafNode};
use std::collections::HashSet;
use std::hash::Hash;

pub fn serial_longest_extent_axis(aabb_indices: &[usize], all_aabbs:&[AABB])->(usize, f64){
    let mut min_v =  V3::new(f64::INFINITY,  f64::INFINITY,  f64::INFINITY);
    let mut max_v = V3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &i in aabb_indices {
        let bb = &all_aabbs[i];
        min_v = min_v.inf(&bb.min_coords);
        max_v = max_v.sup(&bb.max_coords);
    }
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

pub fn serial_split_at_axis<'a>(aabb_indices: &'a mut [usize], all_aabbs:&[AABB], axis:usize, midpoint: f64)->(&'a mut [usize], &'a mut [usize]){
    let mid = aabb_indices.into_iter().partition_in_place(
        |&idx| {
            all_aabbs[idx].center[axis] < midpoint
        }
    );
    aabb_indices.split_at_mut(mid)
}

pub fn serial_build_bvh(aabb_indices: &mut [usize], all_aabbs:&[AABB], cut_off_size:usize)->Box<dyn BVHNode>{
    if aabb_indices.len() <= cut_off_size{
        return Box::new(BVHLeafNode::new( aabb_indices.to_vec(), all_aabbs));
    }
    let (axis,midpoint) = serial_longest_extent_axis(aabb_indices, all_aabbs);
    let (indices_left, indices_right) = serial_split_at_axis(aabb_indices, all_aabbs, axis, midpoint);
    // safeguard for degenerate cases where all bounding boxes equal
    if indices_left.is_empty() || indices_right.is_empty() {
        return Box::new(BVHLeafNode::new(aabb_indices.to_vec(), all_aabbs));
    }
    // do recursive calls
    let left_tree = serial_build_bvh(indices_left, all_aabbs, cut_off_size);
    let right_tree = serial_build_bvh(indices_right, all_aabbs, cut_off_size);
    let aabb = left_tree.union_aabb(&*right_tree);
    Box::new(BVHInternalNode::new(aabb, left_tree, right_tree))
}

pub fn serial_broad_phase_check(
    s1: &dyn BVHNode,
    s2: &dyn BVHNode,
)->Vec<(usize, usize)> {
    // No intersection ⇒ no contacts
    if !s1.intersects(s2) {
        return Vec::new();
    }

    match (s1.is_leaf(), s2.is_leaf()) {
        // both leaves ⇒ collect leaf–leaf pairs
        (true, true) => {
            let is1 = s1.leaf_indices().unwrap();
            let is2 = s2.leaf_indices().unwrap();
            let mut out = Vec::with_capacity(is1.len() * is2.len());
            for &i in is1 {
                for &j in is2 {
                    if i < j {
                        out.push((i, j));
                    }
                }
            }
            out
        }

        // one side is leaf ⇒ recurse on the other side and concatenate
        (true, false) => {
            let (l, r) = s2.children();
            let mut left  =serial_broad_phase_check(s1, l.unwrap());
            let right = serial_broad_phase_check(s1, r.unwrap());
            left.extend(right);
            left
        }
        (false, true) => {
            let (l, r) = s1.children();
            let mut left  = serial_broad_phase_check(l.unwrap(), s2);
            let right = serial_broad_phase_check(r.unwrap(), s2);
            left.extend(right);
            left
        }

        // both internal ⇒ spawn two tasks and merge their results
        (false, false) => {
            let (s1l, s1r) = s1.children();
            let (s2l, s2r) = s2.children();
                let mut v = Vec::new();
                v.extend(serial_broad_phase_check(s1l.unwrap(), s2l.unwrap(), ));
                v.extend(serial_broad_phase_check(s1l.unwrap(), s2r.unwrap(),));
                v.extend(serial_broad_phase_check(s1r.unwrap(), s2l.unwrap(),));
                v.extend(serial_broad_phase_check(s1r.unwrap(), s2r.unwrap(),));
                v
            
        }
    }
}