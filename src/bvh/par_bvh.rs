use std::sync::Mutex;
use apollo_rust_spatial::vectors::V3;
use rayon::prelude::*;
use super::structs::{AABB, BVHNode, BVHInternalNode, BVHLeafNode};
use rayon::slice::ParallelSlice;
use rayon::slice::ParallelSliceMut;
use std::collections::HashSet;
use crate::bvh::srl_bvh::{serial_longest_extent_axis, serial_split_at_axis};
const MAX_DEPTH: usize = 4;

// essentially divide-and-conquer in parallel
fn parallel_longest_extent_axis(aabb_indices: &[usize], all_aabbs:&[AABB])->(usize, f64){
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
fn parallel_split_at_axis<'a>(aabb_indices: &'a mut [usize], all_aabbs:&[AABB], axis:usize, midpoint: f64)->(&'a mut [usize], &'a mut [usize]){
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

const BUILD_PARALLEL_THRESHOLD: usize = 1000;

pub fn parallel_build_bvh(
    aabb_indices: &mut [usize],
    all_aabbs:      &[AABB],
    cut_off_size:   usize,
) -> Box<dyn BVHNode> {
    // 1) check size up front
    let n = aabb_indices.len();
    if n <= cut_off_size {
        return Box::new(BVHLeafNode::new(aabb_indices.to_vec(), all_aabbs));
    }

    // decide in advance whether we’ll parallelize
    let do_parallel = n > BUILD_PARALLEL_THRESHOLD;

    // 2) split path by size
    if do_parallel {
        // – compute axis/median in parallel
        let (axis, midpoint) =
            parallel_longest_extent_axis(aabb_indices, all_aabbs);

        // – this mutably borrows `aabb_indices`
        let (left, right) =
            parallel_split_at_axis(aabb_indices, all_aabbs, axis, midpoint);
        if left.is_empty() || right.is_empty() {
            return Box::new(BVHLeafNode::new(aabb_indices.to_vec(), all_aabbs));
        }
        // now spawn the two big recursive tasks
        let (l, r) = rayon::join(
            || parallel_build_bvh(left,  all_aabbs, cut_off_size),
            || parallel_build_bvh(right, all_aabbs, cut_off_size),
        );
        let node_aabb = l.union_aabb(&*r);
        Box::new(BVHInternalNode::new(node_aabb, l, r))
    } else {
        // serial fallback: no mutable/immutable conflict
        let (axis, midpoint) =
            serial_longest_extent_axis(aabb_indices, all_aabbs);
        let (left, right) =
            serial_split_at_axis(aabb_indices, all_aabbs, axis, midpoint);
        if left.is_empty() || right.is_empty() {
            return Box::new(BVHLeafNode::new(aabb_indices.to_vec(), all_aabbs));
        }
        let l = parallel_build_bvh(left,  all_aabbs, cut_off_size);
        let r = parallel_build_bvh(right, all_aabbs, cut_off_size);
        let node_aabb = l.union_aabb(&*r);
        Box::new(BVHInternalNode::new(node_aabb, l, r))
    }
}

pub fn parallel_broad_phase_check(
    s1: &dyn BVHNode,
    s2: &dyn BVHNode,
    contacts: &Mutex<Vec<(usize, usize)>>,
) {
    parallel_dfs_bvh_pairs(s1, s2, contacts, 0);
}

fn parallel_dfs_bvh_pairs(
    s1: &dyn BVHNode ,
    s2: &dyn BVHNode,
    contacts: & Mutex<Vec<(usize, usize)>>,
    depth: usize
) {
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
                    if i<j{
                        buf.push((i,j));
                    }
                }
            }
        }

        (true, false) => {
            let (left, right) = s2.children();
            parallel_dfs_bvh_pairs(s1, left.unwrap(), contacts, depth + 1);
            parallel_dfs_bvh_pairs(s1, right.unwrap(), contacts, depth + 1);
        }
        (false, true) => {
            let (left, right) = s1.children();
            parallel_dfs_bvh_pairs(left.unwrap(), s2, contacts, depth + 1);
            parallel_dfs_bvh_pairs(right.unwrap(), s2, contacts, depth+1);
        }

        (false, false) => {
            let (s1_left, s1_right) = s1.children();
            let (s2_left, s2_right) = s2.children();

            // Create local variables to avoid multiple unwraps
            let s1l = s1_left.unwrap();
            let s1r = s1_right.unwrap();
            let s2l = s2_left.unwrap();
            let s2r = s2_right.unwrap();

            if depth < MAX_DEPTH {
            // Spawn one pair in a new task
                rayon::join(
                    || {
                        parallel_dfs_bvh_pairs(s1l, s2l, contacts, depth+1);
                        parallel_dfs_bvh_pairs(s1l, s2r, contacts, depth+1);
                    },
                    || {
                        parallel_dfs_bvh_pairs(s1r, s2l, contacts,  depth+1);
                        parallel_dfs_bvh_pairs(s1r, s2r, contacts, depth+1);
                    },
                );
            } else{
                parallel_dfs_bvh_pairs(s1l, s2l, contacts,  depth + 1);
                parallel_dfs_bvh_pairs(s1l, s2r, contacts,  depth+1);
                parallel_dfs_bvh_pairs(s1r, s2l, contacts,  depth+1);
                parallel_dfs_bvh_pairs(s1r, s2r, contacts, depth+1);
            }
        }
    }
}