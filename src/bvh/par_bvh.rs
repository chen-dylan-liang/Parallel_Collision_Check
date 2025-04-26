use std::sync::Mutex;
use apollo_rust_spatial::vectors::V3;
use rayon::prelude::*;
use super::structs::{AABB, BVHNode, BVHInternalNode, BVHLeafNode};
use rayon::slice::ParallelSlice;
use rayon::slice::ParallelSliceMut;




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

pub fn parallel_build_bvh(aabb_indices: &mut [usize], all_aabbs:&[AABB], cut_off_size:usize)->Box<dyn BVHNode>{
    // S=O(1)
    if aabb_indices.len() <= cut_off_size{
        return Box::new(BVHLeafNode::new( aabb_indices.to_vec(), all_aabbs));
    }
    // S=O(logn)
    let (axis,midpoint) = parallel_longest_extent_axis(aabb_indices, all_aabbs);
    // S=O(logn)
    let (indices_left, indices_right) = parallel_split_at_axis(aabb_indices, all_aabbs, axis, midpoint);
    // safeguard for degenerate cases where all bounding boxes equal
    if indices_left.is_empty() || indices_right.is_empty() {
        return Box::new(BVHLeafNode::new(aabb_indices.to_vec(), all_aabbs));
    }
    // do recursive calls
    let (left_tree, right_tree)=rayon::join(
        || parallel_build_bvh(indices_left, all_aabbs, cut_off_size),
        || parallel_build_bvh(indices_right, all_aabbs, cut_off_size));
    // S=O(1)
    let aabb = left_tree.union_aabb(&*right_tree);
    Box::new(BVHInternalNode::new(aabb, left_tree, right_tree))
}

pub fn parallel_broad_phase_check<'a>(
    s1: &'a dyn BVHNode,
    s2: &'a dyn BVHNode,
    contacts: &'a Mutex<Vec<(usize, usize)>>,
) {
    rayon::scope(|scope| {
        // Call the recursive function directly inside the scope
        parallel_dfs_bvh_pairs(s1, s2, contacts, scope);
    });
}

fn parallel_dfs_bvh_pairs<'scope>(
    s1: &'scope (dyn BVHNode + 'scope),
    s2: &'scope (dyn BVHNode + 'scope),
    contacts: &'scope Mutex<Vec<(usize, usize)>>,
    scope: &rayon::Scope<'scope>,
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
                    if i < j { buf.push((i, j)); }       // skip self/dup
                }
            }
        }

        (true, false) => {
            let (left, right) = s2.children();
            parallel_dfs_bvh_pairs(s1, left.unwrap(), contacts, scope);
            parallel_dfs_bvh_pairs(s1, right.unwrap(), contacts, scope);
        }
        (false, true) => {
            let (left, right) = s1.children();
            parallel_dfs_bvh_pairs(left.unwrap(), s2, contacts, scope);
            parallel_dfs_bvh_pairs(right.unwrap(), s2, contacts, scope);
        }

        (false, false) => {
            let (s1_left, s1_right) = s1.children();
            let (s2_left, s2_right) = s2.children();

            // Create local variables to avoid multiple unwraps
            let s1l = s1_left.unwrap();
            let s1r = s1_right.unwrap();
            let s2l = s2_left.unwrap();
            let s2r = s2_right.unwrap();

            // Spawn one pair in a new task
            let s1l_ref = s1l;
            let s2l_ref = s2l;
            scope.spawn(move |subscope| {
                parallel_dfs_bvh_pairs(s1l_ref, s2l_ref, contacts, subscope);
            });

            // Spawn another pair
            let s1l_ref = s1l;
            let s2r_ref = s2r;
            scope.spawn(move |subscope| {
                parallel_dfs_bvh_pairs(s1l_ref, s2r_ref, contacts, subscope);
            });

            // Spawn a third pair
            let s1r_ref = s1r;
            let s2l_ref = s2l;
            scope.spawn(move |subscope| {
                parallel_dfs_bvh_pairs(s1r_ref, s2l_ref, contacts, subscope);
            });

            // Current thread handles the fourth pair
            parallel_dfs_bvh_pairs(s1r, s2r, contacts, scope);
        }
    }
}