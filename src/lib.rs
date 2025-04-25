#![feature(iter_partition_in_place)]

use std::sync::Mutex;
use apollo_rust_spatial::lie::se3_implicit_quaternion::LieGroupISE3q;
use crate::bvh::par_bvh::{parallel_broad_phase_check, parallel_build_bvh};
use crate::bvh::srl_bvh::{serial_broad_phase_check, serial_build_bvh};
use crate::bvh::structs::AABB;
use crate::gjk::gjk::*;
use crate::shape::shape::ShapeTrait;

pub mod shape;
pub mod gjk;
pub mod bvh;

pub fn parallel_double_phase_collision_check(shapes: &[&dyn ShapeTrait],
                                             poses: &[LieGroupISE3q],
                                             cut_off: usize)->Vec<Contact>{
    // construct aabbs
    let aabbs:Vec<AABB>  = shapes.iter()
        .zip(poses.iter()).
        map(|(&shape, pose)|{ let (min,max)=shape.aabb(pose);
    AABB::new(min,max)}).collect();
    // broad phase
    let mut indices: Vec<usize> = (0..aabbs.len()).collect();
    let bvh = parallel_build_bvh(&mut indices, &aabbs, cut_off);
    let potential_pairs: Mutex<Vec<(usize, usize)>> = Mutex::new(Vec::new());
    parallel_broad_phase_check(&bvh, &bvh, &potential_pairs);
    // narrow phase
    let pairs=potential_pairs.into_inner().unwrap();
    parallel_narrow_phase_check(&pairs, shapes, poses)
}

pub fn serial_double_phase_collision_check(shapes: &[&dyn ShapeTrait],
                                           poses: &[LieGroupISE3q],
                                           cut_off: usize)->Vec<Contact>{
    // construct aabbs
    let aabbs:Vec<AABB>  = shapes.iter()
        .zip(poses.iter()).
        map(|(&shape, pose)|{ let (min,max)=shape.aabb(pose);
            AABB::new(min,max)}).collect();
    // broad phase
    let mut indices: Vec<usize> = (0..aabbs.len()).collect();
    let bvh = serial_build_bvh(&mut indices, &aabbs, cut_off);
    let mut potential_pairs = Vec::<(usize, usize)>::new();
    serial_broad_phase_check(&bvh, &bvh, &mut potential_pairs);
    // narrow phase
    serial_narrow_phase_check(&potential_pairs, shapes, poses)
}


