#![feature(iter_partition_in_place)]

use std::collections::HashSet;
use std::hash::Hash;
use std::sync::Mutex;
use std::time::Instant;
use apollo_rust_spatial::lie::se3_implicit_quaternion::LieGroupISE3q;
use apollo_rust_spatial::vectors::V3;
use crate::bvh::par_bvh::{parallel_broad_phase_check, parallel_build_bvh};
use crate::bvh::srl_bvh::{serial_broad_phase_check, serial_build_bvh};
use crate::bvh::structs::AABB;
use crate::gjk::gjk::*;
use crate::shape::shape::ShapeTrait;
use parry3d_f64::shape::{ConvexPolyhedron as ParryConvexHull, TriMesh};
use parry3d_f64::query::{distance as parry_distance, DefaultQueryDispatcher};
use crate::shape::shape::ConvexPolyhedron as ConvexHull;
use parry3d_f64::math::Point as ParryPoint;
use rand::Rng;
use rayon::prelude::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;

pub mod shape;
pub mod gjk;
pub mod bvh;

pub fn generate_random_hulls(n: usize, vn_range: (usize, usize), point_range: (V3, V3)) -> Vec<ConvexHull> {
    let mut rng = rand::thread_rng();
    let mut points: Vec<V3> = Vec::new();
    let mut ret: Vec<ConvexHull> = Vec::new();
    for _ in 0..n {
        let vn: usize = rng.gen_range(vn_range.0..vn_range.1);
        points.clear();
        for _ in 0..vn{
            points.push(V3::new(rng.gen_range(point_range.0.x..point_range.1.x),
            rng.gen_range(point_range.0.y..point_range.1.y),
                                rng.gen_range(point_range.0.z..point_range.1.z)));
        }
        ret.push(ConvexHull::from_points(&points));
    }
    ret
}

pub fn my_hulls_to_parry_hulls(hulls: &[ConvexHull])->Vec<ParryConvexHull>{
    let mut parry_hulls: Vec<ParryConvexHull> = Vec::new();
    for h in hulls.iter(){
        let pts = h.0.points.iter().map(|p| ParryPoint::from_slice(p)).collect();
        let faces: Vec<_>  = h.0.indices.iter().map(|f| [f[0] as u32, f[1] as u32, f[2] as u32]).collect();
        parry_hulls.push(ParryConvexHull::from_convex_mesh(pts, &faces).unwrap());
    }
    parry_hulls
}

pub fn parallel_double_phase_collision_check(shapes: &[ConvexHull],
                                             poses: &[LieGroupISE3q],
                                             cut_off: usize)->Vec<Contact>{
    // construct aabbs
    //let t=Instant::now();
    let aabbs:Vec<AABB>  = shapes.par_iter()
        .zip(poses.par_iter()).
        map(|(shape, pose)|{ let (min,max)=shape.aabb(pose);
    AABB::new(min,max)}).collect();
    //println!("para build AABBs {:?}",t.elapsed());
    
    // broad phase
    let mut indices: Vec<usize> = (0..aabbs.len()).collect();
    //let t=Instant::now();
    let bvh = parallel_build_bvh(&mut indices, &aabbs, cut_off);
   // println!("para build BVH {:?}", t.elapsed());

    //let potential_pairs: Mutex<Vec<(usize, usize)>> = Mutex::new(Vec::new());
    //let t = Instant::now();
    let pairs=parallel_broad_phase_check(&*bvh, &*bvh);
    //println!("para broad check {:?}", t.elapsed());
    // narrow phase
    //let pairs:Vec<_>=potential_pairs.into_inner().unwrap();
    //let t=Instant::now();
   let ret=parallel_narrow_phase_check(&pairs, shapes, poses);
   // println!("para narrow check {:?}", t.elapsed());
    ret

}

pub fn serial_double_phase_collision_check(shapes: &[ConvexHull],
                                           poses: &[LieGroupISE3q],
                                           cut_off: usize)->Vec<Contact>{
    // construct aabbs
    //let t= Instant::now();
    let aabbs:Vec<AABB>  = shapes.iter()
        .zip(poses.iter()).
        map(|(shape, pose)|{ let (min,max)=shape.aabb(pose);
            AABB::new(min,max)}).collect();
   // println!("serial build AABBs {:?}",t.elapsed());
    
    // broad phase
    let mut indices: Vec<usize> = (0..aabbs.len()).collect();
    //let t=Instant::now();
    let bvh = serial_build_bvh(&mut indices, &aabbs, cut_off);
   // println!("serial build BVH {:?}", t.elapsed());

    
    //let t = Instant::now();
    let pairs = serial_broad_phase_check(&*bvh, &*bvh);
    //println!("serial broad check {:?}", t.elapsed());
    // narrow phase
   // let pairs:Vec<_> = potential_pairs.into_iter().collect();
    //let t = Instant::now();
   let ret= serial_narrow_phase_check(&pairs, shapes, poses);
    //println!("serial narrow check {:?}", t.elapsed());
    ret

}

pub fn serial_parry_gjk(pairs: &[(usize, usize)], hulls: &[ParryConvexHull], poses: &[LieGroupISE3q])->Vec<Contact>{
        pairs.iter().filter_map(
            |&(i, j)
            | {
                let dist = parry_distance(&poses[i].0, &hulls[i], &poses[j].0, &hulls[j]).unwrap();
                (dist == 0.0).then(|| Contact { i, j})
            }).collect()
}

pub fn parallel_parry_gjk(pairs: &[(usize, usize)], hulls: &[ParryConvexHull], poses: &[LieGroupISE3q])->Vec<Contact>{
    pairs.par_iter().filter_map(
        |&(i, j)
        | {
            let dist = parry_distance(&poses[i].0, &hulls[i], &poses[j].0, &hulls[j]).unwrap();
            (dist ==0.0).then(|| Contact { i, j})
        }).collect()
}








