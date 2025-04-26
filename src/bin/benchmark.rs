use std::rt::panic_count::increase;
use apollo_rust_spatial::lie::se3_implicit_quaternion::LieGroupISE3q;
use apollo_rust_spatial::vectors::V3;
use rand::Rng;
use parallel_collision_detection::{shape, bvh, gjk};
use parallel_collision_detection::gjk::gjk::{parallel_narrow_phase_check, serial_narrow_phase_check};
use shape::shape::ConvexPolyhedron as ConvexHull;
use parallel_collision_detection::parallel_double_phase_collision_check;
use parallel_collision_detection::serial_double_phase_collision_check;
use parallel_collision_detection::shape::shape::ShapeTrait;
use parallel_collision_detection::generate_random_hulls;

fn main() {
    // random shapes
    let hull_vec: Vec<Box<dyn ShapeTrait>> =
        generate_random_hulls(
            100,
            (4, 10),
            (V3::new(-1.0, -1.0, -1.0), V3::new(1.0, 1.0, 1.0)),
        )
            .into_iter()
            .map(|h| Box::new(h) as Box<dyn ShapeTrait>)
            .collect();

    let hulls: Vec<&dyn ShapeTrait> =
        hull_vec.iter().map(Box::as_ref).collect();
    // random poses
    let poses=(0..hulls.len())
    .map(|_| LieGroupISE3q::new_random()).collect();
    // 1. serial narrow phase only
    let mut indices = Vec::new();
    for i in 0..hulls.len() {
        for j in i + 1..hulls.len() {         // avoid (i,i) and dupes
            indices.push((i, j));
        }
    }
    let c1 = serial_narrow_phase_check(&indices, &hulls, &poses);
    // 2. parallel narrow phase only
    let c2 = parallel_narrow_phase_check(&indices, &hulls, &poses);
    // 3. serial double phase
    let c3 = serial_double_phase_collision_check(&hulls, &poses, 8);
    // 4. parallel double phase
    let c4=parallel_double_phase_collision_check(&hulls, &poses, 8);
}