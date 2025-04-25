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

fn generate_random_hulls(n: usize, vn_range: (usize, usize), point_range: (V3, V3)) -> Vec<ConvexHull> {
    let mut rng = rand::thread_rng();
    let mut points: Vec<V3> = Vec::new();
    let mut ret: Vec<ConvexHull> = Vec::new();
    for _ in 0..n {
        let vn: usize = rng.gen_range(vn_range[0]..vn_range[1]);
        points.clear();
        for _ in 0..vn{
            points.push(rng.gen_range(point_range[0] ..point_range[1]));
        }
       ret.push(ConvexHull::from_points(&points));
    }
   ret
}

fn main() {
    // random shapes
    let hulls = generate_random_hulls(100,
                                      (4,10),
                                      (V3::new(-1.0,-1.0,-1.0), V3::new(1.0,1.0,1.0)))
        .into_iter()
        .map(|hull| Box::new(hull) as Box<dyn ShapeTrait>)
        .collect() // Vec<Box<dyn ShapeTrait>>
        .iter()
        .map(|h| h)        // h has type &Box<dyn ShapeTrait>
        .collect::<Vec<&dyn ShapeTrait>>();
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