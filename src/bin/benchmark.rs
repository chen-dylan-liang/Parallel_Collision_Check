use apollo_rust_spatial::lie::se3_implicit_quaternion::LieGroupISE3q;
use apollo_rust_spatial::vectors::V3;
use rand::Rng;
use parallel_collision_detection::{shape, bvh, gjk, my_hulls_to_parry_hulls, serial_parry_gjk, parallel_parry_gjk};
use parallel_collision_detection::gjk::gjk::{parallel_narrow_phase_check, serial_narrow_phase_check};
use shape::shape::ConvexPolyhedron as ConvexHull;
use parallel_collision_detection::parallel_double_phase_collision_check;
use parallel_collision_detection::serial_double_phase_collision_check;
use parallel_collision_detection::shape::shape::ShapeTrait;
use parallel_collision_detection::generate_random_hulls;
use std::time::Instant;
use std::fs::File;
use std::io::Write;

fn main() ->std::io::Result<()> {
    for cut_off  in [1,2,4,8,16,32,64,128,256,512,1024] {
        let file_name = "./res/benchmark_".to_owned() + &*cut_off.to_string() +".txt";
        let mut file = File::create(file_name)?;
        for n in [10, 20, 50, 100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000] {
            writeln!(file, "n={}", n)?;
            let hulls = generate_random_hulls(n, (4, 20), (V3::new(-1.0, -1.0, -1.0), V3::new(1.0, 1.0, 1.0)));
            let parry_hulls = my_hulls_to_parry_hulls(&hulls);
            let poses: Vec<_> = (0..hulls.len()).map(|_| LieGroupISE3q::new_random()).collect();
            let mut indices = Vec::new();
            for i in 0..hulls.len() {
                for j in i + 1..hulls.len() {         // avoid (i,i) and dupes
                    indices.push((i, j));
                }
            }
            let t = Instant::now();
            serial_parry_gjk(&indices, &parry_hulls, &poses);
            writeln!(file, "serial_parry_narrow={:?}", t.elapsed())?;

            let t = Instant::now();
            serial_narrow_phase_check(&indices, &hulls, &poses);
            writeln!(file, "serial_our_narrow={:?}", t.elapsed())?;

            let t = Instant::now();
            parallel_parry_gjk(&indices, &parry_hulls, &poses);
            writeln!(file, "parallel_parry_narrow={:?}", t.elapsed())?;

            let t = Instant::now();
            parallel_narrow_phase_check(&indices, &hulls, &poses);
            writeln!(file, "parallel_our_narrow={:?}", t.elapsed())?;

            let t = Instant::now();
            serial_double_phase_collision_check(&hulls, &poses, 4);
            writeln!(file, "serial_double={:?}", t.elapsed())?;

            let t = Instant::now();
            parallel_double_phase_collision_check(&hulls, &poses,4);
            writeln!(file, "parallel_double={:?}", t.elapsed())?;
        }
    }
    Ok(())
}