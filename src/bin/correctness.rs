use apollo_rust_spatial::lie::se3_implicit_quaternion::LieGroupISE3q;
use apollo_rust_spatial::vectors::V3;
use parallel_collision_detection::{serial_parry_gjk, serial_double_phase_collision_check, parallel_double_phase_collision_check, generate_random_hulls, my_hulls_to_parry_hulls};
use parallel_collision_detection::gjk::gjk::{Contact, serial_narrow_phase_check, parallel_narrow_phase_check};
use parallel_collision_detection::shape::shape::ShapeTrait;

fn main() {
    let hulls =  generate_random_hulls(1000, (4, 10), (V3::new(-1.0, -1.0, -1.0), V3::new(1.0, 1.0, 1.0)));
    let parry_hulls = my_hulls_to_parry_hulls(&hulls);
    let poses: Vec<_> = (0..hulls.len()).map(|_| LieGroupISE3q::new_random()).collect();
    let mut indices = Vec::new();
    for i in 0..hulls.len() {
        for j in i + 1..hulls.len() {         // avoid (i,i) and dupes
            indices.push((i, j));
        }
    }
    let c1 = serial_parry_gjk(&indices, &parry_hulls, &poses);
    let c2 = serial_narrow_phase_check(&indices, &hulls, &poses);
    println!("{:?}, {:?}", c1.len(), c2.len());
    if c1.len()!=c2.len() {panic!("Size mismatch! c1_size={}, c2_size={}", c1.len(), c2.len());}
    for (p,(contact1, contact2)) in c1.iter().zip(c2.iter()).enumerate() {
        if !(((contact1.i==contact2.i) && (contact1.j==contact2.j)) || ((contact1.i==contact2.j) && (contact1.j==contact2.i))) {
            panic!("{}-th pair mismatch! c1:{},{}; c2:{},{}", p, contact1.i, contact1.j, contact2.i, contact2.j);
        }
    }
}