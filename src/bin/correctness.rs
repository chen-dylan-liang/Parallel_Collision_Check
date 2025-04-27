use apollo_rust_spatial::lie::se3_implicit_quaternion::LieGroupISE3q;
use apollo_rust_spatial::vectors::V3;
use parallel_collision_detection::{serial_parry_gjk, serial_double_phase_collision_check, parallel_double_phase_collision_check, generate_random_hulls, my_hulls_to_parry_hulls, parallel_parry_gjk};
use parallel_collision_detection::gjk::gjk::{Contact, serial_narrow_phase_check, parallel_narrow_phase_check};
use parallel_collision_detection::shape::shape::ShapeTrait;

fn check(ground_truth: &[Contact], res: &[Contact], name: &str){
    println!("Checking {}",name);
    println!("{}",res.len());
    let mut c1 = ground_truth.iter().map(|contact| (contact.i, contact.j)).collect::<Vec<_>>();
    let mut c2 = res.iter().map(|contact| (contact.i, contact.j)).collect::<Vec<_>>();
    c1.sort();
    c2.sort();
    if ground_truth.len()!=res.len() {panic!("Size mismatch! ground_truth_size={}, res={}", ground_truth.len(), res.len());}
    for (p,(contact1, contact2)) in c1.iter().zip(c2.iter()).enumerate() {
        // i < j always holds
        if !((contact1.0==contact2.0) && (contact1.1==contact2.1)) {
            panic!("{}-th pair mismatch! ground_truth:{},{}; res:{},{}", p, contact1.0, contact1.1, contact2.0, contact2.1);
        }
    }
    println!("{} passed",name);
}
fn main() {
    let mut hulls = generate_random_hulls(100, (50, 100), (V3::new(-1.0, -1.0, -1.0), V3::new(0.0, 0.0, 0.0)));
    let mut hull2 = generate_random_hulls(100, (50, 100), (V3::new(0.0, 0.0, 0.0), V3::new(1.0, 1.0, 1.0)));
    hulls.append(&mut hull2);
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
    let c3 = parallel_parry_gjk(&indices, &parry_hulls, &poses);
    let c4 = parallel_narrow_phase_check(&indices, &hulls, &poses);
    let c5 = serial_double_phase_collision_check(&hulls, &poses,1);
    let c6 = parallel_double_phase_collision_check(&hulls, &poses,1);
    check(&c1, &c2, "my serial narrow");
    check(&c1, &c3, "parry's parallel narrow");
    check(&c1, &c4, "my parallel narrow");
    check(&c1, &c5, "my serial double");
    check(&c1, &c6, "my parallel double");
}