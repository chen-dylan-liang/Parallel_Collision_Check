use apollo_rust_spatial::vectors::V3;
use std::ops::{Add, Div, Neg, Sub};
use crate::shape::shape::ShapeTrait;
use apollo_rust_spatial::lie::se3_implicit_quaternion::LieGroupISE3q;
use rayon::prelude::*;

const _PROXIMITY_TOL: f64 =1e-6;
const _PROXIMITY_MAX_ITERS: usize = 100;

#[derive(Clone)]
struct ThreeSimplex{
    arr: [V3;4],
    len: usize,
}

#[derive(Clone)]
struct GJKFeature{
    pub v: V3, // closet point to origin on the simplex
    pub simplex: ThreeSimplex, // the simplex
    pub d: f64, // distance from the simplex to origin
}


impl GJKFeature{
    pub fn new(v: V3, arr: [V3; 4], len: usize)->Self{
        Self{v, simplex: ThreeSimplex::new_with_data(arr,len), d:v.norm()}
    }
    pub fn min(self, other: GJKFeature) -> GJKFeature{
        if self.d < other.d { self } else { other }
    }
}

impl ThreeSimplex {
    pub fn new() -> Self {
        Self{arr:[V3::zeros();4], len:0}
    }
    pub fn new_with_data(arr: [V3;4], len: usize) -> Self {
        Self{arr,len}
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn add(&mut self, point: V3) {
        self.arr[self.len] = point;
        self.len += 1;
    }

    pub fn find_and_reduce(&mut self) -> (V3, f64) {
        match self.len {
            1 => {  (self.arr[0], self.arr[0].norm()) },
            2 => {let f=closest_to_origin_on_line(&self.arr[0], &self.arr[1]); self.arr=f.simplex.arr;self.len=f.simplex.len; (f.v, f.d)},
            3 =>{ let f=closest_to_origin_on_triangle(&self.arr[0], &self.arr[1], &self.arr[2]); self.arr=f.simplex.arr;self.len=f.simplex.len; (f.v, f.d)},
            _ =>{ if tetrahedron_contains_origin(&self.arr[0],&self.arr[1], &self.arr[2], &self.arr[3]) {(V3::zeros(),0.0)} else {let f = closest_to_origin_on_tetrahedron(&self.arr[0],&self.arr[1], &self.arr[2], &self.arr[3]); self.arr=f.simplex.arr;self.len=f.simplex.len; (f.v, f.d)}},
        }
    }


}

fn closest_to_origin_on_line(a:&V3, b:&V3) -> GJKFeature {
    let ab = b.sub(a);
    let t = -(a.dot(&ab)).div(&(ab.dot(&ab)));
    if t<0.0{
        GJKFeature::new(a.clone(),[a.clone(),b.clone(),V3::zeros(), V3::zeros()],2)
    }
    else if t>1.0{
        GJKFeature::new(b.clone(),[a.clone(),b.clone(),V3::zeros(), V3::zeros()],2)
    } else{
        let closest = ab.scale(t).add(a);
        GJKFeature::new(closest,[a.clone(),b.clone(),V3::zeros(), V3::zeros()],2)
    }
}

fn closest_to_origin_on_triangle(a:&V3, b:&V3,c:&V3) -> GJKFeature {
    // compute the barycentric coordinates of the projection of the origin.
    let ab = b.sub(a);
    let ac = c.sub(a);
    let ao = a.neg();
    let d1 = ab.dot(&ao);
    let d2 = ac.dot(&ao);
    let d00 = ab.dot(&ab);
    let d01 = ab.dot(&ac);
    let d11 = ac.dot(&ac);
    let denom = d00 * d11 - d01 * d01;
    let u = (d11 * d1 - d01 * d2) / denom;
    let v = (d00 * d2 - d01 * d1) / denom;

    // inside
    if u > 0.0 && v > 0.0 && (u + v) < 1.0 {
        let closest = u*b+v*c+(1.0-u-v)*a;
        return GJKFeature::new(closest, [a.clone(), b.clone(), c.clone(), V3::zeros()],3);
    }

    // on edge
    let e1 = closest_to_origin_on_line(a, b);
    let e2 = closest_to_origin_on_line(b, c);
    let e3 = closest_to_origin_on_line(a, c);
    e3.min(e2).min(e1)
}



fn closest_to_origin_on_tetrahedron(a: &V3, b: &V3, c: &V3, d: &V3) -> GJKFeature {
    let f1 = closest_to_origin_on_triangle(a, b, c);
    let f2 = closest_to_origin_on_triangle(a, b, d);
    let f3 = closest_to_origin_on_triangle(a, c, d);
    let f4 = closest_to_origin_on_triangle(b, c, d);
    f4.min(f3).min(f2).min(f1)
}

fn signed_volume(a: &V3, b: &V3, c: &V3, d: &V3) -> f64 {
    (b - a).cross(&(c - a)).dot(&(d - a))
}

fn tetrahedron_contains_origin(a: &V3, b: &V3, c: &V3, d: &V3) -> bool {
    let vol = signed_volume(a, b, c, d);
    if vol.abs() < _PROXIMITY_TOL {
        return false;
    }

    // compute the signed volumes of tetrahedra that replace one vertex with the origin.
    let v1 = signed_volume(&V3::zeros(), b, c, d);
    let v2 = signed_volume(a,&V3::zeros(), c, d);
    let v3 = signed_volume(a, b, &V3::zeros(), d);
    let v4 = signed_volume(a, b, c, &V3::zeros());

    if vol > 0.0 {
        v1 > -_PROXIMITY_TOL && v2 > -_PROXIMITY_TOL && v3 > -_PROXIMITY_TOL && v4 > -_PROXIMITY_TOL
    } else {
        v1 < _PROXIMITY_TOL && v2 < _PROXIMITY_TOL && v3 < _PROXIMITY_TOL && v4 < _PROXIMITY_TOL
    }
}


pub fn gjk_contact<S1: ShapeTrait, S2: ShapeTrait>(shape1: &S1, pose1: &LieGroupISE3q, shape2: &S2, pose2:&LieGroupISE3q) -> (V3, f64) {
    let mut simplex = ThreeSimplex::new();
    let mut dir = pose1.0.translation.vector.sub(&pose2.0.translation.vector);
    if dir.norm_squared() > 1e-6 {dir=dir.normalize()} else {dir=V3::new(1.0, 0.0, 0.0)};
    let mut support = shape1.support(&dir, pose1).sub(shape2.support(&dir.neg(), pose2));
    simplex.add(support);
    let mut dist=support.norm();
    let mut iter=0;
    while iter < _PROXIMITY_MAX_ITERS {
        (dir, dist) = simplex.find_and_reduce();
        // intersected
        if dist < _PROXIMITY_TOL && simplex.len()==4 {return (V3::zeros(), 0.0);}
        dir = dir.normalize();
        support = shape1.support(&dir.neg(), pose1).sub(shape2.support(&dir, pose2));
        let proj = support.dot(&dir);
        //the simplex closet to the origin was found
        if dist < proj+_PROXIMITY_TOL {
            return (dir, dist);
        }
        // proceed to origin
        simplex.add(support);
        iter+=1;
    }
    (dir, dist)
}

#[derive(Debug)]
pub struct Contact {
    pub i: usize,
    pub j: usize,
    pub normal: V3,
    pub depth:  f64,
}

pub fn serial_narrow_phase(
    pairs:   &[(usize, usize)],
    shapes:  &[&dyn ShapeTrait],
    poses:   &[LieGroupISE3q],
)-> Vec<Contact>{
    assert_eq!(shapes.len(), poses.len(),
               "shapes and poses slices must have the same length");
    pairs.iter().filter_map(
        |&(i, j)
        | {
            let (p, d) = gjk_contact(
                shapes[i], &poses[i],
                shapes[j], &poses[j],
            );
            (d > 0.0).then(|| Contact { i, j, normal: p, depth: d })
        }).collect()
}

// embarrassingly parallelized narrow phase using Rayon parallel iterator
pub fn parallel_narrow_phase(
    pairs:   &[(usize, usize)],
    shapes:  &[&dyn ShapeTrait],
    poses:   &[LieGroupISE3q],
) -> Vec<Contact>
{
    assert_eq!(shapes.len(), poses.len(),
               "shapes and poses slices must have the same length");

    pairs.par_iter()
        .filter_map(
            |&(i, j)
            | {
            let (p, d) = gjk_contact(
                shapes[i], &poses[i],
                shapes[j], &poses[j],
            );
            (d > 0.0).then(|| Contact { i, j, normal: p, depth: d })
        }).collect()
}