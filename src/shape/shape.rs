use apollo_rust_lie::{EuclideanSpaceElement, LieGroupElement};
use apollo_rust_mesh_utils::trimesh::TriMesh;
use apollo_rust_spatial::lie::se3_implicit_quaternion::LieGroupISE3q;
use apollo_rust_spatial::vectors::{ApolloVector3Trait, V3};
use std::ops::{Add, Div, Neg, Sub};

pub trait ShapeTrait {
    fn support(&self, dir: &V3, shape_pose: &LieGroupISE3q) -> V3;
    fn aabb(&self, shape_pose: &LieGroupISE3q) -> (V3, V3);
}

pub struct ConvexPolyhedron(pub TriMesh);

impl ConvexPolyhedron {
    pub fn new(input_mesh: &TriMesh)->Self{
        Self(input_mesh.to_convex_hull())
    }
}

pub struct Sphere{
    // suppose the center is always the origin
    pub radius: f64,
}

impl Sphere {
    pub fn new(radius:f64)->Self{
        Self{radius}
    }
}

pub struct Cuboid{
    pub half_extents: V3
}

impl Cuboid{
    pub fn new(x:f64,y:f64,z:f64)->Self{
        Self{half_extents: V3::new(x,y,z)}
    }
}

impl ShapeTrait for ConvexPolyhedron {
    fn support(&self, dir: &V3, shape_pose: &LieGroupISE3q) -> V3 {
        let local_dir = shape_pose.0.rotation.inverse() * dir;
        let mut max_point = V3::from_column_slice(&self.0.points[0]);
        let mut max_proj = max_point.dot(&local_dir);
        for point in self.0.points.iter().skip(1) {
            let cur_point = V3::from_column_slice(point);
            let proj =cur_point.dot(&local_dir);
            if proj > max_proj {
                max_proj = proj;
                max_point = cur_point;
            }
        }
        shape_pose.0.rotation*max_point+shape_pose.0.translation.vector
    }
}

impl ShapeTrait for Sphere {
    fn support(&self, dir: &V3, shape_pose: &LieGroupISE3q) -> V3 {
        shape_pose.0.translation.vector+self.radius*dir.normalize()
    }
}

impl ShapeTrait for Cuboid {
    fn support(&self, dir: &V3, shape_pose: &LieGroupISE3q) -> V3 {
        let local_dir = shape_pose.0.rotation.inverse() * dir;
        let mut max_point = V3::new(0.0,0.0,0.0);
        let mut max_proj = f64::NEG_INFINITY;
        for i in 0..8{
            let mut point=self.half_extents.clone();
            let mut _i=i.clone();
            for j in 0..3{
                if (_i%2)==1 {point[j]=point[j].neg()}
                _i=_i.div(2);
            }
            let proj =point.dot(&local_dir);
            if proj > max_proj {
                max_proj = proj;
                max_point = point;
            }
        }
        shape_pose.0.rotation*max_point+shape_pose.0.translation.vector
    }
}

