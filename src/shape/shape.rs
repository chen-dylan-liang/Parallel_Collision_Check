use apollo_rust_lie::{EuclideanSpaceElement, LieGroupElement};
use apollo_rust_mesh_utils::trimesh::TriMesh;
use apollo_rust_spatial::lie::se3_implicit_quaternion::LieGroupISE3q;
use apollo_rust_spatial::vectors::{ApolloVector3Trait, V3};
use parry3d_f64::math::Point;
use parry3d_f64::transformation::convex_hull;

pub trait ShapeTrait {
    fn support(&self, dir: &V3, shape_pose: &LieGroupISE3q) -> V3;
    fn aabb(&self, shape_pose: &LieGroupISE3q) -> (V3, V3);
}

pub struct ConvexPolyhedron(pub TriMesh);

impl ConvexPolyhedron {
    pub fn new(input_mesh: &TriMesh)->Self{
        Self(input_mesh.to_convex_hull())
    }
    pub fn from_points(points: &[V3])->Self{
        let (ch_points, ch_indices) = convex_hull(points.iter().map(|x| Point::from_slice(x)).collect());
        let points: Vec<[f64; 3]> = ch_points.iter().map(|x| [x[0], x[1], x[2]] ).collect();
        let indices: Vec<[usize; 3]> = ch_indices.iter().map(|x| [ x[0] as usize, x[1] as usize, x[2] as usize]).collect();
        Self(TriMesh {
            points,
            indices
        })
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

    fn aabb(&self, shape_pose: &LieGroupISE3q) -> (V3, V3) {
       let mut min_v = V3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
       let mut max_v = V3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
       for point in self.0.points.iter(){
           let cur_point = shape_pose.0*V3::from_column_slice(point);
           min_v = min_v.inf(&cur_point);
           max_v = max_v.sup(&cur_point);
       }
        (min_v, max_v)
    }
}


