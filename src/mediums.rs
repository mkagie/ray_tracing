//! Mediums

use rand::Rng;

use crate::{
    materials::Isotropic,
    objects::{HitRecord, Hittable},
    prelude::HittableObj,
    textures::{SolidColor, Texture},
    Color, Material, Vec3,
};

/// Constant Medium
pub struct ConstantMedium {
    boundary: HittableObj,
    phase_function: Material,
    neg_inv_density: f64,
}
impl ConstantMedium {
    pub fn new(b: HittableObj, d: f64, a: Texture) -> Self {
        Self {
            boundary: b,
            phase_function: Box::new(Isotropic::new(a)),
            neg_inv_density: -1.0 / d,
        }
    }

    pub fn from_color(b: HittableObj, d: f64, c: Color) -> Self {
        let t = Box::new(SolidColor::new(c));
        Self::new(b, d, t)
    }
}
impl Hittable for ConstantMedium {
    fn try_hit(
        &self,
        ray: &crate::Ray,
        t_min: f64,
        t_max: f64,
    ) -> Option<crate::objects::HitRecord> {
        let mut rng = rand::thread_rng();
        if let Some(mut rec1) = self.boundary.try_hit(ray, f64::NEG_INFINITY, f64::INFINITY) {
            if let Some(mut rec2) = self.boundary.try_hit(ray, rec1.t + 0.0001, f64::INFINITY) {
                if rec1.t < t_min {
                    rec1.t = t_min;
                }
                if rec2.t > t_max {
                    rec2.t = t_max;
                }
                if rec1.t >= rec2.t {
                    return None;
                }

                if rec1.t < 0.0 {
                    rec1.t = 0.0;
                }

                let ray_length = ray.dir.norm();
                let distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
                let hit_distance = self.neg_inv_density * rng.gen::<f64>().log(10.0);

                if hit_distance > distance_inside_boundary {
                    return None;
                }
                let t = rec1.t + hit_distance / ray_length;
                let p = ray.get(t);
                let normal = Vec3::new(1.0, 0.0, 0.0);
                Some(HitRecord::new(
                    p,
                    t,
                    ray,
                    &normal,
                    dyn_clone::clone_box(&*self.phase_function),
                    0.0,
                    0.0,
                ))
            } else {
                None
            }
        } else {
            None
        }
    }

    fn try_bounding_box(&self, time0: f64, time1: f64) -> Option<crate::bvh::Aabb> {
        self.boundary.try_bounding_box(time0, time1)
    }
}
