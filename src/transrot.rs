//! Translation and rotation

use crate::{
    aabb::Aabb,
    objects::{HitRecord, Hittable, HittableObj},
    Point, Ray, Vec3,
};

/// Translate an object
pub struct Translate {
    obj: HittableObj,
    offset: Vec3,
}
impl Translate {
    pub fn new(obj: HittableObj, offset: Vec3) -> Self {
        Self { obj, offset }
    }
}
impl Hittable for Translate {
    fn try_hit(
        &self,
        ray: &crate::Ray,
        t_min: f64,
        t_max: f64,
    ) -> Option<crate::objects::HitRecord> {
        let moved_ray = Ray::new(ray.orig - self.offset, ray.dir, ray.time);

        self.obj.try_hit(&moved_ray, t_min, t_max).map(|mut hr| {
            hr.p += self.offset;
            hr
        })
    }

    fn try_bounding_box(&self, time0: f64, time1: f64) -> Option<Aabb> {
        self.obj
            .try_bounding_box(time0, time1)
            .map(|output_box| Aabb::new(output_box.min + self.offset, output_box.max + self.offset))
    }
}

// TODO(mkagie) Revisit this and make it more generic
pub struct RotateY {
    obj: HittableObj,
    sin_theta: f64,
    cos_theta: f64,
    bbox: Option<Aabb>,
}
impl RotateY {
    pub fn new(obj: HittableObj, angle_deg: f64) -> Self {
        let radians = angle_deg.to_radians();
        let sin_theta = radians.sin();
        let cos_theta = radians.cos();
        // TODO(mkagie) A little unclear what to do here if there is no bounding box
        let bbox = obj.try_bounding_box(0.0, 1.0).unwrap();

        let mut min = Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut max = Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);

        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let x = i as f64 * bbox.max[0] + (1 - i) as f64 * bbox.min[0];
                    let y = j as f64 * bbox.max[1] + (1 - j) as f64 * bbox.min[1];
                    let z = k as f64 * bbox.max[2] + (1 - k) as f64 * bbox.min[2];

                    let new_x = cos_theta * x + sin_theta * z;
                    let new_z = -sin_theta * x + cos_theta * z;

                    let tester = Vec3::new(new_x, y, new_z);

                    for c in 0..3 {
                        min[c] = min[c].min(tester[c]);
                        max[c] = max[c].max(tester[c]);
                    }
                }
            }
        }
        let bbox = Aabb::new(min, max);
        Self {
            obj,
            sin_theta,
            cos_theta,
            bbox: Some(bbox),
        }
    }
}
impl Hittable for RotateY {
    fn try_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<crate::objects::HitRecord> {
        let mut origin = ray.orig;
        let mut direction = ray.dir;

        origin[0] = self.cos_theta * ray.orig[0] - self.sin_theta * ray.orig[2];
        origin[2] = self.sin_theta * ray.orig[0] + self.cos_theta * ray.orig[2];

        direction[0] = self.cos_theta * ray.dir[0] - self.sin_theta * ray.dir[2];
        direction[2] = self.sin_theta * ray.dir[0] + self.cos_theta * ray.dir[2];

        let rotated_r = Ray::new(origin, direction, ray.time);

        self.obj.try_hit(&rotated_r, t_min, t_max).map(|rec| {
            let mut p = rec.p;
            let mut normal = rec.normal;

            p[0] = self.cos_theta * rec.p[0] + self.sin_theta * rec.p[2];
            p[2] = -self.sin_theta * rec.p[0] + self.cos_theta * rec.p[2];

            normal[0] = self.cos_theta * rec.normal[0] + self.sin_theta * rec.normal[2];
            normal[2] = -self.sin_theta * rec.normal[0] + self.cos_theta * rec.normal[2];

            HitRecord::new(p, rec.t, &rotated_r, &normal, rec.material, rec.u, rec.v)
        })
    }

    fn try_bounding_box(&self, _time0: f64, _time1: f64) -> Option<Aabb> {
        self.bbox.clone()
    }
}
