//! Objects
use crate::{Material, Point, Ray, Vec3};

pub trait Hittable {
    fn try_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}

#[derive(Default)]
pub struct HittableList(Vec<Box<dyn Hittable + Send + Sync>>);
impl HittableList {
    pub fn add(&mut self, boxed_obj: Box<dyn Hittable + Send + Sync>) {
        self.0.push(boxed_obj)
    }

    pub fn clear(&mut self) {
        self.0.clear()
    }
}
impl Hittable for HittableList {
    // TODO(mkagie) Refactor this more rusty
    fn try_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut closest_so_far = t_max;
        let mut hr_final = None;

        for obj in &self.0 {
            if let Some(hr) = obj.try_hit(ray, t_min, closest_so_far) {
                closest_so_far = hr.t;
                hr_final = Some(hr)
            }
        }
        hr_final
    }
}

/// Represents a hit
pub struct HitRecord {
    /// Point of intersection
    pub p: Point,
    /// Normali vector
    pub normal: Vec3,
    /// roots of intersection
    pub t: f64,
    /// Wither we are facing the normal
    pub front_face: bool,
    /// Material
    pub material: Material,
}
impl HitRecord {
    pub fn new(p: Point, t: f64, ray: &Ray, outward_normal: &Vec3, material: Material) -> Self {
        let front_face = ray.dir.dot(outward_normal) < 0.0;
        let mut normal = outward_normal.to_owned();
        if !front_face {
            normal = -normal;
        }
        Self {
            p,
            normal,
            t,
            front_face,
            material,
        }
    }
}

pub struct Sphere {
    pub center: Point,
    pub radius: f64,
    pub material: Material,
}
impl Sphere {
    pub fn new(center: Point, radius: f64, material: Material) -> Self {
        Self {
            center,
            radius,
            material,
        }
    }
}
impl Hittable for Sphere {
    fn try_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = ray.orig - self.center;
        let a = ray.dir.norm().powi(2);
        let half_b = oc.dot(&ray.dir);
        let c = oc.norm().powi(2) - self.radius.powi(2);
        let discriminant = half_b.powi(2) - a * c;
        if discriminant < 0.0 {
            return None;
        }

        // Find the nearest root that lies in the acceptable range
        let sqrtd = discriminant.sqrt();
        let mut root = (-half_b - sqrtd) / a;
        if root < t_min || t_max < root {
            root = (-half_b + sqrtd) / a;
            if root < t_min || t_max < root {
                return None;
            }
        }
        let p = ray.get(root);
        let t = root;
        let outward_normal = ((p - self.center) / self.radius).normalize();
        // NOTE -- we use dyn_clone here because self.material is a trait object -- you cannot
        // clone a trait object
        Some(HitRecord::new(
            p,
            t,
            ray,
            &outward_normal,
            dyn_clone::clone_box(&*self.material),
        ))
    }
}
