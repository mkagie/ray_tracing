//! Objects
use std::cmp::Ordering;

use crate::{
    aabb::Aabb,
    materials::{self, MaterialConfig},
    utils::SerdeVector,
    Material, Point, Ray, Vec3,
};
use serde::{Deserialize, Serialize};

pub type HittableObj = Box<dyn Hittable + Send + Sync>;

pub trait Hittable {
    fn try_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;

    fn try_bounding_box(&self, time0: f64, time1: f64) -> Option<Aabb>;
}

#[derive(Default)]
pub struct HittableList(pub Vec<HittableObj>);
impl HittableList {
    pub fn add(&mut self, boxed_obj: HittableObj) {
        self.0.push(boxed_obj)
    }

    pub fn clear(&mut self) {
        self.0.clear()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn sort_by<F>(&mut self, f: F)
    where
        F: FnMut(&HittableObj, &HittableObj) -> Ordering,
    {
        self.0.sort_by(f);
    }

    pub fn from_config(config: HittableListConfig) -> Self {
        let mut s = Self::default();
        for obj_cfg in config.objects {
            let obj = Sphere::from_config(obj_cfg);
            s.add(Box::new(obj));
        }
        s
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

    fn try_bounding_box(&self, time0: f64, time1: f64) -> Option<Aabb> {
        if self.0.is_empty() {
            return None;
        }
        let mut output_box: Option<Aabb> = None;

        for obj in &self.0 {
            match obj.try_bounding_box(time0, time1) {
                Some(tmp_box) => {
                    output_box = if let Some(output_box) = output_box {
                        Some(output_box.surrounding_box(&tmp_box))
                    } else {
                        Some(tmp_box)
                    };
                }
                None => return None,
            }
        }
        output_box
    }
}

// TODO(mkagie) extend out of spheres once we have more
/// Hittable List Config
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct HittableListConfig {
    objects: Vec<SphereConfig>,
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

    pub fn from_config(config: SphereConfig) -> Self {
        Self::new(
            config.center.into(),
            config.radius,
            materials::Generator::from_config(config.material),
        )
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

    fn try_bounding_box(&self, _time0: f64, _time1: f64) -> Option<Aabb> {
        let v = Vec3::new(self.radius, self.radius, self.radius);
        Some(Aabb::new(self.center - v, self.center + v))
    }
}

/// Sphere config
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SphereConfig {
    pub center: SerdeVector,
    pub radius: f64,
    pub material: MaterialConfig,
}

// TODO(mkagie) Merge this with sphere
/// Moving Sphere
pub struct MovingSphere {
    initial_center: Point,
    final_center: Point,
    initial_time: f64,
    final_time: f64,
    radius: f64,
    material: Material,
}
impl MovingSphere {
    pub fn new(
        initial_center: Point,
        final_center: Point,
        initial_time: f64,
        final_time: f64,
        radius: f64,
        material: Material,
    ) -> Self {
        Self {
            initial_center,
            final_center,
            initial_time,
            final_time,
            radius,
            material,
        }
    }

    pub fn center(&self, time: f64) -> Point {
        self.initial_center
            + ((time - self.initial_time) / (self.final_time - self.initial_time))
                * (self.final_center - self.initial_center)
    }
}
impl Hittable for MovingSphere {
    fn try_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = ray.orig - self.center(ray.time);
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
        let outward_normal = ((p - self.center(ray.time)) / self.radius).normalize();
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

    fn try_bounding_box(&self, time0: f64, time1: f64) -> Option<Aabb> {
        let v = Vec3::new(self.radius, self.radius, self.radius);

        let box0 = Aabb::new(self.center(time0) - v, self.center(time0) + v);
        let box1 = Aabb::new(self.center(time1) - v, self.center(time1) + v);

        Some(box0.surrounding_box(&box1))
    }
}
