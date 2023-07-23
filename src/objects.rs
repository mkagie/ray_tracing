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
    /// U,V surface coordinates
    pub u: f64,
    /// U,V surface coordinates
    pub v: f64,
}
impl HitRecord {
    pub fn new(
        p: Point,
        t: f64,
        ray: &Ray,
        outward_normal: &Vec3,
        material: Material,
        u: f64,
        v: f64,
    ) -> Self {
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
            u,
            v,
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

    pub fn get_uv(p: &Point) -> (f64, f64) {
        // p: a given point on the sphere of radius one, centered at the origin.
        // u: returned value [0,1] of angle around the Y axis from X=-1.
        // v: returned value [0,1] of angle from Y=-1 to Y=+1.
        //     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
        //     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
        //     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

        let theta = (-p[1]).acos();
        let phi = (-p[2]).atan2(p[0]) + std::f64::consts::PI;

        let u = phi / (2.0 * std::f64::consts::PI);
        let v = theta / std::f64::consts::PI;
        (u, v)
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
        let (u, v) = Self::get_uv(&outward_normal);
        // NOTE -- we use dyn_clone here because self.material is a trait object -- you cannot
        // clone a trait object
        Some(HitRecord::new(
            p,
            t,
            ray,
            &outward_normal,
            dyn_clone::clone_box(&*self.material),
            u,
            v,
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
        let (u, v) = Sphere::get_uv(&outward_normal);
        // NOTE -- we use dyn_clone here because self.material is a trait object -- you cannot
        // clone a trait object
        Some(HitRecord::new(
            p,
            t,
            ray,
            &outward_normal,
            dyn_clone::clone_box(&*self.material),
            u,
            v,
        ))
    }

    fn try_bounding_box(&self, time0: f64, time1: f64) -> Option<Aabb> {
        let v = Vec3::new(self.radius, self.radius, self.radius);

        let box0 = Aabb::new(self.center(time0) - v, self.center(time0) + v);
        let box1 = Aabb::new(self.center(time1) - v, self.center(time1) + v);

        Some(box0.surrounding_box(&box1))
    }
}

/// XY Rectangle
pub struct XyRectangle {
    material: Material,
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    k: f64,
}
impl XyRectangle {
    pub fn new(material: Material, x0: f64, x1: f64, y0: f64, y1: f64, k: f64) -> Self {
        Self {
            material,
            x0,
            x1,
            y0,
            y1,
            k,
        }
    }
}
impl Hittable for XyRectangle {
    fn try_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - ray.orig[2]) / ray.dir[2];
        if t < t_min || t > t_max {
            return None;
        }
        let x = ray.orig[0] + t * ray.dir[0];
        let y = ray.orig[1] + t * ray.dir[1];
        if x < self.x0 || x > self.x1 || y < self.y0 || y > self.y1 {
            return None;
        }
        let u = (x - self.x0) / (self.x1 - self.x0);
        let v = (y - self.y0) / (self.y1 - self.y0);
        let outward_normal = Vec3::new(0.0, 0.0, 1.0);
        let p = ray.get(t);
        Some(HitRecord::new(
            p,
            t,
            ray,
            &outward_normal,
            dyn_clone::clone_box(&*self.material),
            u,
            v,
        ))
    }

    fn try_bounding_box(&self, _time0: f64, _time1: f64) -> Option<Aabb> {
        Some(Aabb::new(
            Point::new(self.x0, self.y0, self.k - 0.0001),
            Point::new(self.x1, self.y1, self.k + 0.0001),
        ))
    }
}
