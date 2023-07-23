//! Ray Tracing Library

use nalgebra::Vector3;

pub mod objects;
use objects::Hittable;

pub mod materials;
use materials::Scatterable;

pub mod utils;

pub mod cameras;

pub mod aabb;

pub type Vec3 = Vector3<f64>;
pub type Point = Vec3;
pub type Color = Vec3;
pub type Material = Box<dyn Scatterable + Send + Sync>;

/// Prelude
pub mod prelude {
    pub use crate::cameras::{Camera, CameraConfig};
    pub use crate::materials::{Dielectric, Lambertian, Metal};
    pub use crate::objects::{HittableList, HittableListConfig, MovingSphere, Sphere};
    pub use crate::{Color, Material, Point, Vec3};
}

/// The ray in ray tracing
#[derive(Debug)]
pub struct Ray {
    pub orig: Point,
    pub dir: Vec3,
    pub time: f64,
}
impl Ray {
    pub fn new(orig: Point, dir: Vec3, time: f64) -> Self {
        Self { orig, dir, time }
    }

    pub fn get(&self, t: f64) -> Point {
        self.orig + t * self.dir
    }

    /// Linearly blends white and blue depending on height of y
    pub fn get_color(&self, obj: &impl Hittable, depth: u32) -> Color {
        // If we have exceeded the ray bounce limit, no more light is gathered
        if depth == 0 {
            return Color::zeros();
        }

        // Put a minimum of 0.001 to reduce shadow acne
        if let Some(hr) = obj.try_hit(self, 0.001, f64::MAX) {
            if let Some(sr) = hr.material.try_scatter(self, &hr) {
                return sr
                    .attenuation
                    .component_mul(&sr.scattered.get_color(obj, depth - 1));
            }
            return Color::zeros();
        }
        // TODO(mkagie) Have a default color passed in -- although this matches the sky nicely
        let unit_direction = self.dir.normalize();
        let t = 0.5 * (unit_direction[1] + 1.0);
        (1.0 - t) * Color::new(1.0, 1.0, 1.0) + t * Color::new(0.5, 0.7, 1.0)
    }
}
