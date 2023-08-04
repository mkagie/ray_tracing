//! Ray Tracing Library

use nalgebra::Vector3;

pub mod bvh;
pub mod cameras;
pub mod materials;
pub mod mediums;
pub mod objects;
pub mod textures;
pub mod transrot;
pub mod utils;

use materials::Scatterable;
use objects::Hittable;

pub type Vec3 = Vector3<f64>;
pub type Point = Vec3;
pub type Color = Vec3;
pub type Material = Box<dyn Scatterable + Send + Sync>;

/// Prelude
pub mod prelude {
    pub use crate::bvh::BvhNode;
    pub use crate::cameras::{Camera, CameraConfig};
    pub use crate::materials::{Dielectric, DiffuseLight, Lambertian, Metal};
    pub use crate::mediums::ConstantMedium;
    pub use crate::objects::{
        BoxObj, HittableList, HittableListConfig, HittableObj, MovingSphere, Rectangle,
        RectangleType, Sphere,
    };
    pub use crate::textures::{Checker, Noise};
    pub use crate::transrot::{RotateY, Translate};
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

    /// Get the location of the ray at a certain time
    pub fn get(&self, t: f64) -> Point {
        self.orig + t * self.dir
    }

    /// Get the color of the ray given a background color, list of things to hit, and a number of
    /// bounces allowed
    pub fn get_color(&self, background_color: &Color, obj: &impl Hittable, depth: u32) -> Color {
        // If we have exceeded the ray bounce limit, no more light is gathered
        if depth == 0 {
            return Color::zeros();
        }

        // Put a minimum of 0.001 to reduce shadow acne
        if let Some(hr) = obj.try_hit(self, 0.001, f64::MAX) {
            let emitted = hr.material.color_emitted(hr.u, hr.v, &hr.p);
            if let Some(sr) = hr.material.try_scatter(self, &hr) {
                return emitted
                    + sr.attenuation.component_mul(&sr.scattered.get_color(
                        background_color,
                        obj,
                        depth - 1,
                    ));
            }
            return emitted;
        }
        // Did not hit anything
        background_color.to_owned()
    }
}
