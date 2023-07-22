//! Ray Tracing Library

use nalgebra::Vector3;
use rand::Rng;

pub mod objects;
use objects::Hittable;

pub mod materials;
use materials::Scatterable;

pub mod utils;

pub type Vec3 = Vector3<f64>;
pub type Point = Vec3;
pub type Color = Vec3;
pub type Material = Box<dyn Scatterable + Send + Sync>;

/// Prelude
pub mod prelude {
    pub use crate::materials::{Dielectric, Lambertian, Metal};
    pub use crate::objects::{HittableList, Sphere};
    pub use crate::{Camera, Color, Material, Point, Vec3};
}

/// The ray in ray tracing
#[derive(Debug)]
pub struct Ray {
    pub orig: Point,
    pub dir: Vec3,
}
impl Ray {
    pub fn new(orig: Point, dir: Vec3) -> Self {
        Self { orig, dir }
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

// TODO(mkagie) Convert the camera from look at and look from to look_from and quaternion
/// Camera and related tasks
#[derive(Debug)]
pub struct Camera {
    origin: Point,
    lower_left_corner: Point,
    horizontal: Vec3,
    vertical: Vec3,
    u: Vec3,
    v: Vec3,
    _w: Vec3,
    lens_radius: f64,
}
impl Camera {
    // TODO(mkagie) use uom for angle
    pub fn new(
        look_from: Point,
        look_at: Point,
        v_up: Vec3,
        vertical_fov_deg: f64,
        aspect_ratio: f64,
        aperture: f64,
        focus_dist: f64,
    ) -> Self {
        // Establish the viewport
        let theta = vertical_fov_deg.to_radians();
        let h = (theta / 2.0).tan();
        let viewport_height = 2.0 * h;
        let viewport_width = aspect_ratio * viewport_height;

        // Calculate the viewing vectors
        let w = (look_from - look_at).normalize();
        let u = (v_up.cross(&w)).normalize();
        let v = w.cross(&u);

        let origin = look_from;
        let horizontal = focus_dist * viewport_width * u;
        let vertical = focus_dist * viewport_height * v;
        let lower_left_corner = origin - horizontal / 2.0 - vertical / 2.0 - focus_dist * w;

        let lens_radius = aperture / 2.0;

        Self {
            origin,
            lower_left_corner,
            horizontal,
            vertical,
            u,
            v,
            _w: w,
            lens_radius,
        }
    }

    pub fn get_ray(&self, s: f64, t: f64) -> Ray {
        let rd = self.lens_radius * Self::random_in_unit_disk();
        let offset = self.u * rd[0] + self.v * rd[1];

        Ray::new(
            self.origin + offset,
            self.lower_left_corner + s * self.horizontal + t * self.vertical - self.origin - offset,
        )
    }

    /// Generate a random vector inside a unit disk
    /// This simulates defocus blur
    fn random_in_unit_disk() -> Vec3 {
        let mut rng = rand::thread_rng();
        loop {
            let mut random_vec: Vec<f64> = (0..2).map(|_| rng.gen_range(-1.0..1.0)).collect();
            random_vec.push(0.0);
            let p = Vec3::from_vec(random_vec);
            if p.norm().powi(2) < 1.0 {
                return p;
            }
        }
    }
}
