//! Cameras and configs for cameras
use crate::utils::SerdeVector;
use crate::{Point, Ray, Vec3};
use rand::Rng;
use serde::{Deserialize, Serialize};

// TODO(mkagie) Spend the time to deserialize this directly into Points
/// Camera Config
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CameraConfig {
    pub look_from: SerdeVector,
    pub look_at: SerdeVector,
    pub v_up: SerdeVector,
    pub vertical_fov_deg: f64,
    pub aspect_ratio: f64,
    pub aperture: f64,
    pub focus_distance: f64,
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
