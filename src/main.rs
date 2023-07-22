//! Code to generate a ray tracer, working through the examples
use std::sync::{Arc, RwLock};

use indicatif::ProgressBar;
use nalgebra::Vector3;
use rand::Rng;
use std::sync::mpsc::channel;
use threadpool::ThreadPool;

type Vec3 = Vector3<f64>;
type Point = Vec3;
type Color = Vec3;

fn write_color(color: &Color, samples_per_pixel: usize) {
    let scale = 1.0 / samples_per_pixel as f64;
    let r = (256.0 * (color[0] * scale).min(0.999).max(0.0)) as u64;
    let g = (256.0 * (color[1] * scale).min(0.999).max(0.0)) as u64;
    let b = (256.0 * (color[2] * scale).min(0.999).max(0.0)) as u64;
    println!("{r} {g} {b}");
}

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
        if depth <= 0 {
            return Color::zeros();
        }

        if let Some(HitRecord {
            p,
            normal,
            t: _,
            front_face: _,
        }) = obj.try_hit(self, 0.0, f64::MAX)
        {
            // Redirect randomly because of material
            let target = p + normal + utils::random_in_unit_sphere();
            let new_ray = Ray::new(p, target - p);
            return 0.5 * new_ray.get_color(obj, depth - 1);
        }
        let unit_direction = self.dir.normalize();
        let t = 0.5 * (unit_direction[1] + 1.0);
        (1.0 - t) * Color::new(1.0, 1.0, 1.0) + t * Color::new(0.5, 0.7, 1.0)
    }
}

#[derive(Debug)]
pub struct Sphere {
    pub center: Point,
    pub radius: f64,
}
impl Sphere {
    pub fn new(center: Point, radius: f64) -> Self {
        Self { center, radius }
    }

    pub fn hit(&self, ray: &Ray) -> f64 {
        let oc = ray.orig - self.center;
        let a = ray.dir.norm().powi(2);
        let half_b = oc.dot(&ray.dir);
        let c = oc.norm().powi(2) - self.radius.powi(2);
        let discriminant = half_b.powi(2) - a * c;
        if discriminant < 0.0 {
            -1.0
        } else {
            (-half_b - discriminant.sqrt()) / a
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
        Some(HitRecord::new(p, t, ray, &outward_normal))
    }
}

/// Represents a hit
#[derive(Debug)]
pub struct HitRecord {
    /// Point of intersection
    pub p: Point,
    /// Normali vector
    pub normal: Vec3,
    /// roots of intersection
    pub t: f64,
    /// Wither we are facing the normal
    front_face: bool,
}
impl HitRecord {
    pub fn new(p: Point, t: f64, ray: &Ray, outward_normal: &Vec3) -> Self {
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
        }
    }
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

pub trait Hittable {
    fn try_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}

/// Camera and related tasks
#[derive(Debug)]
pub struct Camera {
    origin: Point,
    lower_left_corner: Point,
    horizontal: Vec3,
    vertical: Vec3,
}
impl Default for Camera {
    fn default() -> Self {
        let aspect_ratio = 16.0 / 9.0;
        let viewport_height = 2.0;
        let viewport_width = aspect_ratio * viewport_height;
        let focal_length = 1.0;

        let origin = Point::zeros();
        let horizontal = Vec3::new(viewport_width, 0.0, 0.0);
        let vertical = Vec3::new(0.0, viewport_height, 0.0);
        let lower_left_corner =
            origin - horizontal / 2.0 - vertical / 2.0 - Vec3::new(0.0, 0.0, focal_length);
        Self {
            origin,
            lower_left_corner,
            horizontal,
            vertical,
        }
    }
}
impl Camera {
    pub fn get_ray(&self, u: f64, v: f64) -> Ray {
        Ray::new(
            self.origin,
            self.lower_left_corner + u * self.horizontal + v * self.vertical - self.origin,
        )
    }
}

/// Module for random utils
pub mod utils {
    use nalgebra::Vector3;
    use rand::Rng;
    type Vec3 = Vector3<f64>;

    /// Compute a random vector inside the unit circle
    ///
    /// Randomly generate vectors. If the norm is < 1, it is inside the unit circle
    pub fn random_in_unit_sphere() -> Vec3 {
        let mut rng = rand::thread_rng();
        loop {
            let p = Vec3::from_vec(
                (0..3)
                    .into_iter()
                    .map(|_| rng.gen_range(-1.0..1.0))
                    .collect(),
            );
            if p.norm().powi(2) < 1.0 {
                return p;
            }
        }
    }
}

fn main() {
    // Setup Image
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    let image_width: usize = 400;
    let image_height: usize = (image_width as f64 / ASPECT_RATIO).round() as usize;
    let samples_per_pixel = 20;
    let max_depth = 50;

    // Create World
    let mut world = HittableList::default();
    world.add(Box::new(Sphere::new(Point::new(0.0, 0.0, -1.0), 0.5)));
    world.add(Box::new(Sphere::new(Point::new(0.0, -100.5, -1.0), 100.0)));
    let protected_world = Arc::new(RwLock::new(world));

    // Create the camera
    let cam = Camera::default();

    // Random generator
    let mut rng = rand::thread_rng();

    // Create the threadpool
    let n_workers = 20;
    let pool = ThreadPool::new(n_workers);

    // Render
    print!("P3\n{image_width} {image_height}\n255\n");
    let bar = ProgressBar::new((image_width * image_height) as u64);
    for j in (0..=image_height - 1).rev() {
        for i in 0..image_width {
            let (worker_tx, worker_rx) = channel();
            for _ in 0..samples_per_pixel {
                let tx = worker_tx.clone();
                let u = (i as f64 + rng.gen::<f64>()) / (image_width - 1) as f64;
                let v = (j as f64 + rng.gen::<f64>()) / (image_height - 1) as f64;
                let ray = cam.get_ray(u, v);
                let protected_world = Arc::clone(&protected_world);
                pool.execute(move || {
                    tx.send(ray.get_color(&*protected_world.read().unwrap(), max_depth))
                        .unwrap();
                });
            }
            let pixel_color = worker_rx
                .iter()
                .take(samples_per_pixel)
                .reduce(|a, b| a + b)
                .unwrap();

            bar.inc(1);
            write_color(&pixel_color, samples_per_pixel);
        }
    }
    bar.finish()
}
