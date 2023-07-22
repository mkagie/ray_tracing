//! Code to generate a ray tracer, working through the examples
use indicatif::ProgressBar;
use nalgebra::Vector3;
use rand::Rng;
use rayon::prelude::*;

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
    pub fn get_color(&self, obj: &impl Hittable) -> Color {
        if let Some(HitRecord {
            p,
            normal,
            t,
            front_face,
        }) = obj.try_hit(&self, 0.0, f64::MAX)
        {
            return 0.5 * Color::new(normal[0] + 1.0, normal[1] + 1.0, normal[2] + 1.0);
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
pub struct HittableList(Vec<Box<dyn Hittable>>);
impl HittableList {
    pub fn add(&mut self, boxed_obj: Box<dyn Hittable>) {
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
            self.origin.clone(),
            self.lower_left_corner + u * self.horizontal + v * self.vertical - self.origin,
        )
    }
}

fn main() {
    // Setup Image
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    let image_width: usize = 400;
    let image_height: usize = (image_width as f64 / ASPECT_RATIO).round() as usize;
    let samples_per_pixel = 20;

    // Setup Camera
    let viewport_height = 2.0;
    let viewport_width = ASPECT_RATIO * viewport_height;
    let focal_length = 1.0;

    // Scene
    let origin = Point::new(0.0, 0.0, 0.0);
    let horizontal = Vec3::new(viewport_width, 0.0, 0.0);
    let vertical = Vec3::new(0.0, viewport_height, 0.0);
    let lower_left_corner =
        origin - horizontal / 2.0 - vertical / 2.0 - Vec3::new(0.0, 0.0, focal_length);

    // Create World
    let mut world = HittableList::default();
    world.add(Box::new(Sphere::new(Point::new(0.0, 0.0, -1.0), 0.5)));
    world.add(Box::new(Sphere::new(Point::new(0.0, -100.5, -1.0), 100.0)));

    // Create the camera
    let cam = Camera::default();

    // Random sampler
    let mut rng = rand::thread_rng();

    // Render
    print!("P3\n{image_width} {image_height}\n255\n");
    let bar = ProgressBar::new((image_width * image_height * samples_per_pixel) as u64);
    for j in (0..=image_height - 1).rev() {
        for i in 0..image_width {
            let mut pixel_color = Color::zeros();
            for s in 0..samples_per_pixel {
                let u = (i as f64 + rng.gen::<f64>()) / (image_width - 1) as f64;
                let v = (j as f64 + rng.gen::<f64>()) / (image_height - 1) as f64;
                let ray = cam.get_ray(u, v);
                pixel_color += ray.get_color(&world);
                bar.inc(1);
            }
            write_color(&pixel_color, samples_per_pixel);
        }
    }
    bar.finish()
}
