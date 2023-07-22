//! Code to generate a ray tracer, working through the examples
use std::sync::{Arc, RwLock};

use dyn_clone::DynClone;
use indicatif::ProgressBar;
use nalgebra::Vector3;
use rand::Rng;
use std::sync::mpsc::channel;
use threadpool::ThreadPool;

type Vec3 = Vector3<f64>;
type Point = Vec3;
type Color = Vec3;
type Material = Box<dyn Scatterable + Send + Sync>;

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

        // Put a minimum of 0.001 to reduce shadow acne
        if let Some(hr) = obj.try_hit(self, 0.001, f64::MAX) {
            if let Some(sr) = hr.material.try_scatter(self, &hr) {
                return sr
                    .attenuation
                    .component_mul(&sr.scattered.get_color(obj, depth - 1));
            }
            return Color::zeros();
        }
        let unit_direction = self.dir.normalize();
        let t = 0.5 * (unit_direction[1] + 1.0);
        (1.0 - t) * Color::new(1.0, 1.0, 1.0) + t * Color::new(0.5, 0.7, 1.0)
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

/// Represents a hit
pub struct HitRecord {
    /// Point of intersection
    pub p: Point,
    /// Normali vector
    pub normal: Vec3,
    /// roots of intersection
    pub t: f64,
    /// Wither we are facing the normal
    front_face: bool,
    // TODO(mkagie) I'm not sure I like having the material in a HitRecord
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

/// Scatter Result
#[derive(Debug)]
pub struct ScatterResult {
    /// Attenuation Color
    attenuation: Color,
    /// Resulting Scattered Ray
    scattered: Ray,
}

/// Material
pub trait Scatterable: DynClone {
    fn try_scatter(&self, ray_in: &Ray, hit_record: &HitRecord) -> Option<ScatterResult>;
}

/// Lmabertian Scatterer
#[derive(Debug, Clone)]
pub struct Lambertian {
    albedo: Color,
}
impl Lambertian {
    pub fn new(albedo: Color) -> Self {
        Self { albedo }
    }
}
impl Scatterable for Lambertian {
    fn try_scatter(&self, ray_in: &Ray, hit_record: &HitRecord) -> Option<ScatterResult> {
        let mut scatter_direction = hit_record.normal + utils::random_in_unit_sphere();

        // Protect against if hit_record.normal and the random_in_unit_sphere as exact opposites
        if scatter_direction.norm() < 1e-8 {
            scatter_direction = hit_record.normal;
        }
        let scattered = Ray::new(hit_record.p, scatter_direction);
        let attenuation = self.albedo;
        Some(ScatterResult {
            attenuation,
            scattered,
        })
    }
}

/// Metal Scatterer
#[derive(Debug, Clone)]
pub struct Metal {
    albedo: Color,
    fuzz: f64,
}
impl Metal {
    pub fn new(albedo: Color, fuzz: f64) -> Self {
        Self { albedo, fuzz }
    }

    fn reflect(v: &Vec3, n: &Vec3) -> Vec3 {
        v - 2.0 * v.dot(n) * n
    }
}
impl Scatterable for Metal {
    fn try_scatter(&self, ray_in: &Ray, hit_record: &HitRecord) -> Option<ScatterResult> {
        let reflected = Self::reflect(&ray_in.dir.normalize(), &hit_record.normal);
        let scattered = Ray::new(
            hit_record.p,
            reflected + self.fuzz * utils::random_in_unit_sphere(),
        );
        let attenuation = self.albedo;
        if scattered.dir.dot(&hit_record.normal) > 0.0 {
            Some(ScatterResult {
                attenuation,
                scattered,
            })
        } else {
            None
        }
    }
}

/// A Dielectric is a refractive material, such as glass
#[derive(Debug, Clone)]
pub struct Dielectric {
    ir: f64,
}
impl Dielectric {
    pub fn new(ir: f64) -> Self {
        Self { ir }
    }

    fn refract(uv: &Vec3, n: &Vec3, etai_over_etat: f64) -> Vec3 {
        let cos_theta = -uv.dot(&n).min(1.0);
        let r_out_perp = etai_over_etat * (uv + cos_theta * n);
        let r_out_parallel = -(1.0 - r_out_perp.norm().powi(2)).abs().sqrt() * n;
        r_out_perp + r_out_parallel
    }

    fn reflectance(cosine: f64, ref_idx: f64) -> f64 {
        // Use Schlick's approximation for reflectance
        let r0 = ((1.0 - ref_idx) / (1.0 + ref_idx)).powi(2);
        r0 + (1.0 - r0) * (1.0 - cosine).powi(5)
    }
}
impl Scatterable for Dielectric {
    fn try_scatter(&self, ray_in: &Ray, hit_record: &HitRecord) -> Option<ScatterResult> {
        let attenuation = Color::new(1.0, 1.0, 1.0);
        let refraction_ratio = if hit_record.front_face {
            1.0 / self.ir
        } else {
            self.ir
        };

        let unit_direction = ray_in.dir.normalize();
        let cos_theta = -unit_direction.dot(&hit_record.normal).min(1.0);
        let sin_theta = (1.0 - cos_theta.powi(2)).sqrt();

        let cannot_refract = refraction_ratio * sin_theta > 1.0;
        let mut rng = rand::thread_rng();

        let direction = if cannot_refract
            || Self::reflectance(cos_theta, refraction_ratio) > rng.gen::<f64>()
        {
            // TODO(mkagie) Remove reflect and refract from Metal and Dielectric, since they are no
            // longer tied to them
            Metal::reflect(&unit_direction, &hit_record.normal)
        } else {
            Self::refract(&unit_direction, &hit_record.normal, refraction_ratio)
        };

        let scattered = Ray::new(hit_record.p, direction);
        Some(ScatterResult {
            attenuation,
            scattered,
        })
    }
}

/// Module for random utils
pub mod utils {
    use nalgebra::Vector3;
    use rand::Rng;
    type Vec3 = Vector3<f64>;
    use super::Color;

    /// Compute a random vector inside the unit circle
    ///
    /// Randomly generate vectors. If the norm is < 1, it is inside the unit circle.
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
                // We normalize the vector on the output to more exactly represent Lambertian
                return p.normalize();
            }
        }
    }

    /// A more intuitive and easy-to-understand diffuse method
    ///
    /// A uniform scatter direction for all angles away from the hit point
    pub fn random_in_hemisphere(normal: &Vec3) -> Vec3 {
        let in_unit_sphere = random_in_unit_sphere();
        if in_unit_sphere.dot(normal) > 0.0 {
            in_unit_sphere
        } else {
            -in_unit_sphere
        }
    }

    pub fn write_color(color: &Color, samples_per_pixel: usize) {
        let scale = 1.0 / samples_per_pixel as f64;

        // Divide the color by the number of samples and gamma-correct for gamma = 2.0
        let r = scale_color((scale * color[0]).sqrt());
        let g = scale_color((scale * color[1]).sqrt());
        let b = scale_color((scale * color[2]).sqrt());

        println!("{r} {g} {b}");
    }

    /// scale the color to between 0 and 255
    fn scale_color(val: f64) -> u64 {
        (256.0 * val.min(0.999).max(0.0)) as u64
    }
}

fn main() {
    // Setup Image
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    let image_width: usize = 400;
    let image_height: usize = (image_width as f64 / ASPECT_RATIO).round() as usize;
    let samples_per_pixel = 20;
    let max_depth = 20;

    // Create World
    let mut world = HittableList::default();

    // Define materials
    let material_ground = Box::new(Lambertian::new(Color::new(0.8, 0.8, 0.0)));
    // let material_center = Box::new(Lambertian::new(Color::new(0.7, 0.3, 0.3)));
    let material_center = Box::new(Lambertian::new(Color::new(0.1, 0.2, 0.5)));
    // let material_center = Box::new(Dielectric::new(1.5));
    // let material_left = Box::new(Metal::new(Color::new(0.8, 0.8, 0.8), 0.3));
    let material_left = Box::new(Dielectric::new(1.5));
    let material_right = Box::new(Metal::new(Color::new(0.8, 0.6, 0.2), 0.0));

    world.add(Box::new(Sphere::new(
        Point::new(0.0, -100.5, -1.0),
        100.0,
        material_ground,
    )));
    world.add(Box::new(Sphere::new(
        Point::new(0.0, 0.0, -1.0),
        0.5,
        material_center,
    )));
    world.add(Box::new(Sphere::new(
        Point::new(-1.0, 0.0, -1.0),
        0.5,
        material_left.clone(),
    )));
    world.add(Box::new(Sphere::new(
        Point::new(-1.0, 0.0, -1.0),
        -0.4,
        material_left,
    )));
    world.add(Box::new(Sphere::new(
        Point::new(1.0, 0.0, -1.0),
        0.5,
        material_right,
    )));
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
            utils::write_color(&pixel_color, samples_per_pixel);
        }
    }
    bar.finish()
}
