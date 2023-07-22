//! Implementation of materials

use crate::{objects::HitRecord, utils, Color, Ray, Vec3};
use dyn_clone::DynClone;
use rand::Rng;

/// Material
pub trait Scatterable: DynClone {
    fn try_scatter(&self, ray_in: &Ray, hit_record: &HitRecord) -> Option<ScatterResult>;
}

/// Scatter Result
#[derive(Debug)]
pub struct ScatterResult {
    /// Attenuation Color
    pub attenuation: Color,
    /// Resulting Scattered Ray
    pub scattered: Ray,
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
    fn try_scatter(&self, _ray_in: &Ray, hit_record: &HitRecord) -> Option<ScatterResult> {
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
