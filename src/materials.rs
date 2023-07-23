//! Implementation of materials

use crate::{
    objects::HitRecord,
    texture::{SolidColor, Texture},
    utils::{self, SerdeVector},
    Color, Material, Ray,
};
use dyn_clone::DynClone;
use rand::Rng;
use serde::{Deserialize, Serialize};

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

/// Config for materials
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum MaterialConfig {
    Lambertian(LambertianConfig),
    Metal(MetalConfig),
    Dielectric(DielectricConfig),
}

/// Generator from config
pub struct Generator;
impl Generator {
    pub fn from_config(config: MaterialConfig) -> Material {
        match config {
            MaterialConfig::Lambertian(c) => Box::new(Lambertian::from_config(c)),
            MaterialConfig::Metal(c) => Box::new(Metal::from_config(c)),
            MaterialConfig::Dielectric(c) => Box::new(Dielectric::from_config(c)),
        }
    }
}

/// Lmabertian Scatterer
pub struct Lambertian {
    albedo: Texture,
}
impl Lambertian {
    pub fn new(albedo: Color) -> Self {
        Self {
            albedo: Box::new(SolidColor::new(albedo)),
        }
    }

    pub fn from_texture(texture: Texture) -> Self {
        Self { albedo: texture }
    }

    pub fn from_config(config: LambertianConfig) -> Self {
        Self::new(config.albedo.into())
    }
}
impl Scatterable for Lambertian {
    fn try_scatter(&self, ray_in: &Ray, hit_record: &HitRecord) -> Option<ScatterResult> {
        let mut scatter_direction = hit_record.normal + utils::random_in_unit_sphere();

        // Protect against if hit_record.normal and the random_in_unit_sphere as exact opposites
        if scatter_direction.norm() < 1e-8 {
            scatter_direction = hit_record.normal;
        }
        let scattered = Ray::new(hit_record.p, scatter_direction, ray_in.time);
        let attenuation = self.albedo.value(hit_record.u, hit_record.v, &hit_record.p);
        Some(ScatterResult {
            attenuation,
            scattered,
        })
    }
}
impl Clone for Lambertian {
    fn clone(&self) -> Self {
        Self {
            albedo: dyn_clone::clone_box(&*self.albedo),
        }
    }
}

/// Lambertian Config
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LambertianConfig {
    pub albedo: SerdeVector,
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

    pub fn from_config(config: MetalConfig) -> Self {
        Self {
            albedo: config.albedo.into(),
            fuzz: config.fuzz,
        }
    }
}
impl Scatterable for Metal {
    fn try_scatter(&self, ray_in: &Ray, hit_record: &HitRecord) -> Option<ScatterResult> {
        let reflected = utils::reflect(&ray_in.dir.normalize(), &hit_record.normal);
        let scattered = Ray::new(
            hit_record.p,
            reflected + self.fuzz * utils::random_in_unit_sphere(),
            ray_in.time,
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

/// Metal Config
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetalConfig {
    pub albedo: SerdeVector,
    fuzz: f64,
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

    pub fn from_config(config: DielectricConfig) -> Self {
        Self { ir: config.ir }
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
            utils::reflect(&unit_direction, &hit_record.normal)
        } else {
            utils::refract(&unit_direction, &hit_record.normal, refraction_ratio)
        };

        let scattered = Ray::new(hit_record.p, direction, ray_in.time);
        Some(ScatterResult {
            attenuation,
            scattered,
        })
    }
}

/// Dielectric Config
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DielectricConfig {
    pub ir: f64,
}
