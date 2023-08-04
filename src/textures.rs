//! Textures/colors

use dyn_clone::DynClone;
use rand::Rng;

use crate::{utils, Color, Point, Vec3};

pub type Texture = Box<dyn Textured + Send + Sync>;

pub trait Textured: DynClone {
    fn value(&self, u: f64, v: f64, p: &Point) -> Color;
}

/// Solid Color
#[derive(Debug, Clone)]
pub struct SolidColor {
    color_value: Color,
}
impl SolidColor {
    pub fn new(color: Color) -> Self {
        Self { color_value: color }
    }
}
impl Textured for SolidColor {
    fn value(&self, _u: f64, _v: f64, _p: &Point) -> Color {
        self.color_value
    }
}

/// Checker Texture
pub struct Checker {
    odd: Texture,
    even: Texture,
}
impl Checker {
    pub fn new(even: Texture, odd: Texture) -> Self {
        Self { odd, even }
    }

    pub fn from_solid_colors(c1: Color, c2: Color) -> Self {
        Self {
            even: Box::new(SolidColor::new(c1)),
            odd: Box::new(SolidColor::new(c2)),
        }
    }
}
impl Clone for Checker {
    fn clone(&self) -> Self {
        Self {
            odd: dyn_clone::clone_box(&*self.odd),
            even: dyn_clone::clone_box(&*self.even),
        }
    }
}
impl Textured for Checker {
    fn value(&self, u: f64, v: f64, p: &Point) -> Color {
        let sines = (10.0 * p[0]).sin() * (10.0 * p[1]).sin() * (10.0 * p[2]).sin();
        if sines < 0.0 {
            self.odd.value(u, v, p)
        } else {
            self.even.value(u, v, p)
        }
    }
}

/// Noise Texture
#[derive(Default, Clone)]
pub struct Noise {
    noise: Perlin,
    scale: f64,
}
impl Noise {
    pub fn new(scale: f64) -> Self {
        Self {
            noise: Default::default(),
            scale,
        }
    }
}
impl Textured for Noise {
    fn value(&self, _u: f64, _v: f64, p: &Point) -> Color {
        // Color::new(1.0, 1.0, 1.0) * 0.5 * (1.0 + self.noise.noise(&(self.scale * p)))
        // Color::new(1.0, 1.0, 1.0) * self.noise.turbulence(&(self.scale * p), None)
        Color::new(1.0, 1.0, 1.0)
            * 0.5
            * (1.0
                + (self.scale * p[2] + 10.0 * self.noise.turbulence(&(self.scale * p), None)).sin())
    }
}

/// Perlin Texture
#[derive(Debug, Clone)]
struct Perlin {
    ranvec: Vec<Vec3>,
    perm_x: Vec<i32>,
    perm_y: Vec<i32>,
    perm_z: Vec<i32>,
}
impl Default for Perlin {
    fn default() -> Self {
        let mut ranfloat = Vec::with_capacity(Self::POINT_COUNT);
        for _ in 0..Self::POINT_COUNT {
            ranfloat.push(utils::gen_random(3, Some(-1.0), Some(1.0)));
        }

        let perm_x = Self::perlin_generate_perm();
        let perm_y = Self::perlin_generate_perm();
        let perm_z = Self::perlin_generate_perm();

        Self {
            ranvec: ranfloat,
            perm_x,
            perm_y,
            perm_z,
        }
    }
}
impl Perlin {
    pub fn turbulence(&self, p: &Point, depth: Option<u32>) -> f64 {
        let depth = depth.unwrap_or(7);
        let mut accum = 0.0;
        let mut temp_p = *p;
        let mut weight = 1.0;

        for _ in 0..depth {
            accum += weight * self.noise(&temp_p);
            weight *= 0.5;
            temp_p *= 2.0;
        }
        accum.abs()
    }

    pub fn noise(&self, p: &Point) -> f64 {
        let u = p[0] - p[0].floor();
        let v = p[1] - p[1].floor();
        let w = p[2] - p[2].floor();

        let i = p[0].floor() as i32;
        let j = p[1].floor() as i32;
        let k = p[2].floor() as i32;

        // Initialize as zeros
        let mut c: [[[Vec3; 2]; 2]; 2] = [[[Vec3::zeros(); 2]; 2]; 2];
        for (di, c0) in c.iter_mut().enumerate() {
            for (dj, c1) in c0.iter_mut().enumerate() {
                for (dk, c2) in c1.iter_mut().enumerate() {
                    *c2 = self.ranvec[(self.perm_x[((i + di as i32) & 255) as usize]
                        ^ self.perm_y[((j + dj as i32) & 255) as usize]
                        ^ self.perm_z[((k + dk as i32) & 255) as usize])
                        as usize]
                }
            }
        }

        // Self::trilinear_interp(c, u, v, w)
        Self::perlin_interp(c, u, v, w)
    }

    const POINT_COUNT: usize = 256;

    fn perlin_generate_perm() -> Vec<i32> {
        let mut p = Vec::with_capacity(Self::POINT_COUNT);

        for i in 0..Self::POINT_COUNT {
            p.push(i as i32);
        }
        Self::permute(&mut p, Self::POINT_COUNT);
        p
    }

    fn permute(p: &mut [i32], n: usize) {
        let mut rng = rand::thread_rng();
        for i in (1..n).rev() {
            let target = rng.gen_range(0..=i);
            p.swap(i, target);
        }
    }

    fn _trilinear_interp(c: [[[f64; 2]; 2]; 2], u: f64, v: f64, w: f64) -> f64 {
        let mut accum = 0.0;
        for i in 0..2 {
            let i = i as f64;
            for j in 0..2 {
                let j = j as f64;
                for k in 0..2 {
                    let k = k as f64;
                    accum += (i * u + (1.0 - i) * (1.0 - u))
                        * (j * v + (1.0 - j) * (1.0 - v))
                        * (k * w + (1.0 - k) * (1.0 - w))
                        * c[i as usize][j as usize][k as usize];
                }
            }
        }
        accum
    }

    fn perlin_interp(c: [[[Vec3; 2]; 2]; 2], u: f64, v: f64, w: f64) -> f64 {
        // Hermitian Smoothing
        let uu = u.powi(2) * (3.0 - 2.0 * u);
        let vv = v.powi(2) * (3.0 - 2.0 * v);
        let ww = w.powi(2) * (3.0 - 2.0 * w);

        let mut accum = 0.0;

        for i in 0..2 {
            let i = i as f64;
            for j in 0..2 {
                let j = j as f64;
                for k in 0..2 {
                    let k = k as f64;
                    let weight_v = Vec3::new(u - i, v - j, w - k);
                    accum += (i * uu + (1.0 - i) * (1.0 - uu))
                        * (j * vv + (1.0 - j) * (1.0 - vv))
                        * (k * ww + (1.0 - k) * (1.0 - ww))
                        * c[i as usize][j as usize][k as usize].dot(&weight_v);
                }
            }
        }
        accum
    }
}
