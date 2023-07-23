//! Textures/colors

use dyn_clone::DynClone;

use crate::{Color, Point};

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
