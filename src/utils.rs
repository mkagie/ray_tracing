//! Utils

use nalgebra::Vector3;
use rand::Rng;
type Vec3 = Vector3<f64>;
use super::Color;
use image::Rgb;

/// Compute a random vector inside the unit circle
///
/// Randomly generate vectors. If the norm is < 1, it is inside the unit circle.
pub fn random_in_unit_sphere() -> Vec3 {
    loop {
        let p = gen_random(3, Some(-1.0), Some(1.0));

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

pub fn get_pixel(color: &Color, samples_per_pixel: usize) -> Rgb<u8> {
    let scale = 1.0 / samples_per_pixel as f64;

    // Divide the color by the number of samples and gamma-correct for gamma = 2.0
    let r = scale_color((scale * color[0]).sqrt());
    let g = scale_color((scale * color[1]).sqrt());
    let b = scale_color((scale * color[2]).sqrt());

    Rgb([r, g, b])
}

/// scale the color to between 0 and 255
fn scale_color(val: f64) -> u8 {
    (256.0 * val.min(0.999).max(0.0)) as u8
}

/// Generate Random Vectors
pub fn gen_random(len: usize, min: Option<f64>, max: Option<f64>) -> Vec3 {
    let mut rng = rand::thread_rng();
    Vec3::from_vec(
        (0..len)
            .map(|_| {
                if min.is_some() && max.is_some() {
                    rng.gen_range(min.unwrap()..max.unwrap())
                } else {
                    rng.gen()
                }
            })
            .collect(),
    )
}
