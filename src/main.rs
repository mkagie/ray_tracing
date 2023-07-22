//! Code to generate a ray tracer, working through the examples
use indicatif::ProgressBar;
use nalgebra::Vector3;

type Vec3 = Vector3<f64>;
type Point = Vec3;
type Color = Vec3;

fn write_color(color: &Color) {
    let ir = (255.999 * color[0]).round() as usize;
    let ig = (255.999 * color[1]).round() as usize;
    let ib = (255.999 * color[2]).round() as usize;
    println!("{ir} {ig} {ib}");
}

#[derive(Debug)]
struct Ray {
    orig: Point,
    dir: Vec3,
}
impl Ray {
    pub fn new(orig: &Point, dir: &Vec3) -> Self {
        Self {
            orig: orig.clone(),
            dir: dir.clone(),
        }
    }

    pub fn get(&self, t: f64) -> Point {
        self.orig + t * self.dir
    }

    pub fn get_color(&self) -> Color {
        // Blends white and blue
        let unit_direction = self.dir.normalize();
        let t = 0.5 * (unit_direction[1] + 1.0);
        (1.0 - t) * Color::new(1.0, 1.0, 1.0) + t * Color::new(0.5, 0.7, 1.0)
    }
}

fn main() {
    // Setup Image
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    let image_width: usize = 400;
    let image_height: usize = (image_width as f64 / ASPECT_RATIO).round() as usize;

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

    // Render
    print!("P3\n{image_width} {image_height}\n255\n");

    let bar = ProgressBar::new((image_width * image_height) as u64);
    for j in (0..=image_height - 1).rev() {
        for i in 0..image_width {
            let u = i as f64 / (image_width - 1) as f64;
            let v = j as f64 / (image_height - 1) as f64;
            let dir = lower_left_corner + u * horizontal + v * vertical - origin;
            let ray = Ray::new(&origin, &dir);
            let color = ray.get_color();
            write_color(&color);
            bar.inc(1);
        }
    }
    bar.finish()
}
