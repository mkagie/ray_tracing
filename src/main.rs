//! Code to generate a ray tracer, working through the examples
use std::sync::{Arc, RwLock};

use indicatif::{ProgressBar, ProgressStyle};
use rand::Rng;
use ray_tracing::{prelude::*, utils};
use std::sync::mpsc::channel;
use threadpool::ThreadPool;

fn random_scene() -> HittableList {
    let mut world = HittableList::default();
    let ground_material = Box::new(Lambertian::new(Color::new(0.5, 0.5, 0.5)));
    world.add(Box::new(Sphere::new(
        Point::new(0.0, -1000.0, 0.0),
        1000.0,
        ground_material,
    )));

    let mut rng = rand::thread_rng();
    for a in -11..11 {
        for b in -11..11 {
            let choose_mat: f64 = rng.gen();
            let center = Point::new(
                a as f64 + 0.9 * rng.gen::<f64>(),
                0.2,
                b as f64 + 0.9 * rng.gen::<f64>(),
            );

            if (center - Point::new(4.0, 0.2, 0.0)).norm() > 0.9 {
                if choose_mat < 0.8 {
                    //difuse
                    let albedo = utils::gen_random(3, None, None)
                        .component_mul(&utils::gen_random(3, None, None));
                    let sphere_material = Box::new(Lambertian::new(albedo));
                    world.add(Box::new(Sphere::new(center, 0.2, sphere_material)));
                } else if choose_mat < 0.95 {
                    // Metal
                    let albedo = utils::gen_random(3, Some(0.5), Some(1.0));
                    let fuzz = rng.gen_range(0.0..0.5);
                    let sphere_material = Box::new(Metal::new(albedo, fuzz));
                    world.add(Box::new(Sphere::new(center, 0.2, sphere_material)));
                } else {
                    // Glass
                    let sphere_material = Box::new(Dielectric::new(1.5));
                    world.add(Box::new(Sphere::new(center, 0.2, sphere_material)));
                }
            }
        }
    }
    let mat0 = Box::new(Dielectric::new(1.5));
    world.add(Box::new(Sphere::new(Point::new(0.0, 1.0, 0.0), 1.0, mat0)));
    let mat1 = Box::new(Lambertian::new(Color::new(0.4, 0.2, 0.1)));
    world.add(Box::new(Sphere::new(Point::new(-4.0, 1.0, 0.0), 1.0, mat1)));
    let mat2 = Box::new(Metal::new(Color::new(0.7, 0.6, 0.5), 0.0));
    world.add(Box::new(Sphere::new(Point::new(4.0, 1.0, 0.0), 1.0, mat2)));

    world
}

fn orig_scene() -> HittableList {
    // Create World
    let mut world = HittableList::default();
    // let r = (std::f64::consts::PI / 4.0).cos();
    // Define materials
    let material_ground = Box::new(Lambertian::new(Color::new(0.8, 0.8, 0.0)));
    // let material_center = Box::new(Lambertian::new(Color::new(0.7, 0.3, 0.3)));
    let material_center = Box::new(Lambertian::new(Color::new(0.1, 0.2, 0.5)));
    // let material_center = Box::new(Dielectric::new(1.5));
    // let material_left = Box::new(Metal::new(Color::new(0.8, 0.8, 0.8), 0.3));
    let material_left = Box::new(Dielectric::new(1.5));
    // let material_left = Box::new(Lambertian::new(Color::new(0.0, 0.0, 1.0)));
    // let material_right = Box::new(Lambertian::new(Color::new(1.0, 0.0, 0.0)));
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
    //
    // world.add(Box::new(Sphere::new(
    //     Point::new(-r, 0.0, -1.0),
    //     r,
    //     material_left,
    // )));
    // world.add(Box::new(Sphere::new(
    //     Point::new(r, 0.0, -1.0),
    //     r,
    //     material_right,
    // )));
    // let world = random_scene();
    world
}

fn main() {
    // Setup Image
    const ASPECT_RATIO: f64 = 3.0 / 2.0;
    // const ASPECT_RATIO: f64 = 16.0 / 9.0;
    // let image_width: usize = 400;
    let image_width: usize = 1200;
    let image_height: usize = (image_width as f64 / ASPECT_RATIO).round() as usize;
    let samples_per_pixel = 500;
    let max_depth = 50;

    // Create World
    // let world = orig_scene();
    let world = random_scene();
    let protected_world = Arc::new(RwLock::new(world));

    let look_from = Point::new(13.0, 2.0, 3.0);
    let look_at = Point::new(0.0, 0.0, 0.0);
    let v_up = Vec3::new(0.0, 1.0, 0.0);
    let dist_to_focus = 10.0;
    let aperture = 0.1;
    let cam = Camera::new(
        look_from,
        look_at,
        v_up,
        20.0,
        ASPECT_RATIO,
        aperture,
        dist_to_focus,
    );

    // Random generator
    let mut rng = rand::thread_rng();

    // Create the threadpool
    let n_workers = 20; // Number of cores I have
    let pool = ThreadPool::new(n_workers);

    // Render
    print!("P3\n{image_width} {image_height}\n255\n");
    let bar = ProgressBar::new((image_width * image_height) as u64);
    bar.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise} --- {eta_precise} --- {pos}/{len}] {wide_bar}",
            // "[{elapsed_precise} {eta}] {bar:40.cyan/blue} {pos}/{len} {msg}",
        )
        .unwrap(),
    );
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
