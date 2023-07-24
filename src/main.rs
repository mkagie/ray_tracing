//! Code to generate a ray tracer, working through the examples
use std::sync::{Arc, RwLock};

use clap::Parser;
use image::RgbImage;
use indicatif::{ProgressBar, ProgressStyle};
use rand::Rng;
use ray_tracing::{prelude::*, utils};
use serde::{Deserialize, Serialize};
use std::sync::mpsc::channel;
use threadpool::ThreadPool;

fn random_scene() -> HittableList {
    let mut world = HittableList::default();
    // let ground_material = Box::new(Lambertian::new(Color::new(0.5, 0.5, 0.5)));
    // world.add(Box::new(Sphere::new(
    //     Point::new(0.0, -1000.0, 0.0),
    //     1000.0,
    //     ground_material,
    // )));

    let checker = Box::new(Checker::from_solid_colors(
        Color::new(0.2, 0.3, 0.1),
        Color::new(0.9, 0.9, 0.9),
    ));
    world.add(Box::new(Sphere::new(
        Point::new(0.0, -1000.0, 0.0),
        1000.0,
        Box::new(Lambertian::from_texture(checker)),
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
                    let center2 = center + Vec3::new(0.0, rng.gen_range(0.0..0.5), 0.0);
                    world.add(Box::new(MovingSphere::new(
                        center,
                        center2,
                        0.0,
                        1.0,
                        0.2,
                        sphere_material,
                    )));
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

fn two_spheres() -> HittableList {
    let mut world = HittableList::default();

    let checker = Box::new(Checker::from_solid_colors(
        Color::new(0.2, 0.3, 0.1),
        Color::new(0.9, 0.9, 0.9),
    ));
    world.add(Box::new(Sphere::new(
        Point::new(0.0, -10.0, 0.0),
        10.0,
        Box::new(Lambertian::from_texture(checker.clone())),
    )));
    world.add(Box::new(Sphere::new(
        Point::new(0.0, 10.0, 0.0),
        10.0,
        Box::new(Lambertian::from_texture(checker.clone())),
    )));

    world
}

fn two_perlin_spheres() -> HittableList {
    let mut world = HittableList::default();

    let pertext = Box::new(Noise::new(4.0));
    world.add(Box::new(Sphere::new(
        Point::new(0.0, -1000.0, 0.0),
        1000.0,
        Box::new(Lambertian::from_texture(pertext.clone())),
    )));
    world.add(Box::new(Sphere::new(
        Point::new(0.0, 2.0, 0.0),
        2.0,
        Box::new(Lambertian::from_texture(pertext)),
    )));

    world
}

fn simple_light() -> HittableList {
    let mut world = HittableList::default();

    let pertext = Box::new(Noise::new(4.0));
    world.add(Box::new(Sphere::new(
        Point::new(0.0, -1000.0, 0.0),
        1000.0,
        Box::new(Lambertian::from_texture(pertext.clone())),
    )));
    world.add(Box::new(Sphere::new(
        Point::new(0.0, 2.0, 0.0),
        2.0,
        Box::new(Lambertian::from_texture(pertext)),
    )));

    let difflight = Box::new(DiffuseLight::from_color(Color::new(4.0, 4.0, 4.0)));
    world.add(Box::new(Rectangle::new(
        difflight,
        3.0,
        5.0,
        1.0,
        3.0,
        -2.0,
        RectangleType::Xy,
    )));

    world
}

fn cornell_box() -> HittableList {
    let mut world = HittableList::default();

    let red = Box::new(Lambertian::new(Color::new(0.65, 0.05, 0.05)));
    let white = Box::new(Lambertian::new(Color::new(0.73, 0.73, 0.73)));
    let green = Box::new(Lambertian::new(Color::new(0.12, 0.45, 0.15)));
    let light = Box::new(DiffuseLight::from_color(Color::new(15.0, 15.0, 15.0)));

    world.add(Box::new(Rectangle::new(
        green,
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        RectangleType::Yz,
    )));
    world.add(Box::new(Rectangle::new(
        red,
        0.0,
        555.0,
        0.0,
        555.0,
        0.0,
        RectangleType::Yz,
    )));
    world.add(Box::new(Rectangle::new(
        light,
        213.0,
        343.0,
        227.0,
        332.0,
        554.0,
        RectangleType::Xz,
    )));
    world.add(Box::new(Rectangle::new(
        white.clone(),
        0.0,
        555.0,
        0.0,
        555.0,
        0.0,
        RectangleType::Xz,
    )));
    world.add(Box::new(Rectangle::new(
        white.clone(),
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        RectangleType::Xz,
    )));
    world.add(Box::new(Rectangle::new(
        white.clone(),
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        RectangleType::Xy,
    )));

    let box1 = Box::new(BoxObj::new(
        Point::new(0.0, 0.0, 0.0),
        Point::new(165.0, 330.0, 165.0),
        white.clone(),
    ));
    let box1 = Box::new(RotateY::new(box1, 15.0));
    let box1 = Box::new(Translate::new(box1, Vec3::new(265.0, 0.0, 295.0)));
    world.add(box1);
    let box2 = Box::new(BoxObj::new(
        Point::new(0.0, 0.0, 0.0),
        Point::new(165.0, 165.0, 165.0),
        white,
    ));
    let box2 = Box::new(RotateY::new(box2, -18.0));
    let box2 = Box::new(Translate::new(box2, Vec3::new(130.0, 0.0, 65.0)));
    world.add(box2);
    // world.add(Box::new(BoxObj::new(
    //     Point::new(130.0, 0.0, 65.0),
    //     Point::new(295.0, 165.0, 230.0),
    //     white.clone(),
    // )));
    // world.add(Box::new(BoxObj::new(
    //     Point::new(265.0, 0.0, 295.0),
    //     Point::new(430.0, 330.0, 460.0),
    //     white,
    // )));
    world
}

fn cornell_smoke() -> HittableList {
    let mut world = HittableList::default();

    let red = Box::new(Lambertian::new(Color::new(0.65, 0.05, 0.05)));
    let white = Box::new(Lambertian::new(Color::new(0.73, 0.73, 0.73)));
    let green = Box::new(Lambertian::new(Color::new(0.12, 0.45, 0.15)));
    let light = Box::new(DiffuseLight::from_color(Color::new(7.0, 7.0, 7.0)));

    world.add(Box::new(Rectangle::new(
        green,
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        RectangleType::Yz,
    )));
    world.add(Box::new(Rectangle::new(
        red,
        0.0,
        555.0,
        0.0,
        555.0,
        0.0,
        RectangleType::Yz,
    )));
    world.add(Box::new(Rectangle::new(
        light,
        113.0,
        443.0,
        127.0,
        432.0,
        554.0,
        RectangleType::Xz,
    )));
    world.add(Box::new(Rectangle::new(
        white.clone(),
        0.0,
        555.0,
        0.0,
        555.0,
        0.0,
        RectangleType::Xz,
    )));
    world.add(Box::new(Rectangle::new(
        white.clone(),
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        RectangleType::Xz,
    )));
    world.add(Box::new(Rectangle::new(
        white.clone(),
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        RectangleType::Xy,
    )));

    let box1 = Box::new(BoxObj::new(
        Point::new(0.0, 0.0, 0.0),
        Point::new(165.0, 330.0, 165.0),
        white.clone(),
    ));
    let box1 = Box::new(RotateY::new(box1, 15.0));
    let box1 = Box::new(Translate::new(box1, Vec3::new(265.0, 0.0, 295.0)));
    let box2 = Box::new(BoxObj::new(
        Point::new(0.0, 0.0, 0.0),
        Point::new(165.0, 165.0, 165.0),
        white.clone(),
    ));
    let box2 = Box::new(RotateY::new(box2, -18.0));
    let box2 = Box::new(Translate::new(box2, Vec3::new(130.0, 0.0, 65.0)));

    world.add(Box::new(ConstantMedium::from_color(
        box1,
        0.01,
        Color::new(0.0, 0.0, 0.0),
    )));
    world.add(Box::new(ConstantMedium::from_color(
        box2,
        0.01,
        Color::new(1.0, 1.0, 1.0),
    )));
    world
}

#[derive(Debug, Serialize, Deserialize)]
struct Config {
    /// Aspect ratio (width / height)
    /// Image width (in pixels)
    #[serde(default = "Config::default_image_width")]
    image_width: usize,
    /// Samples per pixel (number of rays sent per pixel)
    #[serde(default = "Config::default_samples_per_pixel")]
    samples_per_pixel: usize,
    /// Max number of bounces for a ray
    #[serde(default = "Config::default_max_depth")]
    max_depth: u32,
    camera: CameraConfig,
    scene: Scene,
}
impl Config {
    fn default_image_width() -> usize {
        400
    }

    fn default_samples_per_pixel() -> usize {
        500
    }

    fn default_max_depth() -> u32 {
        50
    }
}

/// Scene enum for canned or config
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum Scene {
    Random,
    TwoSpheres,
    TwoPerlinSpheres,
    SimpleLight,
    Cornell,
    CornellSmoke,
    List(HittableListConfig),
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    config: String,

    #[arg(short, long, default_value_t = String::from("./image.jpeg"))]
    output: String,
}

fn main() {
    let args = Args::parse();

    // Read in
    let config: Config = {
        let f = std::fs::File::open(args.config).expect("Could not open config file.");
        serde_yaml::from_reader(f).expect("Could not read values.")
    };

    let image_height: usize =
        (config.image_width as f64 / config.camera.aspect_ratio).round() as usize;

    // Create World
    // let world = random_scene();
    let world = match config.scene {
        Scene::Random => random_scene(),
        Scene::TwoSpheres => two_spheres(),
        Scene::TwoPerlinSpheres => two_perlin_spheres(),
        Scene::SimpleLight => simple_light(),
        Scene::Cornell => cornell_box(),
        Scene::CornellSmoke => cornell_smoke(),
        Scene::List(c) => HittableList::from_config(c),
    };
    // BvhNode for the win!
    let world = BvhNode::new(world, 0.0, 1.0);
    let protected_world = Arc::new(RwLock::new(world));
    // let background_color = Color::new(0.70, 0.8, 1.0);
    let background_color = Color::new(0.0, 0.0, 0.0);

    let cam = Camera::new(
        config.camera.look_from.into(),
        config.camera.look_at.into(),
        config.camera.v_up.into(),
        config.camera.vertical_fov_deg,
        config.camera.aspect_ratio,
        config.camera.aperture,
        config.camera.focus_distance,
        // TODO(mkagie) Add this to the config
        0.0,
        1.0,
    );

    // Random generator
    let mut rng = rand::thread_rng();

    // Create an output image
    let mut out_image = RgbImage::new(config.image_width as u32, image_height as u32);

    // Create the threadpool
    let n_workers = 20; // Number of cores I have
    let pool = ThreadPool::new(n_workers);

    // Render
    let bar = ProgressBar::new((config.image_width * image_height) as u64);
    bar.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise} --- {eta_precise} --- {pos}/{len}] {wide_bar}",
        )
        .unwrap(),
    );
    for j in (0..=image_height - 1).rev() {
        for i in 0..config.image_width {
            let (worker_tx, worker_rx) = channel();
            for _ in 0..config.samples_per_pixel {
                let tx = worker_tx.clone();
                let u = (i as f64 + rng.gen::<f64>()) / (config.image_width - 1) as f64;
                let v = (j as f64 + rng.gen::<f64>()) / (image_height - 1) as f64;
                let ray = cam.get_ray(u, v);
                let protected_world = Arc::clone(&protected_world);
                pool.execute(move || {
                    tx.send(ray.get_color(
                        &background_color,
                        &*protected_world.read().unwrap(),
                        config.max_depth,
                    ))
                    .unwrap();
                });
            }
            let pixel_color = worker_rx
                .iter()
                .take(config.samples_per_pixel)
                .reduce(|a, b| a + b)
                .unwrap();

            bar.inc(1);
            out_image.put_pixel(
                i as u32,
                // Invert the height
                (image_height - 1 - j) as u32,
                utils::get_pixel(&pixel_color, config.samples_per_pixel),
            );
        }
    }
    bar.finish();
    out_image.save(args.output).unwrap();
}
