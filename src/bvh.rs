//! Axis-Aligned Bounding Boxes and Bounding Volume Hierarchies
//!
//! For optimizing computations

use std::cmp::Ordering;

use crate::objects::{Hittable, HittableList, HittableObj};
use crate::{Point, Ray};
use rand::Rng;

/// Axis-Aligned Bounding Box
///
/// A data structure to bound many objects to speed up computations
#[derive(Debug, Clone)]
pub struct Aabb {
    pub min: Point,
    pub max: Point,
}
impl Aabb {
    pub fn new(min: Point, max: Point) -> Self {
        Self { min, max }
    }

    /// Whether the box is hit by a ray between the time range
    ///
    /// Original implementation from the paper
    pub fn hit_orig(&self, r: &Ray, mut t_min: f64, mut t_max: f64) -> bool {
        for a in 0..3 {
            let t0 =
                ((self.min[a] - r.orig[a]) / r.dir[a]).min((self.max[a] - r.orig[a]) / r.dir[a]);
            let t1 =
                ((self.min[a] - r.orig[a]) / r.dir[a]).max((self.max[a] - r.orig[a]) / r.dir[a]);

            t_min = t0.max(t_min);
            t_max = t1.min(t_max);
            if t_max <= t_min {
                return false;
            }
        }
        true
    }

    /// Whether the box is hit by a ray between the time range
    ///
    /// Improved implementation
    pub fn hit(&self, r: &Ray, mut t_min: f64, mut t_max: f64) -> bool {
        for a in 0..3 {
            let inv_d = 1.0 / r.dir[a];
            let mut t0 = (self.min[a] - r.orig[a]) * inv_d;
            let mut t1 = (self.max[a] - r.orig[a]) * inv_d;
            if inv_d < 0.0 {
                std::mem::swap(&mut t0, &mut t1);
            }
            t_min = if t0 > t_min { t0 } else { t_min };
            t_max = if t1 < t_max { t1 } else { t_max };
            if t_max <= t_min {
                return false;
            }
        }
        true
    }

    /// Compute the surrounding AABB between this and another
    pub fn surrounding_box(&self, other: &Aabb) -> Aabb {
        let small = Point::new(
            self.min[0].min(other.min[0]),
            self.min[1].min(other.min[1]),
            self.min[2].min(other.min[2]),
        );
        let big = Point::new(
            self.max[0].max(other.max[0]),
            self.max[1].max(other.max[1]),
            self.max[2].max(other.max[2]),
        );
        Aabb::new(small, big)
    }
}

/// Bounding Volume Hierarchy
///
/// Tree structure
pub struct BvhNode {
    left: HittableObj,
    // Going to make right an Option, so we don't have to clone left
    right: Option<HittableObj>,
    bbox: Aabb,
}
impl BvhNode {
    pub fn new(list: HittableList, time0: f64, time1: f64) -> Self {
        Self::split_tree(list, time0, time1)
    }

    /// Split the tree
    ///
    /// Randomly choose an axis, sort the primitives, put half in each subtree
    fn split_tree(mut objects: HittableList, time0: f64, time1: f64) -> Self {
        let mut rng = rand::thread_rng();
        let axis = rng.gen_range(0..=2);

        let (left, right) = if objects.len() == 1 {
            // There is only one object, both left and right must be the same
            let left = objects.0.pop().unwrap();
            (left, None)
        } else if objects.len() == 2 {
            // There are two objects, order and set
            let obj_start_p1 = objects.0.pop().unwrap();
            let obj_start = objects.0.pop().unwrap();
            if matches!(box_compare(&obj_start, &obj_start_p1, axis), Ordering::Less) {
                (obj_start, Some(obj_start_p1))
            } else {
                (obj_start_p1, Some(obj_start))
            }
        } else {
            // There are many objects, order and split and call again
            objects.sort_by(|a, b| box_compare(a, b, axis));
            let mid = objects.len() / 2;
            let right_list = HittableList(objects.0.split_off(mid));

            let left: HittableObj = Box::new(BvhNode::split_tree(objects, time0, time1));
            let right: HittableObj = Box::new(BvhNode::split_tree(right_list, time0, time1));
            (left, Some(right))
        };

        // Now, compute the bbox, which is the surrounding box of both of these
        match right {
            Some(right) => {
                // There is a right
                match (
                    left.try_bounding_box(time0, time1),
                    right.try_bounding_box(time0, time1),
                ) {
                    (Some(box_left), Some(box_right)) => {
                        let bbox = box_left.surrounding_box(&box_right);
                        Self {
                            left,
                            right: Some(right),
                            bbox,
                        }
                    }
                    _ => panic!("No bounding box in bvh constructor"),
                }
            }
            None => {
                // There is no right, just do left
                let bbox = left.try_bounding_box(time0, time1);
                match bbox {
                    Some(bbox) => Self {
                        left,
                        right: None,
                        bbox,
                    },
                    None => panic!("No bounding box in bvh constructor"),
                }
            }
        }
    }
}
impl Hittable for BvhNode {
    fn try_hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<crate::objects::HitRecord> {
        if !self.bbox.hit(ray, t_min, t_max) {
            return None;
        }
        if let Some(rec) = self.left.try_hit(ray, t_min, t_max) {
            // We have hit left
            match &self.right {
                Some(right) => {
                    if let Some(rec) = right.try_hit(ray, t_min, rec.t) {
                        // We hit right and overwrote rec
                        Some(rec)
                    } else {
                        // We did not hit right, return left rec
                        Some(rec)
                    }
                }
                None => {
                    // There is no right, return left
                    Some(rec)
                }
            }
        } else {
            // We did not hit left
            match &self.right {
                Some(right) => right.try_hit(ray, t_min, t_max),
                None => None,
            }
        }
    }

    fn try_bounding_box(&self, _time0: f64, _time1: f64) -> Option<Aabb> {
        Some(self.bbox.clone())
    }
}

/// Compare boxes from two HittableObjs
fn box_compare(a: &HittableObj, b: &HittableObj, axis: usize) -> Ordering {
    let box_a = a.try_bounding_box(0.0, 0.0);
    let box_b = b.try_bounding_box(0.0, 0.0);

    match (box_a, box_b) {
        (Some(box_a), Some(box_b)) => box_a.min[axis].partial_cmp(&box_b.min[axis]).unwrap(),
        _ => panic!("No bounding box in bvh constructor"),
    }
}
