use std::ops::{Add, Sub, Mul, Index};
use std::fmt::Display;

#[derive(Copy, Clone, Debug)]
pub struct Vector3D(pub [f64; 3]);

impl Vector3D {
    /// Create a new Vector3D with components `x`, `y`, `z`
    pub fn new(x: f64, y: f64, z: f64) -> Vector3D {
        Vector3D([x, y, z])
    }
    // /// Return the squared euclidean norm of a Vector3D
    // #[inline] pub fn norm2(&self) -> f64 {
    //     self * self
    // }
    // /// Return the euclidean norm of a Vector3D
    // #[inline] pub fn norm(&self) -> f64 {
    //     f64::sqrt(self.norm2())
    // }
    // /// Normalize a Vector3D
    // #[inline] pub fn normalized(&self) -> Vector3D {
    //     self / self.norm()
    // }
}

impl Index<usize> for Vector3D {
    type Output = f64;

    fn index(& self, _index: usize) -> & f64 {
        println!("Indexing!");
        &self.0[_index]
    }
}

impl Add for Vector3D {
    type Output = Vector3D;
    fn add(self, other: Vector3D) -> Vector3D {
        Vector3D::new(self.0[0] + other.0[0],
                      self.0[1] + other.0[1],
                      self.0[2] + other.0[2])
    }
}
