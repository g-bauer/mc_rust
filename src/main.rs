extern crate rand;
use std::io::{self, BufReader, BufWriter};
use std::io::prelude::*;
use std::fs::File;
use rand::{Rng, SeedableRng, XorShiftRng};

use my_vec::Vector3D;
mod my_vec;
use mymod::{energy, particle_energy, translate};
mod mymod;

fn read_config(filename: String) -> Vec<Vector3D> {
    // let f = try!(File::open(filename));
    let f = BufReader::new(File::open(filename).unwrap());
    let mut particles: Vec<Vector3D> = vec![];
    for line in f.lines() {
        let mut dummy: Vec<f64> = vec![];
        for elem in line.unwrap().split_whitespace() {
            dummy.push(elem.parse().unwrap());
        }
        if dummy.len() == 4 {
            particles.push(
                Vector3D::new(
                    dummy[1], dummy[2], dummy[3])
                )
        }
    }
    particles
}


#[allow(dead_code)]
fn set_lattice(n: &u32, l: &f64) -> Vec<Vector3D> {
    let del = ((*n as f64).powf(1./3.)+1f64) as u32;
    let delt = (del as f64)/l;
    let mut particles: Vec<Vector3D> = vec![];
    for i in 0..del {
        for j in 0..del {
            for k in 0..del {
                particles.push(Vector3D::new(
                    i as f64 * delt,
                    j as f64 * delt,
                    k as f64 * delt));
                if particles.len() >= *n as usize {
                    return particles;
                }
            }
        }
    }
    particles
}

fn main() {
    // rng
    let seed = [1u32, 2u32, 3u32, 4u32];
    let mut rng: XorShiftRng = XorShiftRng::from_seed(seed);
    // side length of simulation box
    let l = 10f64;  // the boxlength
    let rc2 = 9f64; // cut off in multiples of sigma [==1]
    // temperature
    let temperature = 1f64;
    let beta = 1f64/temperature;
    // let particles = set_lattice(&npart, &l); // set to lattice
    let mut particles = read_config("nist_lj.dat".to_owned()); // read from file
    let npart = particles.len(); // number of particles
    // numerics: #moves = npart
    let moves = npart;
    let cycles = 1000;
    let density = npart as f64/l.powi(3); // density

    // print some stuff to stdout
    println!("The system contains {} particles.", npart);
    println!("Density: {} / 1/sig^3", density);
    println!("Volume : {} / sig^3", l.powi(3));

    // prepare output file
    // since we write in every loop -> BufWriter
    let mut buffer = BufWriter::new(File::create("output.dat").unwrap());
    // compute initial system energy
    let mut ener = energy(&particles, &l, &rc2);
    println!("Initial system energy: {}", ener);
    // start simulation
    for i in 0..cycles {
        for j in 0..moves {
            ener += translate(&mut particles, &l, &rc2, &beta, &mut rng);
        }
        if i%10 == 0 {
            write!(buffer, "{:.4}\n", ener).unwrap();
        }
    }
    println!("Final energy: {}", ener);
    println!("Delta (system - running): {:e}", ener-energy(&particles, &l, &rc2));
}
