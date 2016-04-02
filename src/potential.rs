use my_vec::Vector3D;
use mymod::squared_distance;

pub trait PairPotential {
    fn energy(&self, r2: f64) -> f64;
    // fn virial(&self, r2: &f64) -> f64;
}

#[derive(Debug)]
pub struct LennardJones {
    sigma: f64,
    epsilon: f64,
}

impl LennardJones {
    pub fn new(sigma: f64, epsilon: f64) -> LennardJones {
        LennardJones{sigma: sigma, epsilon: epsilon}
    }
}

impl PairPotential for LennardJones {
    fn energy(&self, r2: f64) -> f64 {
        let inv2 = self.sigma*self.sigma/r2;
        let inv6 = inv2.powi(3);
        4f64*self.epsilon*(inv6*inv6 - inv6)
    }
}

#[derive(Debug)]
pub struct Mie {
    sigma: f64,
    epsilon: f64,
    n: i32,
    m: i32,
    prefac: f64,
}

impl Mie {
    pub fn new(sigma: f64, epsilon: f64, n: i32, m: i32) -> Mie {

        let nf = n as f64;
        let mf = m as f64;
        let p = nf/(nf-mf)*(nf/mf).powf(mf/(nf-mf));
        Mie{sigma: sigma, epsilon: epsilon, n: n, m: m, prefac: p}
    }
}

impl PairPotential for Mie {
    fn energy(&self, r2: f64) -> f64 {
        let inv = self.sigma/r2.sqrt();
        self.prefac*self.epsilon*(inv.powi(self.n) - inv.powi(self.m))
    }
}

/// Compute system energy using a generic pair potential.
pub fn system_energy<T: PairPotential>(potential: &T, particles: &Vec<Vector3D>,
                                       l: &f64, rc2: &f64) -> f64 {
    let mut energy = 0f64;
    for i in 0..particles.len()-1 {
        let particle_i = particles[i]; // look up of particle i once for all j
        for j in i+1..particles.len() {
        //for j in i+1..npart {
            // compute squared distance including nearest image
            let d2 = squared_distance(&particle_i, &particles[j], &l);
            if d2 <= *rc2 {
                // only compute for distance greater than cut off
                energy += potential.energy(d2);
            } else {
                // else, cycle into nex j iteration
                continue;
            }
        }
    }
    energy
}
