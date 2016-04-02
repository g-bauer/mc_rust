use my_vec::Vector3D;
use rand::{Rng, SeedableRng, XorShiftRng};
use rand::distributions::{IndependentSample, Range};

fn lennard_jones(r2: &f64) -> f64 {
    let inv2 = 1.0/r2;
    let inv6 = inv2*inv2*inv2;
    4.0*1.0*(inv6*inv6 - inv6)
}

/// Compute squared distance between two positions including nearest image convention.
pub fn squared_distance(r1: &Vector3D, r2: &Vector3D, l: &f64) -> f64 {
    let dx = (r1.0[0] - r2.0[0]) - ((r1.0[0] - r2.0[0])/l).round() *l;
    let dy = (r1.0[1] - r2.0[1]) - ((r1.0[1] - r2.0[1])/l).round() *l;
    let dz = (r1.0[2] - r2.0[2]) - ((r1.0[2] - r2.0[2])/l).round() *l;
    dx*dx + dy*dy + dz*dz
}

/// Compute system energy using the Lennard-Jones potential.
pub fn energy(particles: &Vec<Vector3D>, l: &f64, rc2: &f64) -> f64 {
    let mut energy = 0f64;
    for i in 0..particles.len()-1 {
        let particle_i = particles[i]; // look up of particle i once for all j
        for j in i+1..particles.len() {
        //for j in i+1..npart {
            // compute squared distance including nearest image
            let d2 = squared_distance(&particle_i, &particles[j], &l);
            if d2 <= *rc2 {
                // only compute for distance greater than cut off
                energy += lennard_jones(&d2);
            } else {
                // else, cycle into nex j iteration
                continue;
            }
        }
    }
    energy
}

/// Compute system energy using the Lennard-Jones potential.
pub fn particle_energy(particles: &Vec<Vector3D>, idx: &usize, l: &f64, rc2: &f64) -> f64 {
    let mut energy = 0f64;
    let particle_j = particles[*idx];

    // range does not contain idx: 0..idx = 0, 1, ..., idx-1
    for i in 0..*idx {
        // compute squared distance including nearest image
        let d2 = squared_distance(&particle_j, &particles[i], &l);
        if d2 <= *rc2 {
            // only compute for distance greater than cut off
            energy += lennard_jones(&d2);
        } else {
            // else, cycle into nex j iteration
            continue;
        }
    }
    for i in *idx+1..particles.len() {
        // compute squared distance including nearest image
        let d2 = squared_distance(&particle_j, &particles[i], &l);
        if d2 <= *rc2 {
            // only compute for distance greater than cut off
            energy += lennard_jones(&d2);
        } else {
            // else, cycle into nex j iteration
            continue;
        }
    }
    energy
}

pub fn translate(particles: &mut Vec<Vector3D>, l: &f64, rc2: &f64, beta: &f64, rng: &mut XorShiftRng) -> f64 {
    let delta = 0f64;
    let delta_range = Range::new(-0.3f64, 0.3f64);
    // randomly pick a particle
    let select_idx = rng.gen_range(0, particles.len());
    let old_particle = particles[select_idx].clone();

    // compute old energy
    let old_energy = particle_energy(&particles, &select_idx, &l, &rc2);

    let displace = Vector3D::new(
        delta_range.ind_sample(rng),
        delta_range.ind_sample(rng),
        delta_range.ind_sample(rng)
    );
    particles[select_idx] = particles[select_idx] + displace;

    // compute new energy
    let new_energy = particle_energy(&particles, &select_idx, &l, &rc2);
    // acceptance
    let mut delta_energy = new_energy - old_energy;
    let accepted = beta*delta_energy <= 0.0 || rng.next_f64() < f64::exp(-beta*delta_energy);
    if !accepted {
        delta_energy = 0f64;
        particles[select_idx] = old_particle.clone();
    }
    delta_energy
}
