//! Basic usage of rasayan — enzyme kinetics and membrane potential.

use rasayan::enzyme::{self, EnzymeParams};
use rasayan::membrane::{self, IonicState};
use rasayan::metabolism::MetabolicState;
use rasayan::signal;

fn main() {
    // --- Enzyme kinetics: Michaelis-Menten rate curve ---
    let params = EnzymeParams {
        vmax: 10.0,  // M/s
        km: 1.0,     // M (Michaelis constant)
        hill_n: 1.0, // standard (no cooperativity)
        kcat: 100.0,
    };

    println!("=== Enzyme Kinetics ===");
    println!(
        "Catalytic efficiency: {:.0} M^-1 s^-1",
        params.catalytic_efficiency()
    );
    println!();

    println!("[S] (M)    Rate (M/s)   % Vmax");
    println!("-------    ----------   ------");
    for &s in &[0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0] {
        let rate = params.rate(s);
        println!("{s:6.2}     {rate:8.4}     {pct:5.1}%", pct = rate / params.vmax * 100.0);
    }
    println!();

    // --- Hill equation: cooperative binding ---
    println!("=== Cooperative Binding (Hill) ===");
    println!("n=1 (no coop)  vs  n=4 (strong coop) at [S]=Km:");
    let r1 = enzyme::hill_equation(1.0, 10.0, 1.0, 1.0);
    let r4 = enzyme::hill_equation(1.0, 10.0, 1.0, 4.0);
    println!("  n=1: {r1:.2} M/s");
    println!("  n=4: {r4:.2} M/s");
    println!();

    // --- Membrane potential ---
    println!("=== Membrane Potential ===");
    let ions = IonicState::default();
    let vm = ions.resting_potential();
    println!("Resting membrane potential: {vm:.1} mV");

    let e_k = membrane::nernst(310.0, 1, ions.k_out, ions.k_in);
    let e_na = membrane::nernst(310.0, 1, ions.na_out, ions.na_in);
    println!("Nernst K+:  {e_k:.1} mV");
    println!("Nernst Na+: {e_na:.1} mV");
    println!();

    // --- Metabolic state ---
    println!("=== Metabolic State ===");
    let mut met = MetabolicState::default();
    println!("Energy charge: {:.3}", met.energy_charge());
    println!("Aerobic ATP yield: {:.0}", met.aerobic_atp_yield());
    met.consume_atp(3.0);
    println!("After consuming 3 ATP: ATP={:.1}, ADP={:.1}", met.atp, met.adp);
    println!();

    // --- Dose-response ---
    println!("=== Dose-Response Curve ===");
    println!("[L]/EC50    Response");
    println!("--------    --------");
    for &ratio in &[0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0] {
        let r = signal::dose_response(ratio, 1.0, 1.0, 1.0);
        println!("{ratio:8.2}    {r:8.4}");
    }
}
