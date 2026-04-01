//! Validation tests — cross-checking rasayan models against published experimental data.
//!
//! Each test cites the source of the expected value. These tests ensure our
//! implementations produce physiologically/biochemically accurate results.

use rasayan::enzyme;
use rasayan::membrane;
use rasayan::protein;

// ==========================================================================
// Michaelis-Menten — validate at Km gives Vmax/2
// Source: Lehninger Principles of Biochemistry, 8th ed., Ch. 6
// ==========================================================================

#[test]
fn validate_mm_at_km_gives_half_vmax() {
    // Universal property: v(Km) = Vmax/2
    let vmax = 100.0;
    let km = 5.0;
    let rate = enzyme::michaelis_menten(km, vmax, km);
    assert!(
        (rate - vmax / 2.0).abs() < 0.001,
        "At [S]=Km, rate should be Vmax/2, got {rate}"
    );
}

#[test]
fn validate_mm_saturation() {
    // At [S] >> Km, rate → Vmax
    let vmax = 100.0;
    let km = 1.0;
    let rate = enzyme::michaelis_menten(1000.0, vmax, km);
    assert!(
        (rate - vmax).abs() < 0.2,
        "At [S]>>Km, rate should approach Vmax, got {rate}"
    );
}

// ==========================================================================
// Hill equation — hemoglobin O2 binding
// Source: Severinghaus (1979), Voet & Voet Biochemistry 4th ed., Ch. 7
// Hemoglobin: n ≈ 2.8, P50 ≈ 26.8 mmHg
// At pO2=100 mmHg, saturation ≈ 97.5%
// ==========================================================================

#[test]
fn validate_hill_hemoglobin_o2_binding() {
    let n = 2.8;
    let p50 = 26.8; // mmHg (K0.5)
    let vmax = 1.0; // fractional saturation max = 1.0
    let po2 = 100.0; // mmHg, arterial

    let saturation = enzyme::hill_equation(po2, vmax, p50, n);
    assert!(
        (saturation - 0.975).abs() < 0.02,
        "Hemoglobin at pO2=100 should be ~97.5% saturated, got {saturation}"
    );
}

#[test]
fn validate_hill_hemoglobin_at_p50() {
    let n = 2.8;
    let p50 = 26.8;
    let saturation = enzyme::hill_equation(p50, 1.0, p50, n);
    assert!(
        (saturation - 0.5).abs() < 0.001,
        "At P50, saturation should be 50%, got {saturation}"
    );
}

#[test]
fn validate_hill_hemoglobin_venous() {
    // Venous pO2 ≈ 40 mmHg → saturation ≈ 75%
    let n = 2.8;
    let p50 = 26.8;
    let saturation = enzyme::hill_equation(40.0, 1.0, p50, n);
    assert!(
        (saturation - 0.75).abs() < 0.05,
        "Venous pO2=40: saturation should be ~75%, got {saturation}"
    );
}

// ==========================================================================
// Competitive inhibition — succinate dehydrogenase + malonate
// Source: Lehninger 8th ed., Ch. 6; Quastel & Wooldridge (1928)
// Km(succinate) = 0.4 mM, Ki(malonate) = 0.05 mM
// Apparent Km at [I]=1mM: Km_app = Km(1 + [I]/Ki) = 0.4 * 21 = 8.4 mM
// ==========================================================================

#[test]
fn validate_competitive_inhibition_malonate() {
    let vmax = 10.0;
    let km = 0.4; // mM
    let ki = 0.05; // mM
    let inhibitor = 1.0; // mM

    // At [S] = Km_app, rate should be Vmax/2 with inhibitor
    // Km_app = Km * (1 + [I]/Ki) = 0.4 * (1 + 1.0/0.05) = 0.4 * 21 = 8.4
    let km_app = km * (1.0 + inhibitor / ki);
    assert!(
        (km_app - 8.4_f64).abs() < 0.01,
        "Apparent Km should be 8.4, got {km_app}"
    );

    let rate_at_km_app = enzyme::competitive_inhibition(km_app, inhibitor, vmax, km, ki);
    assert!(
        (rate_at_km_app - vmax / 2.0).abs() < 0.1,
        "At [S]=Km_app, rate should be Vmax/2, got {rate_at_km_app}"
    );
}

#[test]
fn validate_competitive_inhibition_zero_inhibitor_equals_mm() {
    let vmax = 10.0;
    let km = 0.4;
    let ki = 0.05;
    let s = 2.0;

    let rate_mm = enzyme::michaelis_menten(s, vmax, km);
    let rate_ci = enzyme::competitive_inhibition(s, 0.0, vmax, km, ki);
    assert!(
        (rate_mm - rate_ci).abs() < 0.001,
        "Zero inhibitor: CI should equal MM, got {rate_ci} vs {rate_mm}"
    );
}

// ==========================================================================
// Arrhenius — temperature dependence
// Source: Stryer Biochemistry 9th ed.; DeLuca & McElroy (1974)
// Ea ~ 50 kJ/mol, rate doubles per 10°C (Q10 ≈ 2) near 25°C
// ==========================================================================

#[test]
fn validate_arrhenius_temperature_sensitivity() {
    let ea = 50_000.0; // J/mol
    let a = 1e10; // pre-exponential (arbitrary scale)
    let t1 = 298.0; // 25°C
    let t2 = 308.0; // 35°C

    let k1 = enzyme::arrhenius(a, ea, t1);
    let k2 = enzyme::arrhenius(a, ea, t2);
    let ratio = k2 / k1;

    // For Ea=50 kJ/mol, 10°C increase should give ~1.9-2.1× rate increase
    assert!(
        (1.8..=2.2).contains(&ratio),
        "10°C increase with Ea=50kJ should ~double rate, got {ratio:.2}×"
    );
}

// ==========================================================================
// Nernst equation — potassium equilibrium potential
// Source: Hodgkin & Huxley (1952), Kandel Principles of Neural Science 6th ed.
// Mammalian neuron (37°C = 310K): [K+]in=140mM, [K+]out=5mM
// E_K = (RT/zF)·ln(5/140) ≈ −89 mV
// ==========================================================================

#[test]
fn validate_nernst_potassium_mammalian() {
    let temp = 310.0; // K (37°C)
    let z = 1;
    let k_out = 5.0; // mM
    let k_in = 140.0; // mM

    let e_k_mv = membrane::nernst(temp, z, k_out, k_in);

    assert!(
        (e_k_mv - (-89.0)).abs() < 3.0,
        "E_K at 37°C should be ~−89 mV, got {e_k_mv:.1} mV"
    );
}

#[test]
fn validate_nernst_sodium_mammalian() {
    // [Na+]in=12mM, [Na+]out=145mM → E_Na ≈ +67 mV
    // Source: Kandel 6th ed., Table 7-1
    let e_na_mv = membrane::nernst(310.0, 1, 145.0, 12.0);

    assert!(
        (e_na_mv - 67.0).abs() < 5.0,
        "E_Na at 37°C should be ~+67 mV, got {e_na_mv:.1} mV"
    );
}

// ==========================================================================
// Goldman equation — resting membrane potential
// Source: Hodgkin & Katz (1949), Hille Ion Channels 3rd ed.
// Default mammalian neuron: V_rest ≈ −60 to −70 mV
// ==========================================================================

#[test]
fn validate_goldman_resting_potential() {
    let ions = membrane::IonicState::default();
    let perm = membrane::MembranePermeability::default();
    let v_mv = membrane::goldman(310.0, &ions, &perm);

    assert!(
        (-80.0..=-50.0).contains(&v_mv),
        "Resting potential should be −50 to −80 mV, got {v_mv:.1} mV"
    );
}

// ==========================================================================
// Isoelectric point — published protein pI values
// Source: ExPASy Compute pI (Swiss-Prot)
// ==========================================================================

#[test]
fn validate_pi_glycine_free_amino_acid() {
    // Glycine: pI = 5.97 (Lehninger)
    let pi = protein::isoelectric_point("G").unwrap();
    assert!(
        (pi - 5.97).abs() < 0.2,
        "Glycine pI should be ~5.97, got {pi}"
    );
}

#[test]
fn validate_pi_basic_amino_acid_arginine() {
    // Arginine: pI ≈ (9.69 + 12.48) / 2 = 11.09 (Lehninger pKa values)
    let pi = protein::isoelectric_point("R").unwrap();
    assert!(
        (pi - 11.09).abs() < 0.15,
        "Arginine pI should be ~11.09, got {pi}"
    );
}

#[test]
fn validate_pi_acidic_amino_acid_glutamate() {
    // Glutamate: pI ≈ (2.34 + 4.25) / 2 = 3.30 (Lehninger pKa values)
    let pi = protein::isoelectric_point("E").unwrap();
    assert!(
        (pi - 3.30).abs() < 0.15,
        "Glutamate pI should be ~3.30, got {pi}"
    );
}

// ==========================================================================
// Extinction coefficient — known protein chromophore contributions
// Source: Pace et al., Protein Science 4:2411 (1995)
// εTrp = 5500, εTyr = 1490, εCystine = 125 M⁻¹cm⁻¹
// ==========================================================================

#[test]
fn validate_extinction_coefficient_pace_values() {
    // A single Trp: ε = 5500
    let ec = protein::extinction_coefficient("W").unwrap();
    assert!(
        (ec.reduced - 5500.0).abs() < f64::EPSILON,
        "εTrp should be 5500, got {}",
        ec.reduced
    );

    // Two Tyr: ε = 2 × 1490 = 2980
    let ec2 = protein::extinction_coefficient("YY").unwrap();
    assert!(
        (ec2.reduced - 2980.0).abs() < f64::EPSILON,
        "2×εTyr should be 2980, got {}",
        ec2.reduced
    );
}

// ==========================================================================
// Enzyme database — validate against published Km/kcat
// Source: Lehninger 8th ed., BRENDA enzyme database
// ==========================================================================

#[test]
fn validate_enzyme_db_carbonic_anhydrase() {
    // Carbonic anhydrase II: kcat ~ 1,000,000 s⁻¹, Km ~ 12-26 mM
    let enz = enzyme::lookup_enzyme("Carbonic anhydrase").unwrap();
    assert!(
        enz.kcat > 500_000.0 && enz.kcat < 2_000_000.0,
        "CA kcat should be ~10^6 s⁻¹, got {}",
        enz.kcat
    );
}

#[test]
fn validate_enzyme_db_catalase() {
    // Catalase: kcat ~ 40,000,000 s⁻¹ (fastest known)
    let enz = enzyme::lookup_enzyme("Catalase").unwrap();
    assert!(
        enz.kcat > 10_000_000.0,
        "Catalase kcat should be >10^7 s⁻¹, got {}",
        enz.kcat
    );
}

#[test]
fn validate_enzyme_db_hexokinase() {
    // Hexokinase: Km(glucose) ~ 0.1 mM = 0.0001 M (high affinity)
    // Enzyme DB stores values in M
    let enz = enzyme::lookup_enzyme("Hexokinase").unwrap();
    assert!(
        enz.km < 0.001 && enz.km > 0.00001,
        "Hexokinase Km should be ~0.0001 M (0.1 mM), got {}",
        enz.km
    );
}

// ==========================================================================
// BLOSUM62 — validate known matrix entries
// Source: Henikoff & Henikoff, PNAS 89:10915 (1992)
// ==========================================================================

#[test]
fn validate_blosum62_known_entries() {
    use rasayan::alignment::{self, Matrix};

    // W-W = 11 (highest self-score)
    assert_eq!(
        alignment::substitution_score('W', 'W', Matrix::Blosum62).unwrap(),
        11
    );
    // C-C = 9 (cysteine, high self-score due to disulfide conservation)
    assert_eq!(
        alignment::substitution_score('C', 'C', Matrix::Blosum62).unwrap(),
        9
    );
    // D-E = 2 (similar acidic residues, positive score)
    assert_eq!(
        alignment::substitution_score('D', 'E', Matrix::Blosum62).unwrap(),
        2
    );
    // W-G = -2 (dissimilar, negative)
    assert!(alignment::substitution_score('W', 'G', Matrix::Blosum62).unwrap() < 0);
}
