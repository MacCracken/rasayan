//! Integration tests for rasayan — cross-module behavior and serde roundtrips.

use rasayan::energy::BioenergyState;
use rasayan::enzyme::{self, EnzymeParams};
use rasayan::membrane::{self, IonicState};
use rasayan::metabolism::MetabolicState;
use rasayan::signal::{self, SecondMessenger};

// --- Serde roundtrip tests ---

#[test]
fn test_metabolic_state_serde_roundtrip() {
    let m = MetabolicState::default();
    let json = serde_json::to_string(&m).unwrap();
    let m2: MetabolicState = serde_json::from_str(&json).unwrap();
    assert!((m2.atp - m.atp).abs() < f64::EPSILON);
    assert!((m2.adp - m.adp).abs() < f64::EPSILON);
    assert!((m2.glucose - m.glucose).abs() < f64::EPSILON);
    assert!((m2.oxygen - m.oxygen).abs() < f64::EPSILON);
    assert!((m2.lactate - m.lactate).abs() < f64::EPSILON);
    assert!((m2.nad_ratio - m.nad_ratio).abs() < f64::EPSILON);
}

#[test]
fn test_enzyme_params_serde_roundtrip() {
    let e = EnzymeParams {
        vmax: 12.5,
        km: 0.003,
        hill_n: 2.8,
        kcat: 450.0,
    };
    let json = serde_json::to_string(&e).unwrap();
    let e2: EnzymeParams = serde_json::from_str(&json).unwrap();
    assert!((e2.vmax - e.vmax).abs() < f64::EPSILON);
    assert!((e2.km - e.km).abs() < f64::EPSILON);
    assert!((e2.hill_n - e.hill_n).abs() < f64::EPSILON);
    assert!((e2.kcat - e.kcat).abs() < f64::EPSILON);
}

#[test]
fn test_second_messenger_serde_roundtrip() {
    let mut sm = SecondMessenger::default();
    sm.activate_gs(0.7);
    let json = serde_json::to_string(&sm).unwrap();
    let sm2: SecondMessenger = serde_json::from_str(&json).unwrap();
    assert!((sm2.camp - sm.camp).abs() < f64::EPSILON);
    assert!((sm2.calcium - sm.calcium).abs() < f64::EPSILON);
    assert!((sm2.ip3 - sm.ip3).abs() < f64::EPSILON);
}

#[test]
fn test_ionic_state_serde_roundtrip() {
    let ions = IonicState::default();
    let json = serde_json::to_string(&ions).unwrap();
    let ions2: IonicState = serde_json::from_str(&json).unwrap();
    assert!((ions2.na_in - ions.na_in).abs() < f64::EPSILON);
    assert!((ions2.na_out - ions.na_out).abs() < f64::EPSILON);
    assert!((ions2.k_in - ions.k_in).abs() < f64::EPSILON);
    assert!((ions2.k_out - ions.k_out).abs() < f64::EPSILON);
    assert!((ions2.cl_in - ions.cl_in).abs() < f64::EPSILON);
    assert!((ions2.cl_out - ions.cl_out).abs() < f64::EPSILON);
}

#[test]
fn test_bioenergy_state_serde_roundtrip() {
    let b = BioenergyState::default();
    let json = serde_json::to_string(&b).unwrap();
    let b2: BioenergyState = serde_json::from_str(&json).unwrap();
    assert!((b2.phosphocreatine - b.phosphocreatine).abs() < f64::EPSILON);
    assert!((b2.glycogen - b.glycogen).abs() < f64::EPSILON);
    assert!((b2.met - b.met).abs() < f64::EPSILON);
    assert!((b2.anaerobic_threshold - b.anaerobic_threshold).abs() < f64::EPSILON);
}

// --- Error display tests ---

#[test]
fn test_error_display_negative_concentration() {
    let err = rasayan::RasayanError::NegativeConcentration {
        name: "ATP".into(),
        value: -1.5,
    };
    let msg = err.to_string();
    assert!(msg.contains("ATP"));
    assert!(msg.contains("-1.5"));
    assert!(msg.contains("concentration"));
}

#[test]
fn test_error_display_invalid_parameter() {
    let err = rasayan::RasayanError::InvalidParameter {
        name: "Km".into(),
        value: -0.01,
        reason: "must be positive".into(),
    };
    let msg = err.to_string();
    assert!(msg.contains("Km"));
    assert!(msg.contains("-0.01"));
    assert!(msg.contains("must be positive"));
}

#[test]
fn test_error_display_unknown_amino_acid() {
    let err = rasayan::RasayanError::UnknownAminoAcid('X');
    assert!(err.to_string().contains('X'));
}

#[test]
fn test_error_display_unknown_pathway() {
    let err = rasayan::RasayanError::UnknownPathway("pentose phosphate".into());
    assert!(err.to_string().contains("pentose phosphate"));
}

// --- Function tests ---

#[test]
fn test_michaelis_menten_half_vmax_at_km() {
    // At [S] = Km, rate = Vmax/2
    let rate = enzyme::michaelis_menten(1.0, 10.0, 1.0);
    assert!((rate - 5.0).abs() < 0.01);
}

#[test]
fn test_michaelis_menten_saturation() {
    // At very high [S], rate -> Vmax
    let rate = enzyme::michaelis_menten(10000.0, 10.0, 1.0);
    assert!((rate - 10.0).abs() < 0.01);
}

#[test]
fn test_michaelis_menten_zero_substrate() {
    let rate = enzyme::michaelis_menten(0.0, 10.0, 1.0);
    assert!((rate - 0.0).abs() < f64::EPSILON);
}

#[test]
fn test_nernst_potassium() {
    // K+: [out]=4, [in]=155, z=1, T=310K -> ~-97 mV
    let e = membrane::nernst(310.0, 1, 4.0, 155.0);
    assert!(e < -80.0 && e > -110.0, "Nernst K+ = {e} mV");
}

#[test]
fn test_nernst_sodium() {
    // Na+: [out]=145, [in]=12, z=1, T=310K -> ~+67 mV
    let e = membrane::nernst(310.0, 1, 145.0, 12.0);
    assert!(e > 50.0 && e < 80.0, "Nernst Na+ = {e} mV");
}

#[test]
fn test_dose_response_at_ec50() {
    // At [L] = EC50, response = Emax/2
    let r = signal::dose_response(1.0, 1.0, 1.0, 1.0);
    assert!((r - 0.5).abs() < 0.01);
}

#[test]
fn test_dose_response_steep_hill() {
    // High Hill coefficient = steeper curve
    let r_low = signal::dose_response(0.5, 1.0, 1.0, 1.0);
    let r_high = signal::dose_response(0.5, 1.0, 1.0, 4.0);
    // Below EC50 with high n -> lower response
    assert!(r_high < r_low);
}

// --- Cross-module integration ---

#[test]
fn test_metabolic_energy_charge_consistency() {
    let m = MetabolicState::default();
    let ec = m.energy_charge();
    // Energy charge should be between 0 and 1
    assert!(ec >= 0.0 && ec <= 1.0);
    // Resting cell should have high energy charge (>0.9)
    assert!(ec > 0.9);
}

#[test]
fn test_bioenergy_and_metabolism_alignment() {
    // Both modules model energy: verify they agree on resting state
    let met = MetabolicState::default();
    let bio = BioenergyState::default();
    // Both should be aerobic at rest
    assert!(!met.is_anaerobic());
    assert!(!bio.is_anaerobic());
}

#[test]
fn test_membrane_potential_physiological_range() {
    let ions = IonicState::default();
    let vm = ions.resting_potential();
    // Resting membrane potential: -60 to -90 mV
    assert!(
        vm < -50.0 && vm > -100.0,
        "Resting potential {vm} mV out of physiological range"
    );
}
