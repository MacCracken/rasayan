//! Integration tests for rasayan — cross-module behavior and serde roundtrips.

use rasayan::energy::BioenergyState;
use rasayan::enzyme::{self, EnzymeParams};
use rasayan::membrane::{self, IonicState, MembranePermeability};
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
fn test_membrane_permeability_serde_roundtrip() {
    let perm = MembranePermeability::default();
    let json = serde_json::to_string(&perm).unwrap();
    let perm2: MembranePermeability = serde_json::from_str(&json).unwrap();
    assert!((perm2.p_na - perm.p_na).abs() < f64::EPSILON);
    assert!((perm2.p_k - perm.p_k).abs() < f64::EPSILON);
    assert!((perm2.p_cl - perm.p_cl).abs() < f64::EPSILON);
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
    let rate = enzyme::michaelis_menten(1.0, 10.0, 1.0);
    assert!((rate - 5.0).abs() < 0.01);
}

#[test]
fn test_michaelis_menten_saturation() {
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
    let e = membrane::nernst(310.0, 1, 4.0, 155.0);
    assert!(e < -80.0 && e > -110.0, "Nernst K+ = {e} mV");
}

#[test]
fn test_nernst_sodium() {
    let e = membrane::nernst(310.0, 1, 145.0, 12.0);
    assert!(e > 50.0 && e < 80.0, "Nernst Na+ = {e} mV");
}

#[test]
fn test_dose_response_at_ec50() {
    let r = signal::dose_response(1.0, 1.0, 1.0, 1.0);
    assert!((r - 0.5).abs() < 0.01);
}

#[test]
fn test_dose_response_steep_hill() {
    let r_low = signal::dose_response(0.5, 1.0, 1.0, 1.0);
    let r_high = signal::dose_response(0.5, 1.0, 1.0, 4.0);
    assert!(r_high < r_low);
}

// --- Cross-module integration ---

#[test]
fn test_metabolic_energy_charge_consistency() {
    let m = MetabolicState::default();
    let ec = m.energy_charge();
    assert!((0.0..=1.0).contains(&ec));
    assert!(ec > 0.9);
}

#[test]
fn test_bioenergy_and_metabolism_alignment() {
    let met = MetabolicState::default();
    let bio = BioenergyState::default();
    assert!(!met.is_anaerobic());
    assert!(!bio.is_anaerobic());
}

#[test]
fn test_membrane_potential_physiological_range() {
    let ions = IonicState::default();
    let vm = ions.resting_potential();
    assert!(
        vm < -50.0 && vm > -100.0,
        "Resting potential {vm} mV out of physiological range"
    );
}

// --- Validation integration ---

#[test]
fn test_all_defaults_validate() {
    assert!(MetabolicState::default().validate().is_ok());
    assert!(BioenergyState::default().validate().is_ok());
    assert!(
        EnzymeParams {
            vmax: 10.0,
            km: 1.0,
            hill_n: 1.0,
            kcat: 100.0,
        }
        .validate()
        .is_ok()
    );
}

#[test]
fn test_goldman_with_custom_permeability() {
    let ions = IonicState::default();
    // All-K+ permeability should give Nernst K+ potential
    let k_only = MembranePermeability {
        p_na: 0.0,
        p_k: 1.0,
        p_cl: 0.0,
    };
    let vm = membrane::goldman(310.0, &ions, &k_only);
    let e_k = membrane::nernst(310.0, 1, ions.k_out, ions.k_in);
    // Goldman with only K+ permeability should approximate Nernst K+
    assert!(
        (vm - e_k).abs() < 1.0,
        "Goldman K+-only ({vm:.1}) should ~ Nernst K+ ({e_k:.1})"
    );
}

#[test]
fn test_nernst_zero_concentrations() {
    assert!(membrane::nernst(310.0, 1, 0.0, 155.0).abs() < f64::EPSILON);
    assert!(membrane::nernst(310.0, 1, 4.0, 0.0).abs() < f64::EPSILON);
    assert!(membrane::nernst(310.0, 0, 4.0, 155.0).abs() < f64::EPSILON);
}
