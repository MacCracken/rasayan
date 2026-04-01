//! Integration tests for rasayan — cross-module behavior and serde roundtrips.

use rasayan::energy::BioenergyState;
use rasayan::enzyme::{self, EnzymeParams};
use rasayan::glycolysis::{GlycolysisConfig, GlycolysisState};
use rasayan::membrane::{self, IonicState, MembranePermeability};
use rasayan::metabolism::MetabolicState;
use rasayan::signal::{self, SecondMessenger};
use rasayan::tca::{TcaConfig, TcaState};

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

// --- 0.2.0 Expanded kinetics integration ---

#[test]
fn test_inhibition_hierarchy() {
    // All inhibition types should reduce rate below uninhibited
    let base = enzyme::michaelis_menten(1.0, 10.0, 1.0);
    let comp = enzyme::competitive_inhibition(1.0, 1.0, 10.0, 1.0, 1.0);
    let uncomp = enzyme::uncompetitive_inhibition(1.0, 1.0, 10.0, 1.0, 1.0);
    let mixed = enzyme::mixed_inhibition(1.0, 1.0, 10.0, 1.0, 1.0, 1.0);
    assert!(comp < base);
    assert!(uncomp < base);
    assert!(mixed < base);
}

#[test]
fn test_reversible_haldane_consistency() {
    // At equilibrium, the rate should be ~zero
    let vf = 10.0;
    let kmf = 1.0;
    let vr = 5.0;
    let kmr = 2.0;
    let keq = enzyme::haldane_keq(vf, kmf, vr, kmr);
    // At equilibrium: [P]/[S] = Keq, so [P] = Keq * [S]
    let s = 1.0;
    let p = keq * s;
    let rate = enzyme::reversible_michaelis_menten(s, p, vf, kmf, vr, kmr);
    assert!(
        rate.abs() < 0.01,
        "Rate at equilibrium should be ~0, got {rate}"
    );
}

#[test]
fn test_enzyme_db_params_produce_valid_rates() {
    // Every enzyme in the database should produce a valid rate
    for enz in enzyme::KNOWN_ENZYMES {
        let params = enz.params(1e-6);
        assert!(params.validate().is_ok(), "Invalid params for {}", enz.name);
        let rate = params.rate(enz.km); // rate at Km
        assert!(rate > 0.0, "Zero rate for {} at Km", enz.name);
        assert!(rate.is_finite(), "Non-finite rate for {}", enz.name);
    }
}

#[test]
fn test_lineweaver_burk_and_eadie_hofstee_agree() {
    let data: Vec<(f64, f64)> = [0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
        .iter()
        .map(|&s| (s, enzyme::michaelis_menten(s, 8.0, 0.5)))
        .collect();
    let lb = enzyme::lineweaver_burk_fit(&data).unwrap();
    let eh = enzyme::eadie_hofstee_fit(&data).unwrap();
    // Both should recover the same Km and Vmax from ideal data
    assert!(
        (lb.km - eh.km).abs() < 0.05,
        "Km mismatch: LB={} EH={}",
        lb.km,
        eh.km
    );
    assert!(
        (lb.vmax - eh.vmax).abs() < 0.05,
        "Vmax mismatch: LB={} EH={}",
        lb.vmax,
        eh.vmax
    );
}

#[test]
fn test_kinetic_fit_serde_roundtrip() {
    let fit = enzyme::lineweaver_burk_fit(&[
        (0.5, enzyme::michaelis_menten(0.5, 10.0, 1.0)),
        (1.0, enzyme::michaelis_menten(1.0, 10.0, 1.0)),
        (5.0, enzyme::michaelis_menten(5.0, 10.0, 1.0)),
    ])
    .unwrap();
    let json = serde_json::to_string(&fit).unwrap();
    let fit2: enzyme::KineticFit = serde_json::from_str(&json).unwrap();
    assert!((fit2.km - fit.km).abs() < f64::EPSILON);
    assert!((fit2.vmax - fit.vmax).abs() < f64::EPSILON);
}

#[test]
fn test_arrhenius_temperature_sensitivity() {
    // Enzyme rate should roughly double for every 10K increase (typical biological range)
    let k1 = enzyme::arrhenius(1e10, 50_000.0, 300.0);
    let k2 = enzyme::arrhenius(1e10, 50_000.0, 310.0);
    let ratio = k2 / k1;
    // For Ea=50kJ/mol, ratio should be ~1.9 over 10K
    assert!(ratio > 1.5 && ratio < 3.0, "Ratio={ratio}");
}

// --- Glycolysis → TCA integration ---

#[test]
fn test_glycolysis_feeds_tca() {
    let glyco_config = GlycolysisConfig::default();
    let tca_config = TcaConfig::default();
    let mut glyco = GlycolysisState::default();
    let mut tca = TcaState::default();

    let mut total_nadh = 0.0;
    let mut total_atp = 0.0;

    // Run 50 seconds: glycolysis produces pyruvate, TCA consumes it
    for _ in 0..500 {
        let gflux = glyco.tick(&glyco_config, 6.0, 0.5, 700.0, 0.1);
        let tflux = tca.tick(&tca_config, glyco.pyruvate, 6.0, 0.5, 700.0, 0.1);

        // TCA consumed some pyruvate — deduct it
        glyco.pyruvate = (glyco.pyruvate - tflux.pyruvate_consumed).max(0.0);

        total_nadh += gflux.nadh_produced + tflux.nadh_produced;
        total_atp += gflux.net_atp + tflux.gtp_produced;
    }

    // Both pathways should have produced meaningful output
    assert!(total_nadh > 0.0, "Combined NADH should be positive");
    assert!(total_atp > 0.0, "Combined ATP+GTP should be positive");
    // TCA should have consumed some pyruvate
    assert!(
        glyco.pyruvate < 5.0,
        "Pyruvate should not accumulate unbounded when TCA is consuming it"
    );
}

#[test]
fn test_glycolysis_tca_pathway_validation() {
    assert!(GlycolysisState::default().validate().is_ok());
    assert!(GlycolysisConfig::default().validate().is_ok());
    assert!(TcaState::default().validate().is_ok());
    assert!(TcaConfig::default().validate().is_ok());
}
