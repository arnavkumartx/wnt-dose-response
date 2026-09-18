// ---------------------------------------------------------------------------
// Parameters for the CHIR99021 -> GSK3 -> Wnt/beta-catenin -> fate model.
//
// Units:  time = hours, small molecules = micromolar (uM),
//         proteins in the beta-catenin module = nanomolar (nM),
//         transcription factors in the GRN = dimensionless 0..1.
//
// Provenance tag on every value:
//   [M] measured  - from a published measurement
//   [D] derived   - computed from measured quantities (e.g. Cheng-Prusoff)
//   [F] fitted    - tuned so the model reproduces a published qualitative result
//   [A] assumed   - plausible placeholder; these are the ones to calibrate first
//
// GRN convention: every TF obeys  dX/dt = k*drive - d*X  with k == d, so each
// TF is bounded on [0,1] and reads directly as "fraction of maximal expression".
// ---------------------------------------------------------------------------

export const P = {

  // === Layer 0: CHIR99021 pharmacology ====================================
  chir: {
    Ki_gsk3b:   0.0099,  // uM   [M] GSK3-beta Ki (IC50 6.7-10 nM in vitro)
    Km_ATP:     15.0,    // uM   [M] GSK3-beta Km for ATP
    ATP_cell:   3000.0,  // uM   [M] free intracellular ATP ~1-5 mM
    Kp:         0.55,    // -    [A] cell:medium partition coefficient
    k_uptake:   6.0,     // /h   [A] equilibration rate into the cell
    hill:       1.0,     // -    [D] ATP-competitive single site -> n = 1

    // off-target / cytotoxicity arm (produces the upper limb of the bell curve)
    IC50_tox:   9.0,     // uM   [F] intracellular (~16 uM in the medium)
    hill_tox:   4.0,     // -    [F] steep: the toxic limb appears abruptly
    k_tox:      0.12,    // /h   [F] max death rate at saturating CHIR
    k_recover:  0.005    // /h   [F] slow regrowth: a toxic pulse is not undone
  },

  // === Layer 1: Wnt / beta-catenin ========================================
  wnt: {
    Kd_wnt:     0.5,     // nM   [A] WNT3A EC50 at Fzd/LRP6
    k_on_R:     2.0,     // /h   [A]
    k_off_R:    1.5,     // /h   [A]
    eta_R:      0.85,    // -    [F] max fractional DC inhibition by the receptor

    Axin0:      1.0,     // a.u. [D] basal Axin1 scaffold, normalised to 1
    k_syn_ax2:  0.45,    // a.u./h [F] AXIN2 is a direct Wnt target -> neg. feedback
    K_ax2:      4.0,     // nM   [F] TCF_B giving half-max AXIN2
    n_ax2:      2.0,     // -    [F]
    k_deg_ax2:  0.55,    // /h   [M] AXIN2 protein half-life ~1.2 h

    k_syn_b:    25.0,    // nM/h [M] beta-catenin synthesis (Lee 2003: 0.423 nM/min)
    k_phos:     254.0,   // nM/h [F] Vmax per unit Axin. Rescaled by 1/0.394 when
                         //          the explicit scaffold assembly replaced the lumped
                         //          term, so calibrated behaviour is unchanged at iwr=0.
    Km_b:       80.0,    // nM   [M] DC Km for beta-catenin
    k_deg_b:    0.09,    // /h   [M] DC-independent turnover (t1/2 ~ 8 h)
    k_deg_bp:   6.0,     // /h   [M] betaTrCP/proteasome clearance of phospho-bcat

    k_imp:      0.90,    // /h   [A] nuclear import
    k_exp:      6.00,    // /h   [A] nuclear export, maximal (fully active DC)
    exp_dc_gain:0.70,    // -    [F] fraction of export that is DC-dependent.
                         //          APC/Axin shuttle beta-catenin out of the
                         //          nucleus, so inhibiting GSK3 also RETAINS it.
    TCF_tot:    15.0,    // nM   [M] total TCF/LEF (Lee 2003)
    k_tcf_on:   0.25,    // /nM/h [A]
    k_tcf_off: 10.00,    // /h   [A]  -> Kd(beta-catenin:TCF) = 40 nM

    // cadherin sink: E-cadherin buffers cytoplasmic beta-catenin in hPSC;
    // SNAI1-driven EMT at the primitive streak releases it (positive feedback)
    CDH1_0:     120.0,   // nM   [A]
    emt_loss:   0.85,    // -    [F] max junctional E-cadherin lost to EMT
    K_snai_cdh: 0.35, n_snai_cdh: 3.0,
    k_mem_on:   0.25,    // /nM/h [A]
    k_mem_off:  1.20,    // /h   [A]

    // Autocrine WNT3/WNT5A. This is a SECRETED species with its own slow
    // clearance, not an instantaneous function of TBXT: after a CHIR pulse it
    // outlives the streak by days, which is the entire reason a Wnt inhibitor
    // step exists in the cardiac protocol.
    k_sec_wnt:    0.32,  // nM/h [F] max secretion rate
    k_clr_wnt:    0.08,  // /h   [F] clearance + dilution, t1/2 ~ 8.7 h
    gata4_wnt:    0.50,  // -    [F] nascent mesoderm keeps secreting after the
                         //          streak TFs are gone
    K_gata4_wnt:  0.40, n_gata4_wnt: 3.0,
    wnt_endo_max: 2.0,   // nM   [F] (legacy scale, kept for reference)
    K_endo_wnt:   0.45,  // -    [F] TBXT level at which WNT3 secretion switches on
    n_endo_wnt:   3.0    // -    [F] (a threshold, not a proportionality: sub-
                         //          threshold TBXT must not self-ignite the loop)
  },

  // === Layer 2: fate gene regulatory network ==============================
  grn: {
    // synthesis == degradation, so each TF saturates at 1.0
    d_oct4:  0.29, d_sox2: 0.29, d_nanog: 0.40,
    d_tbxt:  0.35, d_eomes: 0.30, d_mixl1: 0.33, d_snai1: 0.50,
    d_mesp1: 0.45, d_sox17: 0.20,
    d_gata4: 0.055,  // [F] slow: the stable cardiac-mesoderm latch
    d_nkx25: 0.12,
    d_tnnt2: 0.030,  // [F] sarcomeric protein, accumulates over days

    // --- Wnt (TCF_B, nM) thresholds ---
    K_tcf_tbxt:  3.0, n_tcf_tbxt:  4.0,  // [F] TBXT is the most Wnt-sensitive
    K_tcf_eomes: 3.0, n_tcf_eomes: 3.0,  // [F] EOMES needs more Wnt than TBXT
    K_tcf_mixl:  2.6, n_tcf_mixl:  3.0,
    K_tcf_nkx:   1.4, n_tcf_nkx:   4.0,  // [F] canonical Wnt REPRESSES NKX2-5
    K_tcf_snai:  3.2, n_tcf_snai:  3.0,

    // --- growth-factor thresholds ---
    K_act: 0.35, n_act: 2.0,   // Activin A / Nodal -> SMAD2/3
    K_bmp: 0.40, n_bmp: 2.0,   // BMP4 -> SMAD1/5/8

    // --- TF-TF thresholds ---
    K_os:        0.35, n_os:        2.0,  // OCT4:SOX2 heterodimer
    K_sox2_rep:  0.40, n_sox2_rep:  3.0,  // SOX2 -| TBXT   (Thomson 2011)
    sox2_rep_min: 0.30,                   // [F] SOX2 RAISES the Wnt threshold for
                                          //     TBXT, it does not veto it outright.
                                          //     A hard veto makes the switch unflippable.
    K_tbxt_rep:  0.35, n_tbxt_rep:  3.0,  // TBXT -| SOX2, NANOG
    K_str:       0.30, n_str:       3.0,  // streak TFs -> MESP1 / SOX17
    K_mesp1:     0.25, n_mesp1:     2.0,
    K_gata4:     0.30, n_gata4:     3.0,
    K_nkx25:     0.30, n_nkx25:     2.0,
    nkx_self:    0.85, K_nkx_self: 0.45, n_nkx_self: 3.0,  // [F] NKX2-5 autoregulates:
                                          //     once a cell is a committed cardiac
                                          //     progenitor, returning Wnt no longer
                                          //     reverses it.
    K_sox17_rep: 0.30, n_sox17_rep: 3.0,  // committed endoderm diverts MESP1

    // Cells that have adopted a definitive identity leave the streak state:
    // MESP1/NKX2-5+ cardiac mesoderm and SOX17+ endoderm are both T-negative.
    // Without this, returning autocrine Wnt re-ignites TBXT at day 10.
    commit_rep:  0.90, K_commit: 0.45, n_commit: 3.0,

    nanog_basal: 0.78,  // [F] Nodal-independent NANOG support in E8/mTeSR
    gata4_self:  0.90,  // [F] strength of GATA4 self-activation (the latch)
    // SOX17 and FOXA2 cross-activate into a self-reinforcing definitive-endoderm
    // network. Without this latch SOX17 peaks on day 1 and then collapses when
    // the CHIR spike washes out, so the model can make no endoderm on the real
    // protocol - the same failure the GATA4 latch fixes on the cardiac side.
    sox17_self:  0.85, K_sox17_self: 0.45, n_sox17_self: 3.0,
    tbxt_act:    0.35,  // [F] weight of SMAD2/3 on TBXT relative to Wnt
    leak:        0.02   // basal transcription
  },

  // === Cord blood HSPC module ==============================================
  // Sources (PubMed):
  //  Luis 2011 Cell Stem Cell 10.1016/j.stem.2011.07.017 - canonical Wnt
  //    regulates haematopoiesis in a DOSAGE-dependent way; five Apc allele
  //    levels give lineage-specific optima, explaining why gain-of-function
  //    studies report either HSC expansion or HSC depletion.
  //  Holmes 2008 Stem Cells 10.1634/stemcells.2007-0600 - GSK-3b inhibition in
  //    cord blood CD34+ raises beta-catenin, c-myc and HoxB4; DELAYS expansion
  //    but PRESERVES stem cell activity; raises CXCR4 -> stromal adherence.
  //  Ko 2011 Stem Cells 10.1002/stem.551 - GSK-3b inhibition lengthens cycle
  //    time, raises p57(CDKN1C), lowers cyclin D1, improves engraftment
  //    per repopulating unit rather than unit number.
  //  Shen 2014 Stem Cells Dev 10.1089/scd.2014.0230 - the same inhibitor given
  //    in vivo post-transplant IMPAIRS regeneration and blocks T development.
  hsc: {
    // Two direct Wnt targets with DIFFERENT thresholds. The non-monotonic
    // dose-response is not asserted anywhere - it emerges because HOXB4
    // (self-renewal) switches on below MYC (niche exit / differentiation).
    K_w_hoxb4: 1.9, n_w_hoxb4: 3.0,  // [F] lower threshold
    K_w_myc:   3.6, n_w_myc:   3.0,  // [F] higher threshold
    K_w_p57:   2.1, n_w_p57:   2.0,  // [M] CDKN1C up with GSK3i
    K_w_ccnd1: 7.0, n_w_ccnd1: 2.0,  // [F] cyclin D1 down with GSK3i
    K_w_cxcr4: 2.0, n_w_cxcr4: 2.0,  // [M] CXCR4 up -> stromal adherence
    K_w_lymph: 2.2, n_w_lymph: 3.0,  // [M] high Wnt blocks lymphoid output

    d_hoxb4: 0.25, d_myc: 0.35, d_cdkn1c: 0.20, d_ccnd1: 0.30, d_cxcr4: 0.18,
    d_gata2: 0.20, d_gata1: 0.15, d_spi1: 0.15, d_lymph: 0.12,
    d_cd34: 0.10, d_cd90: 0.10, d_cd45ra: 0.10,

    K_hoxb4:   0.40, n_hoxb4:   2.0,
    K_myc_rep: 0.45, n_myc_rep: 3.0,  // [F] MYC -| self-renewal (Wilson 2004)
    K_p57:     0.90, n_p57:     2.0,
    K_cyt:     0.40, n_cyt:     2.0,  // SCF/TPO/FLT3L drive
    K_diff:    0.40, n_diff:    3.0,
    K_ant:     0.40, n_ant:     3.0,  // GATA1 / PU.1 mutual antagonism
    K_stem:    0.40, n_stem:    2.0,
    K_cxcr4:   0.45, n_cxcr4:   2.0,

    k_prolif:   0.030,  // /h  [F] max division rate (~30 h doubling at full drive)
    d_stem:     0.012,  // /h  [F] loss of repopulating potential in culture
    k_div_loss: 1.50,   // -   [F] each division costs stemness: fast cycling
                        //         converges to the (low) self-renewal set point faster
    homing_floor: 0.30, // -   [F] engraftment with no CXCR4 upregulation
    leak: 0.02
  },

  pop: { V0: 1.0 }
};

// Every TF uses k == d so it is bounded on [0,1] and reads directly as
// "fraction of maximal expression". Filled in here rather than written twice.
for (const group of ['grn', 'hsc']) {
  for (const key of Object.keys(P[group])) {
    if (key.startsWith('d_')) P[group]['k_' + key.slice(2)] = P[group][key];
  }
}

export default P;
