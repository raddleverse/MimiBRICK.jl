using Mimi

# ============================================================================
# glaciers_mengel_component.jl
#
# Temperature-dependent-equilibrium glaciers-and-ice-caps (GIC) module, ported from
# Mengel et al. 2016 (PNAS 113:2597; github.com/matthiasmengel/sealevel,
# contributor_functions.py). This covers the SAME inventory as the component it
# replaces -- mountain glaciers AND small/peripheral ice caps (the ice sheets are
# handled separately by the Greenland and Antarctic components) -- so swapping it in
# does not drop any mass. It is an optional drop-in alternative to BRICK's default
# single-reservoir Wigley-Raper-Bakker `glaciers_small_icecaps` component, fixing its
# commit-everything pathology (any sustained T > teq melts the WHOLE reservoir).
#
#   equilibrium:  S_eq(T) = a * (1 - exp(-b*(T - T_lia)))   (0 at T = T_lia; -> a)
#   TWO-TIMESCALE transient (a single tau cannot be fast-early + slow-modern):
#     dS_fast/dt = (f*S_eq - S_fast)/tau_fast      (committed early melt, short tau)
#     dS_slow/dt = ((1-f)*S_eq - S_slow)/tau_slow  (modern tracking, long tau)
#     S = S_fast + S_slow
#
# T = BRICK's global_surface_temperature (GMT anomaly rel 1850-1900, K = degC).
# The emulator is driven by TOTAL temperature (anthropogenic + natural forcing
# -- they are NOT, and cannot be, disentangled, since the model sees one
# temperature).
#
# `gic_T_lia` = the glacier EQUILIBRIUM temperature ~ the Little-Ice-Age climate
# (negative, rel 1850-1900). It makes glaciers OUT of equilibrium at the
# 1850-1900 baseline (S_eq(0) = a(1 - exp(b*T_lia)) > 0), so they are COMMITTED
# to melt by post-LIA warming that predates the forcing window -- this SIMULATES
# the committed/disequilibrium early-20th-c melt directly (no external "natural"
# budget, no forcing split). It is the physical generalization of the
# single-reservoir `gsic_teq`.
#
# A sustained warming T* commits only S_eq(T*) < a -- a temperature-appropriate
# remnant survives (NO full depletion). The output variable `gsic_sea_level` is
# named identically to the default component so it wires into `global_sea_level`
# unchanged. Mengel published coeffs (19 glacier-model fits, median a ~ 0.47 m,
# b ~ 0.52 /K).
# ============================================================================

# Default Mengel glacier-emulator parameters for the BRICK constructors (used only
# when glacier_model = :mengel). Fixed starting-point values from the MimiBRICK-FM
# workflow, close to the Mengel-2016 19-glacier-model median (a ~ 0.47 m, b ~ 0.52
# /K). The FM workflow calibrates these per posterior draw, so override them via
# update_param! for calibrated runs rather than relying on these as a tuned central.
const MENGEL_GLACIER_DEFAULTS = (a=0.45, b=0.52, T_lia=-0.45, f=0.5, tau_fast=40.0, tau_slow=250.0, sl0=0.0)

# Internal helper: set all Mengel glacier parameters on model `m` (post-replace!).
# Called by every constructor that supports glacier_model=:mengel to avoid triplication.
function _apply_mengel_defaults!(m)
    update_param!(m, :glaciers_small_icecaps, :gic_a,        MENGEL_GLACIER_DEFAULTS.a)
    update_param!(m, :glaciers_small_icecaps, :gic_b,        MENGEL_GLACIER_DEFAULTS.b)
    update_param!(m, :glaciers_small_icecaps, :gic_T_lia,    MENGEL_GLACIER_DEFAULTS.T_lia)
    update_param!(m, :glaciers_small_icecaps, :gic_f,        MENGEL_GLACIER_DEFAULTS.f)
    update_param!(m, :glaciers_small_icecaps, :gic_tau_fast, MENGEL_GLACIER_DEFAULTS.tau_fast)
    update_param!(m, :glaciers_small_icecaps, :gic_tau_slow, MENGEL_GLACIER_DEFAULTS.tau_slow)
    update_param!(m, :glaciers_small_icecaps, :gic_sl0,      MENGEL_GLACIER_DEFAULTS.sl0)
end

@defcomp glaciers_mengel begin

    # --------------------
    # Model Parameters
    # --------------------

    gic_a        = Parameter()             # Asymptotic max contribution from the LIA state (m SLE).
    gic_b        = Parameter()             # Temperature sensitivity (1/K).
    gic_T_lia    = Parameter()             # Glacier equilibrium (LIA) temperature, rel 1850-1900 (°C, <0).
    gic_f        = Parameter()             # Fast-mode fraction of the equilibrium (0-1).
    gic_tau_fast = Parameter()             # Fast (committed) response timescale (yr).
    gic_tau_slow = Parameter()             # Slow (modern) response timescale (yr).
    gic_sl0      = Parameter()             # Initial cumulative sea level contribution at the start year (m).
    global_surface_temperature = Parameter(index=[time]) # Global mean surface temperature anomaly relative to pre-industrial (°C).

    # --------------------
    # Model Variables
    # --------------------

    gsic_fast      = Variable(index=[time])   # Fast-mode (committed) sea level contribution (m).
    gsic_slow      = Variable(index=[time])   # Slow-mode (modern-tracking) sea level contribution (m).
    gsic_sea_level = Variable(index=[time])   # Cumulative sea level rise from glaciers and small ice caps (m).

    # --------------------
    # Model Equations
    # --------------------

    function run_timestep(p, v, d, t)

        if is_first(t)

            # Split the initial sea level contribution between the two reservoirs.
            v.gsic_fast[t]      = p.gic_f * p.gic_sl0
            v.gsic_slow[t]      = (1 - p.gic_f) * p.gic_sl0
            v.gsic_sea_level[t] = p.gic_sl0

        else

            # Temperature-dependent equilibrium contribution (Mengel et al. 2016).
            S_eq = p.gic_a * (1 - exp(-p.gic_b * (p.global_surface_temperature[t-1] - p.gic_T_lia)))

            # Two-timescale relaxation toward the (partitioned) equilibrium.
            v.gsic_fast[t] = v.gsic_fast[t-1] + (p.gic_f * S_eq       - v.gsic_fast[t-1]) / p.gic_tau_fast
            v.gsic_slow[t] = v.gsic_slow[t-1] + ((1 - p.gic_f) * S_eq - v.gsic_slow[t-1]) / p.gic_tau_slow
            v.gsic_sea_level[t] = v.gsic_fast[t] + v.gsic_slow[t]
        end
    end
end
