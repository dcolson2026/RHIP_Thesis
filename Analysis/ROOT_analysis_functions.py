from functools import lru_cache
import ROOT
import numpy as np
from particle import Particle
PI = np.pi

def calphi(y_momentum: float, x_momentum: float) -> float:
    return np.arctan2(y_momentum, x_momentum)

def caldeltaphi(trigger_phi: float, associate_phi: float) -> float:
    deltaphi = trigger_phi - associate_phi
    if deltaphi < -PI / 2:
        deltaphi += 2 * PI
    elif deltaphi >= 3 * PI / 2:
        deltaphi -= 2 * PI
    return deltaphi

def calpT(x_momentum: float, y_momentum: float) -> float:
    return (x_momentum**2 + y_momentum**2) ** 0.5

def caleta(x_momentum: float, y_momentum: float, z_momentum: float) -> float:
    # pseudorapidity is zero if no pz
    three_mom_mag = (x_momentum**2 + y_momentum**2 + z_momentum**2) ** 0.5
    if three_mom_mag == z_momentum:
        return 100 # arbitrarily high, doesn't matter
    three_mom_arg = (three_mom_mag + z_momentum) / (three_mom_mag - z_momentum)
    return np.log(three_mom_arg) / 2

def getParticleName(pdg_code: int) -> str:
    if pdg_code == 90:
        return "Initial process pdg 90"
    return Particle.from_pdgid(pdg_code)

def delta_phi_correlation(trigger_phi_list: list[float], associate_phi_list: list[float], histogram) -> None:
    for trig_phi in trigger_phi_list:
        for asso_phi in associate_phi_list:
            dphi = caldeltaphi(trig_phi, asso_phi)
            histogram.Fill(dphi)

def particle_lineage_dfs(initial_particle_index: int, daughter_1_list: list[int], daughter_2_list: list[int], total_daughters_of_initial_list: list[int]) -> None:
    daughter_1_index = daughter_1_list[initial_particle_index]
    daughter_2_index = daughter_2_list[initial_particle_index]

    if daughter_1_index == daughter_2_index == 0:
        return
    elif daughter_1_index == daughter_2_index:
        total_daughters_of_initial_list.append(daughter_1_index)
        particle_lineage_dfs(daughter_1_index, daughter_1_list, daughter_2_list, total_daughters_of_initial_list)
    else:
        if daughter_1_index != 0:
            total_daughters_of_initial_list. append(daughter_1_index)
            particle_lineage_dfs(daughter_1_index, daughter_1_list, daughter_2_list, total_daughters_of_initial_list)
        if daughter_2_index != 0:
            total_daughters_of_initial_list.append(daughter_2_index)
            particle_lineage_dfs(daughter_2_index, daughter_1_list, daughter_2_list, total_daughters_of_initial_list)

@lru_cache(maxsize=None)
def _pdg_has_charm(pdg_code: int) -> bool:
    if pdg_code == 0:
        return False
    if abs(pdg_code) == 4:
        return True

    try:
        particle = Particle.from_pdgid(pdg_code)
    except Exception:
        return False

    quark_content = getattr(particle, "quarks", None)
    if quark_content and "c" in str(quark_content).lower():
        return True

    flavor_content = getattr(particle, "flavor", None)
    if flavor_content and "c" in str(flavor_content).lower():
        return True

    return False

def following_charm_dfs(initial_particle_index: int, pdg_list: list[int], daughter_1_list: list[int], daughter_2_list: list[int], ordered_charm_daughters_of_initial_list: list[int]) -> None:
    daughter_1_index = daughter_1_list[initial_particle_index]
    daughter_2_index = daughter_2_list[initial_particle_index]
    daughter_1_pdg = pdg_list[daughter_1_index] if daughter_1_index else 0
    daughter_2_pdg = pdg_list[daughter_2_index] if daughter_2_index else 0

    # returns if this is a final state particle that undergoes no more decays or interactions ie no daughters
    if daughter_1_index == daughter_2_index == 0:
        return
    # in this case, there is only one daughter particle. exact copy but with changed momentum eg in shower
    elif 0 < daughter_1_index == daughter_2_index:
        if _pdg_has_charm(daughter_1_pdg):
            ordered_charm_daughters_of_initial_list.append(daughter_1_index)
            following_charm_dfs(daughter_1_index, pdg_list, daughter_1_list, daughter_2_list, ordered_charm_daughters_of_initial_list)
    elif daughter_1_index > 0 and daughter_2_index == 0:
        if _pdg_has_charm(daughter_1_pdg):
            ordered_charm_daughters_of_initial_list.append(daughter_1_index)
            following_charm_dfs(daughter_1_index, pdg_list, daughter_1_list, daughter_2_list, ordered_charm_daughters_of_initial_list)
    # list of daughters case
    elif 0 < daughter_1_index < daughter_2_index:
        for i in range(daughter_1_index, daughter_2_index + 1):
            if _pdg_has_charm(pdg_list[i]):
                ordered_charm_daughters_of_initial_list.append(i)
                following_charm_dfs(i, pdg_list, daughter_1_list, daughter_2_list, ordered_charm_daughters_of_initial_list)
    # two separate decays, not a list. check individually
    elif 0 < daughter_2_index < daughter_1_index:
        if _pdg_has_charm(daughter_1_pdg):
            ordered_charm_daughters_of_initial_list.append(daughter_1_index)
            following_charm_dfs(daughter_1_index, pdg_list, daughter_1_list, daughter_2_list, ordered_charm_daughters_of_initial_list)
        if _pdg_has_charm(daughter_2_pdg):
            ordered_charm_daughters_of_initial_list.append(daughter_2_index)
            following_charm_dfs(daughter_2_index, pdg_list, daughter_1_list, daughter_2_list, ordered_charm_daughters_of_initial_list)

_VON_MISES_MIN = 1e-6
_FIT_RANGE_EPS = 1e-5


def _positive_or_min(value: float, default: float) -> float:
    """Return value if it is reasonably positive, otherwise a safe default."""
    return value if value > _VON_MISES_MIN else default


def create_von_mises_fit(
    name: str,
    x_min: float,
    x_max: float,
    mu_guess: float,
    amplitude_guess: float,
    baseline_guess: float,
):
    """Create a von Mises TF1: f(x) = A * exp(kappa cos(x-mu)) + baseline."""
    import ROOT

    func = ROOT.TF1(
        name,
        "[0]*TMath::Exp([2]*TMath::Cos(x-[1])) + [3]",
        x_min,
        x_max,
    )

    amp_init = _positive_or_min(amplitude_guess, 1.0)
    base_init = _positive_or_min(baseline_guess, 0.1)

    func.SetParNames("Amp", "Mu", "Kappa", "Baseline")

    # Amplitude
    func.SetParameter(0, amp_init)
    func.SetParLimits(0, 0.0, max(amp_init * 10.0, 1.0))

    # Mean angle
    func.SetParameter(1, mu_guess)
    # +- ~40 degrees is a reasonable prior window
    func.SetParLimits(1, mu_guess - 0.7, mu_guess + 0.7)

    # Concentration (kappa)
    func.SetParameter(2, 1.0)
    func.SetParLimits(2, 0.01, 100.0)

    # Baseline (pedestal)
    func.SetParameter(3, base_init)
    func.SetParLimits(3, 0.0, max(base_init * 10.0, 1.0))

    return func


# def fit_von_mises_region(
#     histogram,
#     fit_name: str,
#     x_min: float,
#     x_max: float,
#     mu_guess: float,
#     line_color: int,
#     min_entries: int = 10,
# ):
#     """
#     Fit the provided histogram within a window using a von Mises function.

#     Returns (TF1, chi2, ndf) if the fit succeeds, otherwise None.
#     """
#     import ROOT

#     axis = histogram.GetXaxis()
#     hist_xmin = axis.GetXmin()
#     hist_xmax = axis.GetXmax()

#     fit_xmin = max(x_min, hist_xmin)
#     fit_xmax = min(x_max, hist_xmax)

#     if fit_xmax <= fit_xmin:
#         print(f"[von_mises] Bad fit range: {fit_xmin:.3f} – {fit_xmax:.3f}")
#         return None

#     low_bin = max(1, axis.FindBin(fit_xmin + _FIT_RANGE_EPS))
#     high_bin = min(axis.FindBin(fit_xmax - _FIT_RANGE_EPS), histogram.GetNbinsX())
#     if high_bin < low_bin:
#         print(f"[von_mises] high_bin < low_bin ({high_bin} < {low_bin})")
#         return None

#     entries_in_window = histogram.Integral(low_bin, high_bin)
#     if entries_in_window < min_entries:
#         print(
#             f"[von_mises] Not enough entries in window "
#             f"({entries_in_window} < {min_entries}) for {fit_name}"
#         )
#         return None

#     # Crude amplitude & baseline guesses from the local window
#     local_max = 0.0
#     local_sum = 0.0
#     local_bins = 0
#     for bin_idx in range(low_bin, high_bin + 1):
#         content = histogram.GetBinContent(bin_idx)
#         local_max = max(local_max, content)
#         local_sum += content
#         local_bins += 1

#     baseline_guess = (local_sum / local_bins) if local_bins else 0.1

#     fit = create_von_mises_fit(
#         fit_name,
#         fit_xmin,
#         fit_xmax,
#         mu_guess,
#         local_max,
#         baseline_guess,
#     )
#     fit.SetRange(fit_xmin, fit_xmax)
#     fit.SetLineColor(line_color)
#     fit.SetLineWidth(3)
#     fit.SetNpx(600)

#     # "R" = use range, "Q" = quiet, "S" = return fit result
#     fit_status = histogram.Fit(fit, "RQS+", "", fit_xmin, fit_xmax)

#     # PyROOT 6: Fit returns a TFitResultPtr with Status()
#     if hasattr(fit_status, "Status"):
#         status = fit_status.Status()
#         if status != 0:
#             print(f"[von_mises] Fit failed for {fit_name}, status = {status}")
#             return None
#     else:
#         # Older ROOT: Fit returns an int
#         if fit_status != 0:
#             print(f"[von_mises] Fit failed for {fit_name}, status = {fit_status}")
#             return None

#     chi2 = fit.GetChisquare()
#     ndf = fit.GetNDF()
#     print(
#         f"[von_mises] Fit OK for {fit_name}: "
#         f"chi2/ndf = {chi2:.1f}/{int(ndf) if ndf else 0}"
#     )
#     return fit, chi2, ndf

# def fit_von_mises_region(
#     histogram,
#     fit_name: str,
#     x_min: float,
#     x_max: float,
#     mu_guess: float,
#     line_color: int,
#     min_entries: int = 10,
# ):
#     """
#     Fit the provided histogram within a window using a PERIODIC von Mises shape.

#     The fitted function includes image peaks at mu, mu-2pi, and mu+2pi so that
#     a peak near a boundary still contributes correctly through the periodic wrap.

#     Returns (TF1, chi2, ndf) if the fit succeeds, otherwise None.

#     Parameter convention is kept compatible with the rest of your script:
#         par[0] = amplitude
#         par[1] = baseline
#         par[2] = kappa
#     """
#     import ROOT
#     import math

#     axis = histogram.GetXaxis()
#     hist_xmin = axis.GetXmin()
#     hist_xmax = axis.GetXmax()

#     fit_xmin = max(x_min, hist_xmin)
#     fit_xmax = min(x_max, hist_xmax)

#     if fit_xmax <= fit_xmin:
#         print(f"[von_mises] Bad fit range: {fit_xmin:.3f} – {fit_xmax:.3f}")
#         return None

#     low_bin = max(1, axis.FindBin(fit_xmin + _FIT_RANGE_EPS))
#     high_bin = min(axis.FindBin(fit_xmax - _FIT_RANGE_EPS), histogram.GetNbinsX())
#     if high_bin < low_bin:
#         print(f"[von_mises] high_bin < low_bin ({high_bin} < {low_bin})")
#         return None

#     entries_in_window = histogram.Integral(low_bin, high_bin)
#     if entries_in_window < min_entries:
#         print(
#             f"[von_mises] Not enough entries in window "
#             f"({entries_in_window} < {min_entries}) for {fit_name}"
#         )
#         return None

#     # Crude amplitude & baseline guesses from the local window
#     local_max = 0.0
#     local_sum = 0.0
#     local_bins = 0
#     for bin_idx in range(low_bin, high_bin + 1):
#         content = histogram.GetBinContent(bin_idx)
#         local_max = max(local_max, content)
#         local_sum += content
#         local_bins += 1

#     baseline_guess = (local_sum / local_bins) if local_bins else 0.1
#     amplitude_guess = max(local_max - baseline_guess, 0.1)

#     twopi = 2.0 * math.pi

#     def periodic_von_mises(x, par):
#         phi = x[0]
#         amplitude = par[0]
#         baseline = par[1]
#         kappa = par[2]

#         # Include central peak plus nearest periodic images.
#         # This is the key change: the fit "knows" about wrap-around.
#         vm0  = math.exp(kappa * math.cos(phi - mu_guess))
#         vm_p = math.exp(kappa * math.cos(phi - (mu_guess + twopi)))
#         vm_m = math.exp(kappa * math.cos(phi - (mu_guess - twopi)))

#         return baseline + amplitude * (vm0 + vm_p + vm_m)

#     fit = ROOT.TF1(
#         fit_name,
#         periodic_von_mises,
#         fit_xmin,
#         fit_xmax,
#         3,
#     )

#     fit.SetParNames("Amplitude", "Baseline", "Kappa")
#     fit.SetParameters(amplitude_guess, baseline_guess, 2.0)

#     fit.SetParLimits(0, 0.0, max(10.0 * max(local_max, 1.0), 1e6))
#     fit.SetParLimits(1, 0.0, max(10.0 * max(local_max, 1.0), 1e6))
#     fit.SetParLimits(2, 1e-3, 100.0)

#     fit.SetRange(fit_xmin, fit_xmax)
#     fit.SetLineColor(line_color)
#     fit.SetLineWidth(3)
#     fit.SetNpx(1200)

#     # "R" = use range, "Q" = quiet, "S" = return fit result
#     fit_status = histogram.Fit(fit, "RQS+", "", fit_xmin, fit_xmax)

#     # PyROOT 6: Fit returns a TFitResultPtr with Status()
#     if hasattr(fit_status, "Status"):
#         status = fit_status.Status()
#         if status != 0:
#             print(f"[von_mises] Fit failed for {fit_name}, status = {status}")
#             return None
#     else:
#         # Older ROOT: Fit returns an int
#         if fit_status != 0:
#             print(f"[von_mises] Fit failed for {fit_name}, status = {fit_status}")
#             return None

#     chi2 = fit.GetChisquare()
#     ndf = fit.GetNDF()
#     print(
#         f"[von_mises] Fit OK for {fit_name}: "
#         f"chi2/ndf = {chi2:.1f}/{int(ndf) if ndf else 0}"
#     )
#     return fit, chi2, ndf

def fit_double_von_mises_periodic(
    histogram,
    fit_name: str,
    line_color: int,
    min_entries: int = 20,
):
    """
    Fit the full histogram with a periodic double-von-Mises model:

        f(phi) = baseline
               + bin_width * Y_NS / (2*pi*I0(kappa_NS)) * exp(kappa_NS * cos(phi))
               + bin_width * Y_AS / (2*pi*I0(kappa_AS)) * exp(kappa_AS * cos(phi - pi))

    Physics meaning of parameters:
        par[0] = baseline
        par[1] = near-side yield
        par[2] = near-side kappa
        par[3] = away-side yield
        par[4] = away-side kappa

    Near side is fixed at mu = 0.
    Away side is fixed at mu = pi.

    Returns (TF1, chi2, ndf) if fit succeeds, otherwise None.
    """
    import ROOT
    import math

    axis = histogram.GetXaxis()
    fit_xmin = axis.GetXmin()
    fit_xmax = axis.GetXmax()

    total_entries = histogram.Integral()
    if total_entries < min_entries:
        print(
            f"[double_vm] Not enough entries "
            f"({total_entries} < {min_entries}) for {fit_name}"
        )
        return None

    # Crude initialization from the histogram itself
    nbins = histogram.GetNbinsX()
    bin_width = axis.GetBinWidth(1) if nbins > 0 else 1.0
    hist_max = histogram.GetMaximum()
    hist_mean = total_entries / nbins if nbins > 0 else 0.0
    positive_bin_contents = [
        histogram.GetBinContent(bin_idx)
        for bin_idx in range(1, nbins + 1)
        if histogram.GetBinContent(bin_idx) > 0.0
    ]
    baseline_guess = (
        min(positive_bin_contents)
        if positive_bin_contents
        else max(0.0, min(hist_mean, 0.5 * hist_max))
    )

    # Split the histogram into rough near/away windows only for starting guesses.
    # This does NOT define the fit model; the fit is still global.
    pi = math.pi
    eps = 1e-6

    near_low_bin_1 = axis.FindBin(-pi/2 + eps)
    near_high_bin_1 = axis.FindBin(pi/2 - eps)

    away_low_bin = axis.FindBin(pi/2 + eps)
    away_high_bin = axis.FindBin(3*pi/2 - eps)

    near_guess = histogram.Integral(near_low_bin_1, near_high_bin_1)
    away_guess = histogram.Integral(away_low_bin, away_high_bin)

    # Prevent pathological zero seeds
    near_guess = max(near_guess, max(1.0, 0.25 * total_entries))
    away_guess = max(away_guess, max(1.0, 0.25 * total_entries))

    def double_vm_periodic(x, par):
        phi = x[0]

        baseline = par[0]
        yield_ns = par[1]
        kappa_ns = par[2]
        yield_as = par[3]
        kappa_as = par[4]

        # Histograms are fitted as counts/bin. The bin-width factor keeps
        # yield parameters in true integrated counts over delta phi.
        ns_norm = (
            yield_ns * bin_width
            / (2.0 * math.pi * ROOT.TMath.BesselI0(kappa_ns))
        )
        as_norm = (
            yield_as * bin_width
            / (2.0 * math.pi * ROOT.TMath.BesselI0(kappa_as))
        )

        near = ns_norm * math.exp(kappa_ns * math.cos(phi))
        away = as_norm * math.exp(kappa_as * math.cos(phi - math.pi))

        return baseline + near + away

    fit = ROOT.TF1(
        fit_name,
        double_vm_periodic,
        fit_xmin,
        fit_xmax,
        5,
    )
    fit.SetParNames(
        "Baseline",
        "Yield_NS",
        "Kappa_NS",
        "Yield_AS",
        "Kappa_AS",
    )

    fit.SetParameters(
        baseline_guess,
        near_guess,
        2.0,
        away_guess,
        2.0,
    )

    # Sensible parameter limits
    upper_scale = max(10.0 * max(hist_max, 1.0), 1e6)
    fit.SetParLimits(0, 0.0, upper_scale)
    fit.SetParLimits(1, 0.0, max(10.0 * near_guess, 1e6))
    fit.SetParLimits(2, 1e-3, 100.0)
    fit.SetParLimits(3, 0.0, max(10.0 * away_guess, 1e6))
    fit.SetParLimits(4, 1e-3, 100.0)

    fit.SetRange(fit_xmin, fit_xmax)
    fit.SetLineColor(line_color)
    fit.SetLineWidth(3)
    fit.SetNpx(1200)

    # Full-range fit
    fit_status = histogram.Fit(fit, "RQS+", "", fit_xmin, fit_xmax)

    if hasattr(fit_status, "Status"):
        status = fit_status.Status()
        if status != 0:
            print(f"[double_vm] Fit failed for {fit_name}, status = {status}")
            return None
    else:
        if fit_status != 0:
            print(f"[double_vm] Fit failed for {fit_name}, status = {fit_status}")
            return None

    chi2 = fit.GetChisquare()
    ndf = fit.GetNDF()
    print(
        f"[double_vm] Fit OK for {fit_name}: "
        f"chi2/ndf = {chi2:.1f}/{int(ndf) if ndf else 0}"
    )
    return fit, chi2, ndf

def is_charged_pdg(pdgid: int) -> bool:
    pdgDB = ROOT.TDatabasePDG.Instance()
    part = pdgDB.GetParticle(pdgid)
    if not part:
        return False
    # Charge() returns charge in units of e/3 in ROOT's PDG DB
    return abs(part.Charge()) > 0.0

def get_bin_label(value, bin_defs):
    """
    bin_defs: list of tuples like (label, low_edge, high_edge)
    returns the matching label or None
    """
    for label, lo, hi in bin_defs:
        if lo < value <= hi:
            return label
    return None
