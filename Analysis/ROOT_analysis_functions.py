from functools import lru_cache

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

    # returns if this is a final state particle that undergoes no more decays or interactions
    if daughter_1_index == daughter_2_index == 0:
        return
    # in this case, there is only one daughter particle
    elif daughter_1_index == daughter_2_index:
        if _pdg_has_charm(daughter_1_pdg):
            ordered_charm_daughters_of_initial_list.append(daughter_1_index)
            following_charm_dfs(daughter_1_index, pdg_list, daughter_1_list, daughter_2_list, ordered_charm_daughters_of_initial_list)
    else:
        if daughter_1_index != 0 and _pdg_has_charm(daughter_1_pdg):
            ordered_charm_daughters_of_initial_list. append(daughter_1_index)
            following_charm_dfs(daughter_1_index, pdg_list, daughter_1_list, daughter_2_list, ordered_charm_daughters_of_initial_list)
        if daughter_2_index != 0 and _pdg_has_charm(daughter_2_pdg):
            ordered_charm_daughters_of_initial_list.append(daughter_2_index)
            following_charm_dfs(daughter_2_index, pdg_list, daughter_1_list, daughter_2_list, ordered_charm_daughters_of_initial_list)
