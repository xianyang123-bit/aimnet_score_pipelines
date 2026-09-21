"""Unchanged ligand-only FIRE optimizer; every receptor atom remains fixed."""
import numpy as np
import torch
from rdkit import Chem
from aimnet_safe_calculator import AIMNet2Calculator

EV_TO_KCAL_MOL = 23.060547830619

def ligand_tensors(mol: Chem.Mol, device: str) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    coord = torch.as_tensor(
        np.asarray(mol.GetConformer().GetPositions(), dtype=np.float32),
        device=device,
    ).unsqueeze(0)
    numbers = torch.tensor(
        [[atom.GetAtomicNum() for atom in mol.GetAtoms()]],
        dtype=torch.long,
        device=device,
    )
    charge = torch.tensor([Chem.GetFormalCharge(mol)], dtype=torch.long, device=device)
    return coord, numbers, charge

def gas_energy(
    calc: AIMNet2Calculator,
    coord: torch.Tensor,
    numbers: torch.Tensor,
    charge: torch.Tensor,
) -> float:
    with torch.no_grad():
        result = calc(
            {
                "coord": coord.squeeze(0),
                "numbers": numbers.squeeze(0),
                "charge": charge.float(),
            }
        )
    return float(result["energy"].reshape(-1)[0].detach().cpu())

def minimize_in_pocket(calc, coord, numbers, charge, pocket_coord, pocket_numbers,
                       pocket_charge, *, fmax, max_steps):
    """FIRE on ligand coordinates only; return the bound geometry and diagnostics."""
    x = coord.detach().clone()
    velocity = torch.zeros_like(x)
    dt, alpha, n_positive = 0.1, 0.1, 0
    full_numbers = torch.cat((pocket_numbers, numbers.squeeze(0)))
    total_charge = charge.float() + pocket_charge
    initial_energy = None
    for step in range(max_steps + 1):
        result = calc({
            "coord": torch.cat((pocket_coord, x.squeeze(0))).detach().clone(),
            "numbers": full_numbers.clone(), "charge": total_charge.clone(),
        }, forces=True)
        energy = float(result["energy"].reshape(-1)[0].detach().cpu())
        force = result["forces"].reshape(-1, 3)[len(pocket_numbers):].detach().unsqueeze(0)
        if not np.isfinite(energy) or not torch.isfinite(force).all():
            raise ValueError("Nonfinite complex energy or ligand forces")
        if initial_energy is None:
            initial_energy = energy
        max_force = float(force.norm(dim=-1).max().cpu())
        if max_force < fmax or step == max_steps:
            return x, energy, step, max_force, max_force < fmax, initial_energy
        with torch.no_grad():
            power = float((velocity * force).sum().cpu())
            if power > 0:
                velocity.mul_(1-alpha).add_(
                    alpha * velocity.norm().clamp_min(1e-12)
                    * force / force.norm().clamp_min(1e-12))
                n_positive += 1
                if n_positive > 5:
                    dt = min(dt * 1.2, 0.1)
                    alpha *= 0.99
            else:
                velocity.zero_()
                dt *= 0.8
                alpha, n_positive = 0.1, 0
            velocity.add_(dt * force)
            displacement = dt * velocity
            displacement.mul_(min(1.0, 0.1 / float(displacement.norm().clamp_min(1e-12).cpu())))
            x = (x + displacement).detach()
    raise AssertionError("unreachable")
