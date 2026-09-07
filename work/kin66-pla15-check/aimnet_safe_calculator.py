"""Verified process-local nonperiodic reference neighbors; quadratic memory/time."""
import torch
from aimnet.calculators import AIMNet2Calculator as _BaseCalculator
from aimnet.calculators.calculator import AdaptiveNeighborList

def reference_neighbors(self, positions, cell=None, pbc=None, batch_idx=None, fill_value=None):
    if cell is not None:
        raise ValueError("Reference neighbors support nonperiodic inputs only")
    n = len(positions)
    fill = n if fill_value is None else fill_value
    d = torch.cdist(positions, positions, compute_mode="donot_use_mm_for_euclid_dist")
    valid = d < self.cutoff
    valid.fill_diagonal_(False)
    if batch_idx is not None:
        valid &= batch_idx[:, None] == batch_idx[None, :]
    count = valid.sum(1).to(torch.int32)
    indices = torch.arange(n, device=positions.device, dtype=torch.int32).expand(n, n)
    indices = indices.masked_fill(~valid, fill).sort(dim=1).values
    return indices[:, :max(1, int(count.max()))].contiguous(), count, None

AdaptiveNeighborList.__call__ = reference_neighbors


class AIMNet2Calculator(_BaseCalculator):
    """Reference neighbors with full Coulomb range for legacy embedded models."""
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        if (self.metadata or {}).get("coulomb_mode") == "full_embedded":
            self._coulomb_cutoff = float("inf")
            self._update_lr_nblists()

NEIGHBOR_BACKEND = "explicit_torch_cdist_full_legacy_coulomb_v2"
