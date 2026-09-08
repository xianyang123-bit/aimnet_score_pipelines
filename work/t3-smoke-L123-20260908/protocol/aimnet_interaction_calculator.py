"""Interaction calculator with an explicit AIMNet2(2025) ensemble option."""
import torch
from aimnet_safe_calculator import AIMNet2Calculator

class InteractionCalculator:
    def __init__(self, model="aimnet2-2025", **kwargs):
        names = [f"aimnet2-b973c-2025-d3_{i}" for i in range(4)] if model == "aimnet2-2025-ensemble" else [model]
        self.calculators = [AIMNet2Calculator(name, **kwargs) for name in names]
        self.metadata = self.calculators[0].metadata
        self.model_names = names

    def __call__(self, data, **kwargs):
        # Calculators may prepare/mutate inputs: each member receives fresh tensors/arrays.
        outputs = []
        for calc in self.calculators:
            inputs = {k: v.clone() if isinstance(v, torch.Tensor) else v.copy() if hasattr(v, "copy") else v for k, v in data.items()}
            outputs.append(calc(inputs, **kwargs))
        result = dict(outputs[0])
        for key in ("energy", "forces", "charges"):
            if key in result:
                result[key] = torch.stack([o[key] for o in outputs]).mean(0)
        return result
