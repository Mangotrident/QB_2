
import numpy as np

class NoiseModel:
    """Noise utilities used in the pipeline (static disorder + gamma sweep)."""
    def __init__(self, sigma=0.0):
        self.sigma = float(sigma)

    def add_static_disorder(self, H, seed=None):
        if self.sigma <= 0:
            return H
        rng = np.random.default_rng(seed)
        H2 = H.copy()
        n = H.shape[0]
        diag_noise = rng.normal(0.0, self.sigma, size=n)
        for i in range(n):
            H2[i, i] = H2[i, i].real + diag_noise[i]
        return H2

    def create_gamma_sweep(self, start, end, step):
        # Accept both float step or 'log' dictionary
        if isinstance(step, dict) and step.get('scale') == 'log':
            # e.g., {'scale': 'log', 'num': 50}
            num = int(step.get('num', 50))
            # ensure >0 bounds
            start = max(start, 1e-6)
            end = max(end, start*1.1)
            return np.geomspace(start, end, num=num)
        else:
            return np.arange(start, end + 1e-12, step, dtype=float)
