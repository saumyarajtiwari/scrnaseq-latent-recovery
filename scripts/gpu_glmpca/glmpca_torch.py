import torch
import numpy as np
from scipy.io import mmread
import time
import sys

torch.manual_seed(42)
device = torch.device("cuda")

mtx_path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/glmpca_bridge/test_matrix.mtx"
out_path = sys.argv[2] if len(sys.argv) > 2 else "/tmp/glmpca_bridge/py_factors.csv"
n_epochs = int(sys.argv[3]) if len(sys.argv) > 3 else 5000

Y_np = mmread(mtx_path).toarray()
Y = torch.tensor(Y_np, dtype=torch.float32, device=device)
n_genes, n_cells = Y.shape
L = 30

a = torch.nn.Parameter(torch.zeros(n_genes, 1, device=device))
theta_raw = torch.nn.Parameter(torch.zeros(n_genes, 1, device=device))
U = torch.nn.Parameter(torch.randn(n_genes, L, device=device) * 0.01)
V = torch.nn.Parameter(torch.randn(n_cells, L, device=device) * 0.01)

optimizer = torch.optim.Adam([a, theta_raw, U, V], lr=0.05)

def nb_nll(Y, mu, theta):
    return -(torch.lgamma(Y + theta) - torch.lgamma(theta) - torch.lgamma(Y + 1)
              + theta * torch.log(theta / (theta + mu) + 1e-8)
              + Y * torch.log(mu / (theta + mu) + 1e-8)).sum()

t0 = time.time()
prev_loss = None
converged_at = None
for epoch in range(n_epochs):
    optimizer.zero_grad()
    eta = a + U @ V.T
    mu = torch.exp(eta).clamp(min=1e-6, max=1e6)
    theta = torch.nn.functional.softplus(theta_raw) + 1e-3
    loss = nb_nll(Y, mu, theta)
    loss.backward()
    optimizer.step()

    if epoch % 200 == 0:
        cur = loss.item()
        rel_change = abs(prev_loss - cur) / abs(prev_loss) if prev_loss else float('inf')
        print(f"epoch {epoch} loss {cur:.2f} rel_change {rel_change:.6f}")
        if rel_change < 1e-5 and converged_at is None and prev_loss is not None:
            converged_at = epoch
        prev_loss = cur

elapsed = time.time() - t0
print(f"time(s): {elapsed:.2f}")
print(f"converged_at epoch: {converged_at}")

np.savetxt(out_path, V.detach().cpu().numpy(), delimiter=",")
print("Python factors saved, dim:", V.shape)
