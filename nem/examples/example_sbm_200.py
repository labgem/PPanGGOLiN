"""Generate a SBM dataset: 200 nodes, 3 clusters, bivariate Gaussian emissions."""

from generate import SBMData

sbm = SBMData(
    n=200,
    k=3,
    d=2,
    p_in=0.2,
    p_out=0.02,
    centers=[[0, 0], [4, 4], [0, 4]],
    sigma=1.0,
    seed=42,
)

sbm.export("sbm_200_3")
sbm.plot()

import matplotlib.pyplot as plt
plt.savefig("sbm_200_3.png", dpi=150, bbox_inches="tight")
plt.show()
