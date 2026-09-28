import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("T_cst_delta/out", skiprows=1)
nl = 30
nt = data.shape[0] // nl
fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("viridis", nt)
for t in range(nt):
    layer = data[t * nl : (t + 1) * nl, 2]
    T = data[t * nl : (t + 1) * nl, 3]
    ax.plot(T, layer, color=cmap(t), marker="+", linestyle="-")
ax.set_xlabel("T")
ax.set_ylabel("Layer")
ax.set_title("Temperature profiles (T=cst, dT/dz=0 bot and top)")
# ax.set_xlim([19.999, 20.001])
# ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
# plt.savefig("T_profiles_T_cst.png", dpi=150)
plt.show()

# data = np.loadtxt("simple.cuda/log")
# nl = 30
# nt = data.shape[0] // nl
# fig, ax = plt.subplots(figsize=(8, 6))
# cmap = plt.get_cmap("viridis", nt)
# for t in range(nt):
#     layer = data[t * nl : (t + 1) * nl, 0]
#     T = data[t * nl : (t + 1) * nl, 1]
#     ax.plot(T, layer, color=cmap(t), marker="+", linestyle="-")
# ax.set_xlabel("T")
# ax.set_ylabel("Layer")
# ax.set_title("Temperature profiles (T=cst, dT/dz=0 bot and top)")
# ax.set_xlim([-1, 1])
# # ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
# plt.tight_layout()
# # plt.savefig("T_profiles_T_cst.png", dpi=150)
# plt.show()
