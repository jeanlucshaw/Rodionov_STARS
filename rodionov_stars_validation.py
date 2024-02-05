"""
Reproduces the results from Figure 2 of the original paper (Rodionov, 2004), starting from a
digitized time series of the original figure.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rodionov_stars import regime_shift_rodionov

# ======================
# Import validation data
# ======================
rsi_validation = pd.read_csv('data/rodionov_2004_figure_2_RSI.csv', names=['time', 'RSI'])
pdo_validation = pd.read_csv('data/rodionov_2004_figure_2_PDO.csv', names=['time', 'PDO'])

# =============
# Run algorithm
# =============

rsi, regime_id, regime_mean, diff = regime_shift_rodionov(pdo_validation.PDO.values)

# ========================================================================
# Two-panel plot showing the results from the paper and the implementation
# ========================================================================

# Init
validation_color = 'limegreen'
_, ax = plt.subplots(2, figsize=(10, 6), sharex=True)

# ========
# PDO data

# draw the detected regimes as alternating colored bars
rbcm = np.array(['lightskyblue', 'lightcoral'])
colors = rbcm[np.int64(regime_id % 2)]
ax[0].bar(pdo_validation.time, pdo_validation.PDO, color=colors)

for y_ in rsi_validation.time:
    ax[0].axvline(y_, color='limegreen')
    ax[0].text(y_, 3.1, f'{y_:.0f}', ha='center', color='limegreen')

# ========
# RSI data
ax[1].plot(rsi_validation.time, rsi_validation.RSI, 'o', ms=5, mfc='none', mec='limegreen')

index_change = rsi > 0
ax[1].plot(pdo_validation.time.loc[index_change], rsi[index_change], 'o', mec='r', mfc='r')


# Plot parameters
for a_ in ax:
    a_.tick_params(which='both', top=False, right=False)
    a_.tick_params(which='minor', bottom=False, left=False)
ax[0].set(ylabel='PDO', ylim=(-3, 3))
ax[1].set(ylabel='RSI', ylim=(0, 1.6))

# Save
plt.savefig('rodionov_stars_validation.png', dpi=300)

