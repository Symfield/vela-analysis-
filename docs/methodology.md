# Methodology

## Framework Overview

The Symfield Earth™ framework treats substrates as phase-coherent systems.

## Core Operators

### ⧖ (Tau Update)
```
g_new = g_old + α * ΔΦ
μ_new = μ_old * exp(-λ_eff * |ΔΦ|)
```

### FCI (Field Coherence Index)
```
FCI = mean(μ / (|g - g_mean| + ε))
```
Range 1.96-2.04 indicates coherent substrate.

### ITCI (Information-Theoretic Coherence Index)
```
ITCI = -Σ μ log(μ) / log(N)
```
Based on Tononi's Integrated Information Theory.

## Application to Vela

1. Extract density field from Chandra or synthetic ISM
2. Compute strain field ΔΦ = |∇(density)|
3. Evolve through 100 ⧖-updates from identity
4. Compute coherence metrics

## Validation Strategy

- **Random Noise**: Should give FCI < 1.8
- **Crab Nebula**: Should have < 4 rings (web structure)
- **Uniform Field**: Should give FCI ~2.0

## Assumptions

1. Density gradients proxy field curvature
2. Plasma can maintain long-term coherence
3. Fractal organization indicates non-random structure

## Limitations

- Primarily uses synthetic data
- Parameters derived from Earth system
- No temporal evolution data (snapshot only)
```

4. Click **"Commit new file"**

---

## **What You'll Have When Done:**
```
vela-analysis/
├── README.md ✅ (you have this)
├── LIMITATIONS ✅ (you have this)
├── vela_symfield_analysis.py ✅ (you have this)
├── requirements.txt ← add this
├── .gitignore ← add this
├── data/
│   └── README.md ← add this
├── docs/
│   └── methodology.md ← add this
└── results/
    └── README.md ← add this
