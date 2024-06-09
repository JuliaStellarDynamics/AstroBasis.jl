
# Bases

---
## Basis functions

`AstroBasis` currently offers four different types of bases.
### 3d (Spherical)
1. Clutton-Brock (1973), which is a match to the Plummer (1911) profile.
2. Hernquist & Ostriker (1992), which is a match to the Hernquist (1990) profile.

### 2d (Razor-Thin Discs)
1. Clutton-Brock (1972), which is a 2d match to a Plummer profile.
2. Kalnajs (1976), which is a spiral basis.

--

## Spherical bases

### Clutton-Brock (1973)
```@docs
AstroBasis.CB73Basis
```

### Hernquist (1992)
```@docs
AstroBasis.HernquistBasis
```


## Razor-thin bases

### Clutton-Brock (1972)
```@docs
AstroBasis.CB72Basis
```

### Kalnajs (1976)
```@docs
AstroBasis.K76Basis
```