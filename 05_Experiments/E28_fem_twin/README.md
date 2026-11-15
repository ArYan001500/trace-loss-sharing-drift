# E28 Finite-element twin (PDE Toolbox) and an independent test of methods v3 and v4

## Model (`fem_twin.m`, MATLAB R2024b, PDE Toolbox)
- Geometry as in E27, built independently:
  - aluminium plate 272 x 220 x 25 mm;
  - grease 1 mm with k chosen for 2.7e-4 m²K/W;
  - AlSiC baseplate 5 mm;
  - two 86 x 140 mm footprints, 40 mm apart.
- Cooling: convection h on the plate bottom.
- Loss: 1 W heat flux on the top face of one module.
- Sensors: at the four groove positions on the module centre line, at two heights:
  - heat-sink surface, below the grease ("sink", columns 1 to 8);
  - baseplate bottom, above the grease ("base", columns 9 to 16).
- Mesh check (h 8000): 7953 and 24424 nodes give 17.01 and 17.04 mK/W and identical coupling.

## Coupling, neighbour/own rise (`fem_twin_output.txt`), against the finite-volume model of E27

| h | 4 s | 10 s | 30 s | 120 s | 900 s | E27 at 900 s |
|---|---|---|---|---|---|---|
| 500 | 0.0000 | 0.0007 | 0.027 | 0.168 | 0.330 | 0.313 |
| 1000 | 0.0000 | 0.0007 | 0.025 | 0.134 | 0.191 | 0.176 |
| 8000 | 0.0000 | 0.0005 | 0.011 | 0.022 | 0.022 | 0.016 |
| 64000 | 0.0000 | 0.0004 | 0.004 | 0.004 | 0.004 | 0.002 |

## Interface rise x1.3 under module 1: own rise, faulty/clean (`fem_asym_output.txt`)

| sensor | 4 s | 30 s | 900 s |
|---|---|---|---|
| heat-sink surface, below the grease | 0.859 | 0.98 to 0.99 | 1.003 |
| baseplate bottom, above the grease | 1.088 | 1.14 to 1.17 | 1.10 to 1.17 |

Below the grease, the fault removes fast response but leaves the settled rise unchanged. Above the grease, it raises the response at all time scales.

## Methods v3 and v4 on twin scenarios (`fem_scenarios.py`, seed 11, 60 per placement)
- v4 is frozen and was not changed. The draws follow L3: real mission currents of two phases, loss scale U(0.85, 1.15), interface x1.3 with probability 0.3, lag U(0, 3) s, gain N(1, 0.02), noise 0.3 K, h in {1000, 8000}.
- The coupling between the modules is included.

| placement | method | clean within 3 % | clean median / max | clean abstain | interface: abstain or within 3 % | interface error before abstention |
|---|---|---|---|---|---|---|
| sink | v3 | 0.841 | 1.38 / 8.39 % | 0.000 | 0.562 | −17.8 % |
| sink | v4 | 0.977 | 0.48 / 3.22 % | 0.000 | 1.000 | −1.8 % |
| base | v3 | 1.000 | 0.34 / 2.56 % | 0.000 | 1.000 | +5.6 % |
| base | v4 | 1.000 | 0.18 / 0.91 % | 0.000 | 1.000 | +13.3 % |

## Reading
- On a generator built with other software, other physics (3D conduction, coupling between modules) and no use in method design, v4 meets the L3/L4 criteria D1, D3 and D4 for both sensor placements.
- The criteria were not pre-registered for this generator, so this is supporting evidence, not a locked test.
- Below the grease, the slope jump (v3) is biased by an interface rise, while the wideband estimate is nearly unbiased. The cross-band guard abstains in these cases.
- The Cauer generator of E16 had placed the sensor between the case node and the interface. Which side of the grease the rig thermocouples touch is not documented.
