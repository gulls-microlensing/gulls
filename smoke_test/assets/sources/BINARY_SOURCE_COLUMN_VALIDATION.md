# Binary Source Catalog Column Validation

## Required Columns for Binary Source Runs

Based on analysis of `buildEvent.cpp`, `pllxLightcurveGenerator.cpp`, and `readSLList.cpp`, the following columns are required for binary source simulations when `MULTIPLE_SOURCES > 0`:

### 1. Fixed-Position Columns (accessed via slcat structure)
These columns are mapped by position in the catalog structure:
- **`mul`** (position 0): Proper motion in l direction (mas/yr) - used for parallax and relative motion
- **`mub`** (position 1): Proper motion in b direction (mas/yr) - used for parallax and relative motion
- **`Mass`** (position 12): Stellar mass (solar masses) - used for orbital calculations
- **`Radius`** (position 14): Stellar radius (solar radii) - used for finite source size calculations
- **`Dist`** (position 20): Distance (kpc) - used for Einstein radius and orbital separation calculations

### 2. Named Columns (accessed via datadict)
These columns are accessed by name through the `Sources->datadict` mapping:
- **`Is_Binary`**: Binary status flag
  - `0` = single star
  - `1` = primary star in binary system
  - `2+` = companion star in binary system
- **`ID`**: Unique identifier for each source
- **`primary_ID`**: ID of the primary star (for companions only)
  - Must match the `ID` of the primary star
  - Used to link companion stars to their primaries
- **`combined_logP`**: Log10 of orbital period (days)
  - Used in: `double P = pow(10, Sources->data[sc][Sources->datadict["combined_logP"]])/DAYINYR;`
  - Used to calculate orbital semi-major axis via Kepler's third law

### 3. Standard Magnitude Columns
All magnitude columns defined in the parameter file's `Nfilters` setting are required for photometry.

## Code Usage Details

### In `buildEvent.cpp` (lines 385-488):

#### Binary Identification (line 388):
```cpp
int isbinary = Sources->data[sn][Sources->datadict["Is_Binary"]];
```

#### Companion Linking (lines 398, 417):
```cpp
if(Sources->data[i][Sources->datadict["primary_ID"]] == Sources->data[sn][Sources->datadict["ID"]])
```

#### Orbital Parameter Calculation (lines 467-474):
```cpp
Event->scomp_rs.push_back((Sources->data[sc][Sources->RADIUS] * Rsun / Sources->data[sc][Sources->DIST]) / Event->thE);
double P = pow(10, Sources->data[sc][Sources->datadict["combined_logP"]])/DAYINYR;
double M1 = Sources->data[sn][Sources->datadict["Mass"]];
double M2 = Sources->data[sc][Sources->datadict["Mass"]];
double acomb = pow(P*P*(M1+M2), 1.0/3.0);
Event->scomp_s.push_back(acomb/(Event->thE * Sources->data[sn][Sources->DIST]));
```

### In `pllxLightcurveGenerator.cpp` (lines 104-117):

#### Binary Source Lightcurve Calculation (lines 107-112):
```cpp
double x2off = Event->scomp_s[0] * cos(Event->scomp_phase[0]*TO_RAD);
double y2off = Event->scomp_s[0] * sin(Event->scomp_phase[0]*TO_RAD) * cos(Event->scomp_inc[0]*TO_RAD);
xs2CoM = xsCoM + x2off * cos(Event->scomp_alpha[0]*TO_RAD) - y2off * sin(Event->scomp_alpha[0]*TO_RAD);
ys2Center = ysCenter + x2off * sin(Event->scomp_alpha[0]*TO_RAD) + y2off * cos(Event->scomp_alpha[0]*TO_RAD);
double amp2 = Event->vbm->BinaryMag2(a, q, xs2CoM, ys2Center, Event->scomp_rs[0]);
```

## Validation: smoke_source_binary_catalog.dat

### File Location
`smoke_test/assets/sources/smoke_source_binary_catalog.dat`

### Header Analysis
Header line (first non-comment line):
```
R062 Z087 Y106 J129 W146 H158 F184 K213 2MASS_Ks 2MASS_J 2MASS_H Bessell_U Bessell_B Bessell_V Bessell_I Bessell_R Kepler_Kp TESS DECam_z DECam_u DECam_g DECam_r DECam_i DECam_Y Gaia_G_EDR3 Gaia_BP_EDR3 Gaia_RP_EDR3 VISTA_Z VISTA_Y VISTA_J VISTA_H VISTA_Ks mul mub Vr U V W iMass CL age Teff logg pop Mass Mbol Radius [Fe/H] l b RA2000.0 DEC2000.0 Dist x y z A_Ks [alpha/Fe] ID Is_Binary primary_ID combined_logP
```

### Column Verification

#### ✅ Fixed-Position Columns Present:
- ✅ `mul` (column 33) - proper motion l
- ✅ `mub` (column 34) - proper motion b
- ✅ `Mass` (column 43) - stellar mass
- ✅ `Radius` (column 45) - stellar radius
- ✅ `Dist` (column 51) - distance

#### ✅ Named Columns Present:
- ✅ `ID` (column 58) - unique identifier
- ✅ `Is_Binary` (column 59) - binary status flag
- ✅ `primary_ID` (column 60) - primary star ID
- ✅ `combined_logP` (column 61) - log orbital period

### Sample Data Analysis

The catalog contains 6 sources organized as follows:

**Row 2:** Single source (Is_Binary=0, ID=1)
**Row 3:** Single source (Is_Binary=0, ID=2)  
**Row 4:** Single source (Is_Binary=0, ID=3)
**Row 5:** Single source (Is_Binary=0, ID=4)
**Row 6:** **Primary binary source** (Is_Binary=1, ID=5, combined_logP=2.0)
**Row 7:** **Companion to ID=5** (Is_Binary=2, ID=6, primary_ID=5, combined_logP=2.0)

This structure correctly demonstrates:
1. ✅ Primary marked with `Is_Binary=1`
2. ✅ Companion marked with `Is_Binary=2`
3. ✅ Companion's `primary_ID` (5) matches primary's `ID` (5)
4. ✅ Both have same `combined_logP` value (2.0 = 100-day orbital period)
5. ✅ Companion immediately follows primary in catalog (required by code comment: "companions immediately trail the primary")

## Summary

✅ **ALL REQUIRED COLUMNS ARE PRESENT** in `smoke_source_binary_catalog.dat`

✅ **CATALOG STRUCTURE IS VALID** for binary source testing:
- Contains both single and binary sources
- Binary system correctly configured with primary-companion linkage
- Companion follows primary in catalog order as required by code

The smoke test binary source catalog is properly formatted and contains all necessary columns for binary source microlensing simulations.

