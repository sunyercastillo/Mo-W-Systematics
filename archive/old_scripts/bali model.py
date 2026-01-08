import pandas as pd

# Each row: element; columns: DM, SubdSed, N-MORB, DcpX/fluid, Dgrt/fluid, Drutile/fluid
data = {
    "Element": ["Mo", "Ce", "W", "U", "Th", "Nb", "La"],
    "DM": [0.0250, 0.7720, 0.0024, 0.0047, 0.0137, 0.2100, 0.2340],
    "SubdSed": [5.49, 57.3, 1.69, 1.68, 0.691, 8.94, 28.8],
    "N-MORB": [0.408, 12.0, 0.038, 0.07115, 0.1871, 3.507, 3.895],
    "DcpX/fluid": [
        "40/10*(0.435*log10(fO2) + 0.42*log10(NaCl) - 1.8 - 1000/T + 4.8)",  # Mo
        2.000,  # Ce
        "60/10**(0.07*log10(fO2) - 4.7236*1000/T + 4.4271)",  # W
        "11/10**(2.681 + 0.1433*log10(fO2) + 0.594*log10(mol))",  # U
        1.190,  # Th
        0.172,  # Nb
        1.429,  # La
    ],
    "Dgrt/fluid": [
        "12/10*(0.435*log10(fO2) + 0.42*log10(NaCl) - 1.8 - 1000/T + 4.8)",  # Mo
        0.400,  # Ce
        "12/10**(0.07*log10(fO2) - 4.7236*1000/T + 4.4271)",  # W
        "40/10**(2.681 + 0.1433*log10(fO2) + 0.594*log10(mol))",  # U
        0.610,  # Th
        0.204,  # Nb
        0.204,  # La
    ],
    "Drutile/fluid": [
        "87670/10*(0.435*log10(fO2) + 0.42*log10(NaCl) - 1.8 - 1000/T + 4.8)",  # Mo
        2.000,  # Ce
        1.250,  # W
        "94/10**(1.7954 + 0.1433*log10(fO2) + 0.594*log10(mol))",  # U
        0.100,  # Th
        200.0,  # Nb
        1.250,  # La
    ],
    "is_formula": [
        True, False, True, True, False, False, False,  # DcpX/fluid
        True, False, True, True, False, False, False,  # Dgrt/fluid
        True, False, False, True, False, False, False,  # Drutile/fluid
    ]
}

df = pd.DataFrame(data)
