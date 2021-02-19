import pytest
import numpy as np
import unittest.mock as mock


def test_fastchem_molecules():
    from taurex_fastchem import FastChem
    from taurex.cache import OpacityCache
    nlayers = 10

    T = np.linspace(1000,500, nlayers)

    P = np.logspace(6,1, nlayers)
    
    molecules_active = ['H2O', 'CH4', 'NH3']
    test_name = 'H3N'

    with mock.patch.object(OpacityCache, "find_list_of_molecules") as mock_my_method:
        mock_my_method.return_value = molecules_active
        fc = FastChem()

    fc.initialize_chemistry(nlayers=nlayers, temperature_profile=T, pressure_profile=P)

    gases = fc.gases

    num_gases = len(gases)
    num_active = len(molecules_active)
    num_inactive = num_gases - num_active

    mix_profile = fc.mixProfile

    for m in molecules_active:
        assert m in gases
    for m in molecules_active:
        assert m in fc.activeGases

    for m in molecules_active:
        assert m not in fc.inactiveGases  


    assert test_name not in gases
    assert test_name not in fc.activeGases
    assert test_name not in fc.inactiveGases

    assert fc.mixProfile.shape == (num_gases, nlayers)
    assert fc.activeGasMixProfile.shape == (num_active, nlayers)
    assert fc.inactiveGasMixProfile.shape == (num_inactive, nlayers)







    