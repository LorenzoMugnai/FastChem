from taurex.chemistry import Chemistry
import taurex_fastchem.external.fastchem as fastchem
import tempfile
import os
import shutil
import numpy as np
from taurex.util.util import mass
class FastChem(Chemistry):

    parameter_template = """#element abundance file   
{}

#elements data file         
{}

#species data file    
{}

#accuracy of chemistry iteration
{}

#accuracy of pressure iteration
{}

#max error of Newton's method
{}

#max number of chemistry iterations 
{}  

#max number of pressure iterations                   
{}

#max number of Nelder-Mead iterations
{}
"""

    solar_metallicity = 0.0134

    default_elements = ['H', 'He','O', 'Al', 'Ar', 'C', 'Ca', 'Cl', 'Co',
                        'Cr', 'Cu', 'F', 'Fe', 'Ge', 'K', 'Mg', 'Mn', 'N', 'Na',
                        'Ne', 'Ni',  'P', 'S', 'Si',
                        'Ti', 'V', 'Zn','e-']

    default_abundances = [12.00, 10.93, 8.69, 6.45, 6.40, 8.43, 6.34, 5.50, 
                          4.99, 5.64, 4.19, 4.56, 7.50, 3.65, 
                           5.03, 7.60, 5.43, 7.83, 
                          6.24, 7.93, 6.22, 5.41, 7.12, 7.51, 
                          4.95, 3.93, 4.56, 13.1139]
    default_mass = [mass[el] for el in default_elements[:-1]]

    def __init__(self,H_He_ratio=0.083,elements=None,ratios_to_O=None, metallicity=0, base_profile='solar',
                 elements_datafile=None, species_datafile=None, chem_accuracy=1.0e-4, with_ions=False,
                 pressure_accuracy=1e-4, newton_error=1e-4, max_chem_iter=300,
                 max_press_iter=100, max_nedler_iter=20000):
        super().__init__(self.__class__.__name__)


        self._elements_datafile = elements_datafile
        self._species_datafile = species_datafile
        self._chem_accr = chem_accuracy
        self._press_accuracy = pressure_accuracy
        self._newton_error = newton_error 
        self._max_chem_iter = max_chem_iter
        self._max_press_iter = max_press_iter
        self._max_nedler_iter = max_nedler_iter
        self._with_ions = with_ions
        self._metallicity = metallicity
        self._h_he_ratio = H_He_ratio
        self._nonmetal = np.array(self.default_abundances[:2])
        self._nonmetal[1] = 12 + np.log10(self._h_he_ratio)
        self._O_abund = self.default_abundances[2]
        self._metal_elements = self.default_elements[3:-1]
        self._ratios = 10**(np.array(self.default_abundances[3:]) - self._O_abund)

        self._electron = self.default_abundances[-1]
        self.generate_abundances(elements=elements, ratios_to_O=ratios_to_O)

    

    def generate_abundances(self,elements=None,ratios_to_O=None):
        import math

        nonmetal = self._nonmetal[...]
        if elements is not None and ratios_to_O is not None:
            index = np.array([self._metal_elements.index(el) for el in elements])
            self._ratios[index] = np.array(ratios_to_O)
        
        ratios = np.log10(self._ratios[...])
        O_abund = math.log10(self._metallicity *
                                      (10**(self._O_abund-12.)))+12.

        nonmetal[1] = 12 + np.log10(self._h_he_ratio)

        print(ratios)

        metals = O_abund + ratios

        complete = np.concatenate((nonmetal,[O_abund],metals))


        return zip(self.default_elements[:-1],complete)

    def determine_molecules(self):
        param_file, element_file = self.generate_parameter_file()

        fchem = fastchem.PyDoubleFastChem(param_file,0)
        species = list(fchem.speciesIter())

        available = self.availableActive

        self._active_gases = [s for s in species if s in available]
        self._inactive_gases = [s for s in species if s not in available]


    def generate_element_file(self):
        elem_file = tempfile.NamedTemporaryFile('w',delete=False)

        elem_file.write('#### Cool stuff\n')

        for elm,val in self.generate_abundances():
            elem_file.write(f'{elm} {val}\n')
        elm_filename  = elem_file.name
        elem_file.close()
        return elm_filename


    def generate_parameter_file(self):
        param_file = tempfile.NamedTemporaryFile('w',delete=False)
        element_filename = self.generate_element_file()

        base_data_path = os.path.join(os.path.abspath(os.path.dirname(fastchem.__file__)),'data')
        elm_file = self._elements_datafile if self._elements_datafile is not None else os.path.join(base_data_path,'chemical_elements.dat')
        self.info('Elements data file used is in %s',elm_file)
        
        species_file = 'logK.dat' if self._with_ions else 'logK_wo_ions.dat'

        spec_file = self._species_datafile if self._species_datafile is not None else os.path.join(base_data_path,species_file)

        param_output = self.parameter_template.format(element_filename,
                                                    elm_file, spec_file, self._chem_accr, self._press_accuracy, 
                                                    self._newton_error, self._max_chem_iter, self._max_press_iter,
                                                    self._max_nedler_iter)
        param_file.write(param_output)
        param_filename = param_file.name
        param_file.close()

        return param_filename, element_filename
