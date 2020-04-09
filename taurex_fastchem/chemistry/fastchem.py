from taurex.chemistry import Chemistry
import taurex_fastchem.external.fastchem as fastchem
import tempfile
import os
import shutil
import numpy as np
from taurex.util.util import mass
from taurex.core import fitparam




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

    def __init__(self,H_He_ratio=0.083,selected_elements=None,ratio_elements=None,ratios_to_O=None, metallicity=1.0, base_profile='solar',
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

        self._elements = self.default_elements
        
        if selected_elements is not None:
            self._elements = [elem for elem in self.default_elements if elem in selected_elements]

        
        self._abundances = np.array([self.default_abundances[self.default_elements.index(elem)] 
                                       for elem in self._elements])
        
        self.info('Elements and thier abundances selected: %s',list(zip(self._elements,self._abundances)))



        self._nonmetal = np.array(self._abundances[:2])
        self._nonmetal[1] = 12 + np.log10(self._h_he_ratio)
        self._O_abund = self._abundances[2]
        self._metal_elements = self._elements[3:]
        self._ratios = 10**(np.array(self._abundances[3:]) - self._O_abund)

        self._electron = self.default_abundances[-1]
        self.generate_abundances(elements=ratio_elements, ratios_to_O=ratios_to_O)
        self.determine_molecules()
        self._gases = None
        self.add_ratio_params()
    

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

        #print(ratios)

        metals = O_abund + ratios

        complete = np.concatenate((nonmetal,[O_abund],metals))

        # print('FASTTTTTTTTTTT')
        # print(list(zip(self.default_elements[:-1],complete)))

        return zip(self._elements,complete)

    def determine_molecules(self):
        param_file, element_file = self.generate_parameter_file()

        fchem = fastchem.PyDoubleFastChem(param_file,0)
        species = [s.replace('1','') for s in fchem.speciesIter()]

        available = self.availableActive
        self._active_index = np.array([idx for idx,s in enumerate(species) if s in available])
        self._inactive_index = np.array([idx for idx,s in enumerate(species) if s not in available])

        self._active_gases = [s for s in species if s in available]
        self._inactive_gases = [s for s in species if s not in available]
        os.unlink(param_file)
        os.unlink(element_file)

    @property
    def activeGases(self):
        return self._active_gases

    @property
    def inactiveGases(self):
        return self._inactive_gases

    def initialize_chemistry(self, nlayers=100, temperature_profile=None,
                             pressure_profile=None, altitude_profile=None):
        """
        **Requires implementation**

        Derived classes should implement this to compute the active and
        inactive gas profiles

        Parameters
        ----------
        nlayers: int
            Number of layers in atmosphere

        temperature_profile: :obj:`array`
            Temperature profile in K, must have length ``nlayers``

        pressure_profile: :obj:`array`
            Pressure profile in Pa, must have length ``nlayers``

        altitude_profile: :obj:`array`
            Altitude profile in m, must have length ``nlayers``

        """

        CONST_K = 1.3806504e-16

        density = pressure_profile*10/(CONST_K*temperature_profile)

        param_file, element_file = self.generate_parameter_file()

        fchem = fastchem.PyDoubleFastChem(param_file,0)
        result ,density_out, h_density_out,mean_mol_out = \
            fchem.calcDensities(temperature_profile,pressure_profile*10)

        self._gases = np.array(density_out).T/density

        self.compute_mu_profile(nlayers)

        os.unlink(param_file)
        os.unlink(element_file)

    @fitparam(param_name='metallicity',param_latex='Z',default_bounds=[0.2,2.0],default_fit=False)
    def metallicity(self):
        return self._metallicity
    
    @metallicity.setter
    def metallicity(self,value):
        self._metallicity = value

    def add_ratio_params(self):

        for idx,element in enumerate(self._metal_elements):
            param_name = f'{element}_O_ratio'
            param_tex = f'{element}/O'

            def read_mol(self, idx=idx):
                return self._ratios[idx]

            def write_mol(self, value, idx=idx):
                self._ratios[idx] = value

            fget = read_mol
            fset = write_mol

            bounds = [1.0e-12, 0.1]

            default_fit = False
            self.add_fittable_param(param_name, param_tex, fget,
                                    fset, 'log', default_fit, bounds)




    @property
    def activeGasMixProfile(self):
        """
        **Requires implementation**

        Should return profiles of shape ``(nactivegases,nlayers)``. Active
        refers to gases that are actively absorbing in the atmosphere.
        Another way to put it these are gases where molecular cross-sections
        are used.

        """

        return self._gases[self._active_index]

    @property
    def inactiveGasMixProfile(self):
        """
        **Requires implementation**

        Should return profiles of shape ``(ninactivegases,nlayers)``.
        These general refer to gases: ``H2``, ``He`` and ``N2``


        """
        return self._gases[self._inactive_index]


    def generate_element_file(self):
        elem_file = tempfile.NamedTemporaryFile('w',delete=False)

        elem_file.write('#### Cool stuff\n')

        for elm,val in self.generate_abundances():
            self.info('Writing %s %s',elm,val)
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


    @classmethod
    def input_keywords(cls):
        return ['fastchem',]