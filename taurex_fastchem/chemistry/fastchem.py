from taurex.chemistry import AutoChemistry
import taurex_fastchem.external.fastchem as fastchem
import tempfile
import os
import shutil
import numpy as np
from taurex.util.util import mass
from taurex.core import fitparam


def element_count(mol):
    import re
    s=re.findall('([A-Z][a-z]?)([0-9]*)', mol)
    compoundweight = 0
    elems_dict = {}
    for element, count in s:
        if count == '':
            count = 1
        if element not in elems_dict:
            elems_dict[element] = count
        else:
            elems_dict[element] += count
    return elems_dict


class FastChem(AutoChemistry):

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

    def __init__(self,H_He_ratio=0.083,selected_elements=None, elements_abundance_file=None, ratio_elements=None,ratios_to_O=None, metallicity=1.0,
                 elements_datafile=None, species_datafile=None, chem_accuracy=1.0e-4, with_ions=False,
                 pressure_accuracy=1e-4, newton_error=1e-4, max_chem_iter=300,
                 max_press_iter=100, max_nedler_iter=20000, longdouble=False):
        super().__init__(self.__class__.__name__)

        self._gases = None
        self._mixprofile = None
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
        self._long_double = longdouble

        self._elements = self.default_elements
        
        self._abundances = np.array(self.default_abundances[:])

        if elements_abundance_file is not None:
            elements = []
            abundances = []
            with open(elements_abundance_file,'r') as f:
                for l in f:
                    if l.startswith('#'):
                        continue
                    elem ,val = l.split()
                    elem = elem.strip()
                    val = float(val)
                    elements.append(elem)
                    abundances.append(val)          

            H_index = elements.index('H')
            elements.insert(0,elements.pop(H_index))
            abundances.insert(0,abundances.pop(H_index))

            He_index = elements.index('He')
            elements.insert(1,elements.pop(He_index))
            abundances.insert(1,abundances.pop(He_index))

            O_index = elements.index('O')
            elements.insert(2,elements.pop(O_index))
            abundances.insert(2,abundances.pop(O_index))
            
            e_index = elements.index('e-')
            elements.append(elements.pop(e_index))
            abundances.append(abundances.pop(e_index))




            


            self._elements = elements
            self._abundances = np.array(abundances)

        if selected_elements is not None:
            elements = [elem for elem in self._elements if elem in selected_elements]
            abundances = np.array([self._abundances[self._elements.index(elem)] 
                                        for elem in self._elements])
            self._elements = elements
            self._abundances = abundances
            
        
        self.info('Elements and thier abundances selected: %s',list(zip(self._elements,self._abundances)))




        self._nonmetal = np.array(self._abundances[:2])
        self._nonmetal[1] = 12 + np.log10(self._h_he_ratio)
        self._O_abund = self._abundances[2]
        self._metal_elements = self._elements[3:]
        self._ratios = 10**(np.array(self._abundances[3:]) - self._O_abund)

        self._electron = self.default_abundances[-1]
        self.generate_abundances(elements=ratio_elements, ratios_to_O=ratios_to_O)
        self.determine_molecules()
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
        _new_gas_list = []

        # Since FastChem orders atoms alphabetically
        # We need to find molecules that are the same as
        # the available active ones but ordered different atomically
        # i.e NH3 == H3N

        for s in species:
            found = False
            gas_name = s
            elem_check = element_count(s)
            for av in available:
                if elem_check == element_count(av):
                    gas_name = av
                    break
            
            _new_gas_list.append(gas_name)
        
        self._gases = _new_gas_list
        self.determine_active_inactive()

        os.unlink(param_file)
        os.unlink(element_file)

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
        fchem = None
        if self._long_double:
            fchem = fastchem.PyLongDoubleFastChem(param_file, 0)
        else:
            fchem = fastchem.PyDoubleFastChem(param_file,0)
        result ,density_out, h_density_out,mean_mol_out = \
                fchem.calcDensities(temperature_profile,pressure_profile*10)

        self._mixprofile = (np.array(density_out).T).astype(np.float)/density

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
        import pkg_resources
        param_file = tempfile.NamedTemporaryFile('w',delete=False)
        element_filename = self.generate_element_file()

        base_data_path = pkg_resources.resource_filename('taurex_fastchem',os.path.join('data','input'))
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


    @property
    def gases(self):
        return self._gases
    
    @property
    def mixProfile(self):
        return self._mixprofile

    @classmethod
    def input_keywords(cls):
        return ['fastchem',]

    BIBTEX_ENTRIES = [
        """
        @article{fastchem,
            author = {Stock, Joachim W and Kitzmann, Daniel and Patzer, A Beate C and Sedlmayr, Erwin},
            title = "{FastChem: A computer program for efficient complex chemical equilibrium calculations in the neutral/ionized gas phase with applications to stellar and planetary atmospheres }",
            journal = {Monthly Notices of the Royal Astronomical Society},
            volume = {479},
            number = {1},
            pages = {865-874},
            year = {2018},
            month = {06},
            abstract = "{For the calculation of complex neutral/ionized gas-phase chemical equilibria, we present a semi-analytical, versatile, and efficient computer program, called FastChem. The applied method is based on the solution of a system of coupled non-linear (and linear) algebraic equations, namely the law of mass action and the element conservation equations including charge balance, in many variables. Specifically, the system of equations is decomposed into a set of coupled nonlinear equations in one variable each, which are solved analytically whenever feasible to reduce computation time. Notably, the electron density is determined by using the method of Nelder and Mead at low temperatures. The program is written in object-oriented C++ which makes it easy to couple the code with other programs, although a stand-alone version is provided. FastChem can be used in parallel or sequentially and is available under the GNU General Public License version 3 at https://github.com/exoclime/FastChem together with several sample applications. The code has been successfully validated against previous studies and its convergence behaviour has been tested even for extreme physical parameter ranges down to \\$100\\,\\mathrm\\{K\\}\\$ and up to \\$1000\\,\\mathrm\\{bar\\}\\$. FastChem converges stable and robust in even most demanding chemical situations, which posed sometimes extreme challenges for previous algorithms.}",
            issn = {0035-8711},
            doi = {10.1093/mnras/sty1531},
            url = {https://doi.org/10.1093/mnras/sty1531},
            eprint = {https://academic.oup.com/mnras/article-pdf/479/1/865/25126582/sty1531.pdf},
        }
        """
    ]