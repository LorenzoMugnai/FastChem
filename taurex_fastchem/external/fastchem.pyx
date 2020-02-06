from libcpp.vector cimport vector
from libcpp.string cimport string
import collections

cdef extern from "fastchem.h" namespace "fastchem":
    cdef cppclass FastChem[T]:
        FastChem(string & model_parameter_file, int verbose) except +
        FastChem(const FastChem &obj) except +
        unsigned int getSpeciesNumber()
        unsigned int getElementNumber()
        unsigned int calcDensities(const double temperature, const double pressure,
                                   vector[double] & density_n_out, double& h_density_out, double& mean_molecular_weight_out)
        string getSpeciesName(const unsigned int species_index)
        string getSpeciesSymbol(const unsigned int species_index)

cdef class PyDoubleFastChem(object):
    
    cdef FastChem[double]* fchem

    def __cinit__(self,str filename,int verbose):
        cdef string model_file
        cdef int verb = verbose
        model_file =filename.encode('utf-8')
        self.fchem = new FastChem[double](model_file, verbose)
        #self.fchem = fchem
    
    def calcDensities(self, T, P):
        cdef vector[double] density_out
        cdef double h_density_out
        cdef double mean_mol_out
        cdef unsigned int result = self.fchem.calcDensities(T,P,density_out, h_density_out, mean_mol_out)

        return result ,density_out, h_density_out,mean_mol_out

    @property
    def numSpecies(self):
        return self.fchem.getSpeciesNumber()

    @property
    def numElements(self):
        return self.fchem.getElementNumber()

    def speciesName(self, index):
        cdef string result = self.fchem.getSpeciesName(index)
        return result.c_str().decode('utf-8')

    def speciesSymbol(self, index):
        cdef string result = self.fchem.getSpeciesSymbol(index)
        return result.c_str().decode('utf-8')

    def speciesIter(self):
        itera = range(self.numSpecies)

        for i in itera:
            yield self.speciesSymbol(i)

    def __dealloc__(self):
            del self.fchem