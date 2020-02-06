from taurex.chemistry import Chemistry
from .external.fastchem import PyDoubleFastChem
class FastChem(Chemistry):

    def __init__(self):
        super().__init__(self.__class__.__name__)

    