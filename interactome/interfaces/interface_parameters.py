
from interface import Interface

class InterfaceParameters:
    def __init__(self):
        pass

    def calcParams(self, code):
        interface = Interface(code)
        residues, atoms, contacts = interface.get()
        