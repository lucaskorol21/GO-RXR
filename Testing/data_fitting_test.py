import unittest
import material_structure as ms

class TestDataFitting(unittest.TestCase):
    def test_ChangeSampleParams_element(self):
        sample = ms.slab(2)
        sample.addlayer(0,'SrTiO3',50)
        sample.addlayer(1, 'LaMnO3',10)

        parameters = [[0,'STRUCTURAL','ELEMENT','Sr', 'THICKNESS'], [1,'STRUCTURAL','ELEMENT','La', 'THICKNESS'],
                      [0,'STRUCTURAL','ELEMENT','Ti', 'DENSITY'], [1,'STRUCTURAL','ELEMENT','O', 'DENSITY'],
                      [0,'STRUCTURAL','ELEMENT','O', 'ROUGHNESS'],[1,'STRUCTURAL','ELEMENT','Mn', 'ROUGHNESS'],
                      [0,'STRUCTURAL','ELEMENT','Sr', 'LINKED ROUGHNESS'],[1,'STRUCTURAL','ELEMENT','Sr', 'LINKED ROUGHNESS']]

        x = [10,20,0.025,0.086,0.5,5,2.75,2.8]

        backS = {}
        scaleF ={}
        script = []
        orbitals = {}

        sample_new = 1
    def test_ChangeSampleParams_compound(self):
        pass
    def test_ChangeSampleParams_variation(self):
        pass
    def test_ChangeSampleParams_magnetic(self):
        pass
    def test_ChangeSampleParams_eShift(self):
        pass
    def test_ChangeSampleParams_eShift_mag(self):
        pass
    def test_ChangeSampleParams_orbitals(self):
        pass
    def test_ChangeSampleParams_scaling(self):
        pass
    def test_ChangeSampleParams_shift(self):
        pass
    def test_ChangeSampleParams_script(self):
        pass
    def test_read_script(self):
        pass