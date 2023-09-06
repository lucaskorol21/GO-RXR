import unittest
from GUI_GO import stringcheck, checkbracket

class TestDataFitting(unittest.TestCase):
    def test_checkstring(self):
        strings = ['.000', '5.135481', '6', '5.2131.1551', '15,1531','-5.14','423143#.34','abcd',
                   '12312313112312132122311212123.0','0.0']
        outputs = [False, True, True, False, False,False,False,False, True,True]

        for i,string in enumerate(strings):
            value = stringcheck(string)
            self.assertEqual(value,outputs[i])

    def test_checkbrackets(self):
        tests = ['()', '{}' ,'[]', '(asd)','{15csc}','[bfb%]', '([sdfsd])','([sdfsd]}', '((asdasd)', '(sadasdad))',
                 '(sdasdsadsa{dasdsad})', '[dasdsadsa[adasda]asdasd]','']
        outputs = [True, True, True, True, True, True, True, False, False, False,True,True, True]

        for i,test in enumerate(tests):
            value = checkbracket(test)
            self.assertEqual(value, outputs[i])
