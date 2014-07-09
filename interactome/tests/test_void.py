
import unittest
import pytest

from interactome.interfaces.demo import *

def test_void():
    assert 2+2 == 4


class VoidTestCase(unittest.TestCase):
    def setUp(self):
        self.a = 1
        self.b = 2

    def tearDown(self):
        pass

    def test_void(self):
        self.assertEqual(self.a*2 + self.b, 4)


    def test_demo(self):
        self.assertEqual(demo(), 5)

# class StructureTestCase(unittest.TestCase):
#     def setUp():
