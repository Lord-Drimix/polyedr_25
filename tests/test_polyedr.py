import unittest
from unittest.mock import patch, mock_open

from shadow.polyedr import Polyedr


class TestPolyedr(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        fake_file_content = """200.0	45.0	45.0	30.0
8	4	16
-0.5	-0.5	0.5
-0.5	0.5	0.5
0.5	0.5	0.5
0.5	-0.5	0.5
-0.5	-0.5	-0.5
-0.5	0.5	-0.5
0.5	0.5	-0.5
0.5	-0.5	-0.5
4	5    6    2    1
4	3    2    6    7
4	3    7    8    4
4	1    4    8    5"""
        fake_file_path = 'data/holey_box.geom'
        with patch('shadow.polyedr.open'.format(__name__),
                   new=mock_open(read_data=fake_file_content)) as _file:
            self.polyedr = Polyedr(fake_file_path)
            _file.assert_called_once_with(fake_file_path)

    def test_num_vertexes(self):
        self.assertEqual(len(self.polyedr.vertexes), 8)

    def test_num_facets(self):
        self.assertEqual(len(self.polyedr.facets), 4)

    def test_num_edges(self):
        self.assertEqual(len(self.polyedr.edges), 16)

    def test_perimetr1(self):
        self.assertEqual(self.polyedr.shad(), 0)

    def test_perimeter2(self):
        fake_file_content = """20.0	0.0	0.0	0.0
        8	2	8
        -0.5	-0.5	0.5
        -0.5	0.5	0.5
        0.5	0.5	0.5
        0.5	-0.5	0.5
        -0.25	-0.25	-0.5
        -0.25	0.25	-0.5
        0.25	0.25	-0.5
        0.25	-0.25	-0.5
        4	1    2    3    4
        4   5    6    7    8"""
        fake_file_path = 'data/cccc.geom'
        with patch('shadow.polyedr.open'.format(__name__),
                   new=mock_open(read_data=fake_file_content)) as _file:
            self.polyedr = Polyedr(fake_file_path)
            _file.assert_called_once_with(fake_file_path)
        self.assertAlmostEqual(self.polyedr.shad(), 2)

    def test_perimeter3(self):
        fake_file_content = """100.0	0.0	0.0	0.0
8	6	24
1	1	1
1	-1	1
-1	-1	1
-1	1	1
2	2	5
2	-2	5
-2	-2	5
-2	2	5
4	1    2    3    4
4	5    6    2    1
4	3    2    6    7
4	3    7    8    4
4	1    4    8    5
4	5    6    7    8"""
        fake_file_path = 'data/cccc.geom'
        with patch('shadow.polyedr.open'.format(__name__),
                   new=mock_open(read_data=fake_file_content)) as _file:
            self.polyedr = Polyedr(fake_file_path)
            _file.assert_called_once_with(fake_file_path)
        self.assertAlmostEqual(self.polyedr.shad(), 8)

    def test_perimeter4(self):
        fake_file_content = """100.0	0.0	90.0	0.0
6	6	24
0	-1	-1
0	0	3
0	3	0
1	0	0
1	0	1
1	1	0
3	1    2    3
3	4    5    6
4	1    4    5    2
4	1    4    6    3
4	2    5    6    3"""
        fake_file_path = 'data/cccc.geom'
        with patch('shadow.polyedr.open'.format(__name__),
                   new=mock_open(read_data=fake_file_content)) as _file:
            self.polyedr = Polyedr(fake_file_path)
            _file.assert_called_once_with(fake_file_path)
        self.assertAlmostEqual(self.polyedr.shad(), 3.414, 3)
